# =============================================================================
# genomics.R — Core population genetics computation functions
# popgen-map Phase 1
# =============================================================================
# All functions operate on a standard genotype matrix:
#   rows    = individuals
#   columns = SNP loci
#   values  = 0 (hom ref), 1 (heterozygous), 2 (hom alt), NA (missing)
# =============================================================================

# ── Conservation He thresholds (Frankham et al. 2014 / IUCN) ─────────────────
# Defined here (sourced first) so he_colour() and he_risk_label() reference
# the same constants used across the whole app. global.R loads this file before
# any module, so HE_THRESHOLDS is available everywhere at runtime.
HE_THRESHOLDS <- list(
  critical   = 0.10,
  high_risk  = 0.15,
  vulnerable = 0.20,
  low_risk   = 0.30   # boundary between "lower risk" and "healthy"
)

#' Parse a VCF file and return a named list with genotype matrix + metadata
#'
#' Applies the following filters in order:
#'   1. Retain only biallelic SNPs (FILTER = PASS or ".", no indels)
#'   2. Remove sites with global missing-data rate > max_missing
#'   3. Remove sites with minor allele frequency < maf
#'   4. Remove monomorphic sites remaining after the above
#'
#' @param vcf_path    Path to a VCF file (plain text or .gz)
#' @param maf         Minor allele frequency threshold (default 0.05)
#' @param max_missing Maximum proportion of missing genotypes per site (default 0.20)
#' @return List: $gt (numeric genotype matrix, samples × loci),
#'               $n_loci, $n_samples, $sample_ids,
#'               $filter_log (named integer vector of site counts at each step)
parse_vcf <- function(vcf_path, maf = 0.05, max_missing = 0.20) {
  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("Package 'vcfR' is required to read VCF files. Run: install.packages('vcfR')")
  }

  # ── Read ─────────────────────────────────────────────────────────────────────
  vcf <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  n_raw <- nrow(vcf)

  if (n_raw == 0) stop("VCF file contains no variant records.")

  # ── Filter 1: biallelic SNPs only ────────────────────────────────────────────
  # Keep rows where REF and ALT are each a single nucleotide (A/C/G/T)
  fix <- vcfR::getFIX(vcf)
  ref <- fix[, "REF"]
  alt <- fix[, "ALT"]

  is_biallelic_snp <- nchar(ref) == 1L &
                      nchar(alt) == 1L &
                      grepl("^[ACGTacgt]$", ref) &
                      grepl("^[ACGTacgt]$", alt)

  vcf <- vcf[is_biallelic_snp, ]
  n_biallelic <- nrow(vcf)

  if (n_biallelic == 0) {
    stop(paste0(
      "No biallelic SNPs remain after filtering. ",
      "Your VCF may contain only indels or multi-allelic variants. ",
      "Please pre-filter with bcftools view -v snps -m2 -M2 before uploading."
    ))
  }

  # ── Extract genotype matrix ──────────────────────────────────────────────────
  gt_raw <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
  # gt_raw: loci × samples

  gt_num <- apply(gt_raw, c(1, 2), gt_to_numeric)
  # gt_num: loci × samples → transpose to samples × loci
  gt_num <- t(gt_num)

  n_samples <- nrow(gt_num)
  if (n_samples < 2) stop("VCF must contain genotypes for at least 2 samples.")

  # ── Filter 2: global missing-data rate per site ──────────────────────────────
  miss_rate <- apply(gt_num, 2, function(x) mean(is.na(x)))
  pass_miss  <- miss_rate <= max_missing
  gt_num     <- gt_num[, pass_miss, drop = FALSE]
  n_after_miss <- ncol(gt_num)

  if (n_after_miss == 0) {
    stop(paste0(
      "No sites remain after applying the missing-data filter (max_missing = ",
      max_missing, "). Your dataset may have excessive missing data."
    ))
  }

  # ── Filter 3: minor allele frequency ────────────────────────────────────────
  maf_vals <- apply(gt_num, 2, function(locus) {
    valid <- locus[!is.na(locus)]
    if (length(valid) == 0) return(0)
    p <- sum(valid) / (2 * length(valid))
    min(p, 1 - p)
  })
  pass_maf <- maf_vals >= maf
  gt_num   <- gt_num[, pass_maf, drop = FALSE]
  n_after_maf <- ncol(gt_num)

  if (n_after_maf == 0) {
    stop(paste0(
      "No sites remain after applying the MAF filter (maf = ", maf, "). ",
      "Try lowering the MAF threshold."
    ))
  }

  # ── Filter 4: remove any remaining monomorphic sites ────────────────────────
  is_poly <- apply(gt_num, 2, function(locus) {
    vals <- locus[!is.na(locus)]
    length(unique(vals)) > 1
  })
  gt_num   <- gt_num[, is_poly, drop = FALSE]
  n_final  <- ncol(gt_num)

  filter_log <- c(
    raw            = n_raw,
    biallelic_snps = n_biallelic,
    after_missing  = n_after_miss,
    after_maf      = n_after_maf,
    final_polymorphic = n_final
  )

  list(
    gt          = gt_num,
    n_loci      = n_final,
    n_samples   = n_samples,
    sample_ids  = rownames(gt_num),
    filter_log  = filter_log
  )
}

#' Convert a single VCF genotype string to numeric dosage
#' "0/0" → 0, "0/1" or "1/0" → 1, "1/1" → 2, "./." → NA
#' @param gt_str Character string
#' @return Integer 0/1/2 or NA
gt_to_numeric <- function(gt_str) {
  if (is.na(gt_str) || gt_str %in% c("./.", ".|.")) return(NA_integer_)

  # Handle phased (|) and unphased (/) separators
  alleles <- strsplit(gt_str, "[/|]")[[1]]
  if (any(alleles == ".")) return(NA_integer_)

  tryCatch(
    as.integer(sum(as.integer(alleles))),
    error = function(e) NA_integer_
  )
}

#' Validate that a genotype matrix and metadata table are compatible
#'
#' Checks column presence, sample-ID overlap, coordinate ranges, duplicate IDs,
#' and minimum per-population sample sizes. Returns a structured report rather
#' than stopping on non-fatal issues, so the UI can display warnings alongside
#' errors in a single pass.
#'
#' @param gt_matrix  Numeric genotype matrix (samples × loci), rownames = sample IDs
#' @param metadata   Data frame with columns: sample_id, population, latitude, longitude
#' @param min_n_pop  Minimum samples per population for reliable statistics (default 5)
#' @return List: $errors (character), $warnings (character), $info (character),
#'               $matched_ids (character), $valid (logical)
validate_gt_meta <- function(gt_matrix, metadata, min_n_pop = 5) {
  errors   <- character(0)
  warnings <- character(0)
  info     <- character(0)

  # ── Required columns ─────────────────────────────────────────────────────────
  required_cols <- c("sample_id", "population", "latitude", "longitude")
  missing_cols  <- setdiff(required_cols, names(metadata))
  if (length(missing_cols) > 0) {
    errors <- c(errors, paste0(
      "Metadata CSV is missing required columns: ",
      paste(missing_cols, collapse = ", "), ". ",
      "Required: sample_id, population, latitude, longitude."
    ))
  }

  # ── Duplicate sample IDs in metadata ─────────────────────────────────────────
  meta_ids <- as.character(metadata$sample_id)
  dupe_ids <- meta_ids[duplicated(meta_ids)]
  if (length(dupe_ids) > 0) {
    errors <- c(errors, paste0(
      "Duplicate sample IDs found in metadata CSV: ",
      paste(head(unique(dupe_ids), 5), collapse = ", "),
      if (length(unique(dupe_ids)) > 5) paste0(" … and ", length(unique(dupe_ids)) - 5, " more.") else "."
    ))
  }

  # ── Coordinate range validation ───────────────────────────────────────────────
  if (all(c("latitude", "longitude") %in% names(metadata))) {
    lat <- suppressWarnings(as.numeric(metadata$latitude))
    lon <- suppressWarnings(as.numeric(metadata$longitude))

    bad_lat <- which(!is.na(lat) & (lat < -90 | lat > 90))
    bad_lon <- which(!is.na(lon) & (lon < -180 | lon > 180))

    if (length(bad_lat) > 0) {
      errors <- c(errors, paste0(
        length(bad_lat), " row(s) have latitude values outside [-90, 90]: rows ",
        paste(head(bad_lat, 5), collapse = ", "), "."
      ))
    }
    if (length(bad_lon) > 0) {
      errors <- c(errors, paste0(
        length(bad_lon), " row(s) have longitude values outside [-180, 180]: rows ",
        paste(head(bad_lon, 5), collapse = ", "), "."
      ))
    }

  }

  # ── Sample ID matching ────────────────────────────────────────────────────────
  gt_ids    <- rownames(gt_matrix)
  in_both   <- intersect(gt_ids, meta_ids)
  vcf_only  <- setdiff(gt_ids,   meta_ids)
  meta_only <- setdiff(meta_ids, gt_ids)

  if (length(in_both) == 0) {
    errors <- c(errors, paste0(
      "No sample IDs match between the VCF and metadata CSV. ",
      "First 3 VCF IDs: [", paste(head(gt_ids, 3), collapse = ", "), "] | ",
      "First 3 CSV IDs: [", paste(head(meta_ids, 3), collapse = ", "), "]. ",
      "Check that sample names are spelled identically in both files."
    ))
  } else {
    info <- c(info, paste0(length(in_both), " samples matched between VCF and metadata."))
    if (length(vcf_only) > 0) {
      warnings <- c(warnings, paste0(
        length(vcf_only), " VCF sample(s) have no metadata row and will be excluded: ",
        paste(head(vcf_only, 5), collapse = ", "),
        if (length(vcf_only) > 5) " …" else "."
      ))
    }
    if (length(meta_only) > 0) {
      warnings <- c(warnings, paste0(
        length(meta_only), " metadata row(s) have no matching VCF genotype and will be excluded: ",
        paste(head(meta_only, 5), collapse = ", "),
        if (length(meta_only) > 5) " …" else "."
      ))
    }
  }

  # ── Minimum per-population sample size ───────────────────────────────────────
  if (length(in_both) > 0) {
    matched_meta <- metadata[metadata$sample_id %in% in_both, ]
    pop_counts   <- table(matched_meta$population)
    small_pops   <- names(pop_counts[pop_counts < min_n_pop])

    if (length(small_pops) > 0) {
      warnings <- c(warnings, paste0(
        length(small_pops), " population(s) have fewer than ", min_n_pop,
        " matched samples (statistics may be unreliable): ",
        paste(small_pops, collapse = ", "), "."
      ))
    }

    # Info: population breakdown
    pop_summary <- paste(
      paste0(names(pop_counts), " (n=", as.integer(pop_counts), ")"),
      collapse = ", "
    )
    info <- c(info, paste0("Populations: ", pop_summary, "."))
  }

  valid <- length(errors) == 0

  list(
    errors      = errors,
    warnings    = warnings,
    info        = info,
    matched_ids = in_both,
    valid       = valid
  )
}

#' Compute per-population diversity statistics from a genotype matrix
#'
#' Returns a data frame with one row per population and columns:
#'   population, n, n_loci, Ho, He, pi, inbreeding_f
#'
#' pi uses the Nei (1978) diploid unbiased estimator: (2n / (2n-1)) * He
#'
#' @param gt_matrix Numeric genotype matrix (samples × loci), rownames = sample IDs
#' @param metadata  Data frame with columns: sample_id, population, latitude, longitude
#' @return Data frame of per-population statistics
calc_pop_stats <- function(gt_matrix, metadata) {
  # Align matrix rows with metadata
  common_ids <- intersect(rownames(gt_matrix), as.character(metadata$sample_id))
  gt_matrix  <- gt_matrix[common_ids, , drop = FALSE]
  meta       <- metadata[metadata$sample_id %in% common_ids, ]

  populations <- unique(as.character(meta$population))

  stats_list <- lapply(populations, function(pop) {
    pop_ids <- as.character(meta$sample_id[meta$population == pop])
    pop_gt  <- gt_matrix[rownames(gt_matrix) %in% pop_ids, , drop = FALSE]

    # Remove loci with >50% missing data within this population
    miss_rate <- apply(pop_gt, 2, function(x) mean(is.na(x)))
    pop_gt    <- pop_gt[, miss_rate <= 0.5, drop = FALSE]

    n  <- nrow(pop_gt)
    nl <- ncol(pop_gt)

    if (n < 2 || nl == 0) {
      return(data.frame(
        population = pop, n = n, n_loci = nl,
        Ho = NA, He = NA, pi = NA, inbreeding_f = NA,
        stringsAsFactors = FALSE
      ))
    }

    # Observed heterozygosity (Ho)
    # Proportion of heterozygous genotypes per locus, averaged across loci
    ho_per_locus <- apply(pop_gt, 2, function(locus) {
      valid <- locus[!is.na(locus)]
      if (length(valid) < 2) return(NA_real_)
      sum(valid == 1) / length(valid)
    })
    Ho <- mean(ho_per_locus, na.rm = TRUE)

    # Expected heterozygosity (He) = mean(2*p*(1-p))
    # where p = alt allele frequency per locus
    he_per_locus <- apply(pop_gt, 2, function(locus) {
      valid <- locus[!is.na(locus)]
      if (length(valid) < 2) return(NA_real_)
      p <- sum(valid) / (2 * length(valid))   # alt allele frequency
      2 * p * (1 - p)
    })
    He <- mean(he_per_locus, na.rm = TRUE)

    # Nucleotide diversity (pi) — Nei (1978) unbiased estimator for diploids
    # Correct diploid correction uses haplotype count (2n), not individual count (n):
    #   pi = (2n / (2n - 1)) * He
    # Using n/(n-1) is the haploid correction and underestimates the adjustment.
    pi_val <- (2 * n / (2 * n - 1)) * He

    # Inbreeding coefficient F = (He - Ho) / He
    inbreeding_f <- if (!is.na(He) && He > 0) (He - Ho) / He else NA_real_

    data.frame(
      population   = pop,
      n            = n,
      n_loci       = nl,
      Ho           = round(Ho, 4),
      He           = round(He, 4),
      pi           = round(pi_val, 4),
      inbreeding_f = round(inbreeding_f, 4),
      stringsAsFactors = FALSE
    )
  })

  pop_stats <- do.call(rbind, stats_list)

  # Join population coordinates (centroid of sample coordinates)
  coord_summary <- meta |>
    dplyr::group_by(population) |>
    dplyr::summarise(
      lat = mean(as.numeric(latitude),  na.rm = TRUE),
      lon = mean(as.numeric(longitude), na.rm = TRUE),
      .groups = "drop"
    )

  pop_stats <- dplyr::left_join(pop_stats, coord_summary, by = "population")

  pop_stats
}

#' Generate built-in demo dataset (no VCF file required)
#'
#' Simulates a conservation genetics scenario for a fictional large felid
#' species across a fragmented tropical landscape (Borneo-like).
#' Populations have varying diversity levels reflecting conservation status.
#'
#' @param seed Integer random seed for reproducibility
#' @return List with $gt (genotype matrix), $metadata (data frame), $scenario_info (text)
generate_demo_data <- function(seed = 42) {
  set.seed(seed)

  # --- Population definitions --------------------------------------------------
  pops <- list(
    list(name = "Central Highland",  lat =  1.52, lon = 114.48, n = 20, he = 0.32, label = "Source"),
    list(name = "Northern Forest",   lat =  1.87, lon = 115.18, n = 20, he = 0.27, label = "Connected"),
    list(name = "River Delta",       lat =  0.81, lon = 113.94, n = 18, he = 0.21, label = "Partial"),
    list(name = "Eastern Ridge",     lat =  1.23, lon = 116.05, n = 15, he = 0.14, label = "Fragmented"),
    list(name = "Southern Peatland", lat =  0.31, lon = 114.22, n = 12, he = 0.09, label = "Isolated")
  )

  n_loci <- 350

  # --- Simulate genotypes per population ---------------------------------------
  # For each population, draw per-locus allele frequencies that produce
  # the target expected heterozygosity (He = 2*p*(1-p))
  # Target He h → p ~ Uniform(0.5 - w, 0.5 + w) where w = sqrt(0.25 - h/2)
  # This ensures mean(2*p*(1-p)) ≈ h

  all_gt    <- list()
  meta_rows <- list()

  for (pop in pops) {
    h <- pop$he
    w <- sqrt(max(0, 0.25 - h / 2))
    p_vec <- runif(n_loci, max(0.01, 0.5 - w), min(0.99, 0.5 + w))

    # Draw diploid genotypes: number of alt alleles (0, 1, or 2)
    geno_mat <- sapply(p_vec, function(p) rbinom(pop$n, 2, p))
    rownames(geno_mat) <- paste0(gsub(" ", "_", pop$name), "_ind", seq_len(pop$n))

    all_gt[[pop$name]] <- geno_mat

    # Add small jitter to coordinates so individual points differ slightly
    lat_jitter <- rnorm(pop$n, 0, 0.02)
    lon_jitter <- rnorm(pop$n, 0, 0.02)

    meta_rows[[pop$name]] <- data.frame(
      sample_id  = rownames(geno_mat),
      population = pop$name,
      latitude   = round(pop$lat + lat_jitter, 4),
      longitude  = round(pop$lon + lon_jitter, 4),
      stringsAsFactors = FALSE
    )
  }

  # Stack into combined matrix and metadata data frame
  gt_combined   <- do.call(rbind, all_gt)
  meta_combined <- do.call(rbind, meta_rows)
  rownames(meta_combined) <- NULL

  # Remove monomorphic loci
  is_poly <- apply(gt_combined, 2, function(x) {
    vals <- x[!is.na(x)]
    length(unique(vals)) > 1
  })
  gt_combined <- gt_combined[, is_poly, drop = FALSE]

  # --- Simulate a K=2 Q-matrix consistent with demo populations ---------------
  # Pop structure: Highland (source) and Northern are mostly K1-ancestry;
  # River Delta and Eastern Ridge show increasing K2-admixture;
  # Southern Peatland (isolated, bottlenecked) is predominantly K2.
  # This mirrors a scenario of two historical source populations with recent
  # fragmentation and partial admixture along the central corridor.
  k2_means <- c(
    "Central Highland"  = 0.88,
    "Northern Forest"   = 0.76,
    "River Delta"       = 0.52,
    "Eastern Ridge"     = 0.34,
    "Southern Peatland" = 0.14
  )

  qmat_rows <- lapply(pops, function(pop) {
    ids   <- rownames(all_gt[[pop$name]])
    n_ind <- length(ids)
    k1_mu <- k2_means[pop$name]
    # Draw individual-level K1 values from a Beta centred on the pop mean
    alpha  <- k1_mu * 8; beta_param <- (1 - k1_mu) * 8
    k1_ind <- rbeta(n_ind, shape1 = alpha, shape2 = beta_param)
    k1_ind <- pmax(0.001, pmin(0.999, k1_ind))
    data.frame(
      sample_id  = ids,
      population = pop$name,
      K1         = round(k1_ind, 4),
      K2         = round(1 - k1_ind, 4),
      stringsAsFactors = FALSE
    )
  })
  demo_qmatrix <- do.call(rbind, qmat_rows)
  rownames(demo_qmatrix) <- NULL

  list(
    gt       = gt_combined,
    metadata = meta_combined,
    qmatrix  = demo_qmatrix,
    scenario_info = paste0(
      "Demo dataset: Simulated Sunda Clouded Leopard (Neofelis diardi), ",
      "Borneo conservation landscape. ",
      nrow(gt_combined), " individuals, ",
      ncol(gt_combined), " polymorphic SNPs, ",
      length(pops), " populations. [SIMULATED — not real genotype data]"
    )
  )
}

#' Assign a conservation risk label based on expected heterozygosity
#'
#' @param he Numeric expected heterozygosity value
#' @return Character string label
he_risk_label <- function(he) {
  dplyr::case_when(
    is.na(he)                      ~ "Unknown",
    he < HE_THRESHOLDS$critical    ~ "Critical — very low diversity",
    he < HE_THRESHOLDS$high_risk   ~ "High risk — low diversity",
    he < HE_THRESHOLDS$vulnerable  ~ "Moderate risk",
    he < HE_THRESHOLDS$low_risk    ~ "Low risk — moderate diversity",
    TRUE                           ~ "Healthy — high diversity"
  )
}

#' Assign a colour to an He value for map rendering
#'
#' @param he Numeric expected heterozygosity
#' @return Hex colour string
he_colour <- function(he) {
  dplyr::case_when(
    is.na(he)                      ~ "#AAAAAA",
    he < HE_THRESHOLDS$critical    ~ "#D32F2F",   # red   — Critical
    he < HE_THRESHOLDS$high_risk   ~ "#F57C00",   # orange — High risk
    he < HE_THRESHOLDS$vulnerable  ~ "#FBC02D",   # amber  — Moderate risk
    he < HE_THRESHOLDS$low_risk    ~ "#388E3C",   # green  — Lower risk
    TRUE                           ~ "#1565C0"    # blue   — Healthy
  )
}

# =============================================================================
# Phase 2 — Population Structure: Pairwise FST
# =============================================================================

# ── Conservation FST thresholds (Wright 1978) ─────────────────────────────────
# Developed for allozymes; SNP-based FST tends to run higher.
# Use as indicative guidance, not hard cut-offs.
FST_THRESHOLDS <- list(
  little   = 0.05,   # < 0.05  = little differentiation
  moderate = 0.15,   # 0.05–0.15 = moderate differentiation
  great    = 0.25    # 0.15–0.25 = great; > 0.25 = very great
)

#' Hudson et al. (1992) pairwise FST for one population pair
#'
#' Uses the ratio-of-averages estimator: sum(numerators) / sum(denominators)
#' across all polymorphic loci. This is more robust than average-of-ratios
#' when sample sizes are small or allele frequencies are near 0/1.
#'
#' Formula per locus:
#'   num = (pA - pB)^2 - pA*(1-pA)/(nA-1) - pB*(1-pB)/(nB-1)
#'   den = pA*(1-pB) + pB*(1-pA)
#' where pA, pB = alt-allele frequencies; nA, nB = number of non-missing diploid
#' individuals at that locus.
#'
#' @param gt_A Numeric matrix: individuals × loci for population A (dosage 0/1/2)
#' @param gt_B Numeric matrix: individuals × loci for population B (dosage 0/1/2)
#' @return Numeric FST clamped to [0, 1], or NA_real_ if insufficient data
hudson_fst_pair <- function(gt_A, gt_B) {
  if (nrow(gt_A) < 2 || nrow(gt_B) < 2 || ncol(gt_A) == 0) return(NA_real_)

  result <- vapply(seq_len(ncol(gt_A)), function(l) {
    a  <- gt_A[, l]; a <- a[!is.na(a)]
    b  <- gt_B[, l]; b <- b[!is.na(b)]
    nA <- length(a); nB <- length(b)
    if (nA < 2 || nB < 2) return(c(NA_real_, NA_real_))
    pA  <- sum(a) / (2 * nA)
    pB  <- sum(b) / (2 * nB)
    num <- (pA - pB)^2 - pA * (1 - pA) / (nA - 1) - pB * (1 - pB) / (nB - 1)
    den <- pA * (1 - pB) + pB * (1 - pA)
    c(num, den)
  }, numeric(2))

  nums  <- result[1, ]
  dens  <- result[2, ]
  valid <- !is.na(nums) & !is.na(dens) & dens > 0
  if (sum(valid) == 0) return(NA_real_)
  max(0, min(1, sum(nums[valid]) / sum(dens[valid])))
}

#' Compute the full pairwise FST matrix for all population pairs
#'
#' @param gt_matrix Numeric genotype matrix (samples × loci), rownames = sample IDs
#' @param metadata  Data frame with columns: sample_id, population
#' @return Named numeric matrix (populations × populations).
#'         Diagonal = 0; upper and lower triangles are symmetric.
calc_fst_matrix <- function(gt_matrix, metadata) {
  common_ids <- intersect(rownames(gt_matrix), as.character(metadata$sample_id))
  gt_matrix  <- gt_matrix[common_ids, , drop = FALSE]
  meta       <- metadata[metadata$sample_id %in% common_ids, ]

  populations <- unique(as.character(meta$population))
  n_pop       <- length(populations)

  fst_mat <- matrix(
    NA_real_,
    nrow     = n_pop,
    ncol     = n_pop,
    dimnames = list(populations, populations)
  )
  diag(fst_mat) <- 0

  for (i in seq_len(n_pop - 1)) {
    for (j in seq(i + 1, n_pop)) {
      ids_A <- as.character(meta$sample_id[meta$population == populations[i]])
      ids_B <- as.character(meta$sample_id[meta$population == populations[j]])
      gt_A  <- gt_matrix[rownames(gt_matrix) %in% ids_A, , drop = FALSE]
      gt_B  <- gt_matrix[rownames(gt_matrix) %in% ids_B, , drop = FALSE]
      fst   <- hudson_fst_pair(gt_A, gt_B)
      fst_mat[i, j] <- fst
      fst_mat[j, i] <- fst
    }
  }

  fst_mat
}

#' Interpret an FST value as a plain-language differentiation label
#'
#' @param fst Numeric FST value (scalar or vector)
#' @return Character string label(s)
fst_label <- function(fst) {
  dplyr::case_when(
    is.na(fst)                      ~ "Unknown",
    fst < FST_THRESHOLDS$little     ~ "Little differentiation",
    fst < FST_THRESHOLDS$moderate   ~ "Moderate differentiation",
    fst < FST_THRESHOLDS$great      ~ "Great differentiation",
    TRUE                            ~ "Very great differentiation"
  )
}

# =============================================================================
# Phase 2 — Goal 3: Principal Component Analysis
# =============================================================================

#' Compute PCA on the genotype matrix
#'
#' Missing data are imputed per-locus with the column mean (standard practice
#' for genotype PCA) before calling \code{prcomp}. No scaling is applied
#' (scale. = FALSE) because all columns are on the same 0/1/2 dosage scale.
#'
#' @param gt_matrix  Numeric genotype matrix (samples × loci), rownames = sample IDs
#' @param metadata   Data frame with columns: sample_id, population
#' @param n_pcs      Number of PCs to retain (default 10; actual n may be less)
#' @return List: $scores (data frame with PC1…PCn, sample_id, population columns),
#'               $var_explained (numeric vector, proportion variance per PC),
#'               $n_pcs (integer, actual number of PCs computed)
calc_pca <- function(gt_matrix, metadata, n_pcs = 10L) {
  common_ids <- intersect(rownames(gt_matrix), as.character(metadata$sample_id))
  gt_sub     <- gt_matrix[common_ids, , drop = FALSE]
  meta_sub   <- metadata[metadata$sample_id %in% common_ids, ]

  if (nrow(gt_sub) < 3) {
    stop("PCA requires at least 3 samples. Not enough matched samples found.")
  }
  if (ncol(gt_sub) < 2) {
    stop("PCA requires at least 2 SNP loci.")
  }

  # Mean-impute missing values per locus
  gt_imp <- apply(gt_sub, 2, function(col) {
    mu <- mean(col, na.rm = TRUE)
    col[is.na(col)] <- if (is.na(mu)) 0 else mu
    col
  })

  # Remove invariant loci that slipped through (zero variance → prcomp would warn)
  col_var <- apply(gt_imp, 2, var)
  gt_imp  <- gt_imp[, col_var > 0, drop = FALSE]

  if (ncol(gt_imp) < 2) {
    stop("Too few polymorphic loci remain after imputation for PCA.")
  }

  # Run PCA (centred, not scaled — dosage values are already on a common scale)
  pc_result <- prcomp(gt_imp, center = TRUE, scale. = FALSE)

  n_pcs_actual <- min(as.integer(n_pcs), ncol(pc_result$x))
  scores_mat   <- as.data.frame(pc_result$x[, seq_len(n_pcs_actual), drop = FALSE])
  scores_mat$sample_id <- rownames(scores_mat)

  # Join population labels
  scores <- merge(
    scores_mat,
    meta_sub[, c("sample_id", "population")],
    by       = "sample_id",
    all.x    = TRUE,
    sort     = FALSE
  )

  var_explained <- (pc_result$sdev^2 / sum(pc_result$sdev^2))[seq_len(n_pcs_actual)]

  list(
    scores        = scores,
    var_explained = var_explained,
    n_pcs         = n_pcs_actual
  )
}

# =============================================================================
# Phase 3 — Goal 1: F-ROH (Runs of Homozygosity)
# =============================================================================

#' Compute per-individual F-ROH from a genotype matrix
#'
#' Uses a sliding-window approach to identify runs of homozygosity (ROH)
#' within each individual's genotype sequence. F-ROH is calculated as:
#'   F-ROH = (number of homozygous windows above threshold) / total informative windows
#'
#' This is a SNP-based proxy for the classical ROH statistic. Without physical
#' map positions the method counts consecutive homozygous SNPs; results are
#' comparable across populations within the same dataset but not across
#' datasets with different SNP densities.
#'
#' @param gt_matrix Numeric genotype matrix (samples × loci), rownames = sample IDs.
#'                  Values: 0 (hom ref), 1 (het), 2 (hom alt), NA (missing).
#'                  Same orientation as all other genomic functions: rows = samples, columns = loci.
#' @param metadata  Data frame with columns: sample_id, population
#' @param min_snps  Minimum number of consecutive homozygous SNPs to count as a ROH
#'                  (default 50, roughly equivalent to ~500 kb at 1 SNP/10 kb density)
#' @return Data frame with columns:
#'   sample_id, population, froh, n_roh, total_roh_snps, n_informative_snps
calc_froh <- function(gt_matrix, metadata, min_snps = 50L) {
  # gt_matrix is samples × loci (rows = samples, cols = loci) — canonical convention
  # throughout the app: rownames = sample IDs, colnames = locus identifiers.

  # Guard: duplicate sample IDs would cause gt_sub[, id] to return a matrix
  # instead of a vector, breaking all downstream per-individual operations.
  dup_ids <- as.character(metadata$sample_id[duplicated(as.character(metadata$sample_id))])
  if (length(dup_ids) > 0) {
    stop(paste0(
      "Duplicate sample IDs found in metadata: ",
      paste(head(unique(dup_ids), 5), collapse = ", "),
      ". Each sample must appear exactly once."
    ))
  }

  common_ids <- intersect(rownames(gt_matrix), as.character(metadata$sample_id))
  gt_sub     <- gt_matrix[common_ids, , drop = FALSE]
  meta_sub   <- metadata[metadata$sample_id %in% common_ids, ]

  if (nrow(gt_sub) < 1) stop("No samples matched between genotype matrix and metadata.")

  results <- lapply(common_ids, function(id) {
    geno <- gt_sub[id, ]                         # one individual's dosages (row = one sample)
    # Informative loci: not missing, not heterozygous (already resolved)
    is_hom <- !is.na(geno) & (geno == 0L | geno == 2L)
    is_het <- !is.na(geno) & geno == 1L
    n_informative <- sum(is_hom | is_het, na.rm = TRUE)

    if (n_informative < min_snps) {
      pop <- as.character(meta_sub$population[meta_sub$sample_id == id][1])
      return(data.frame(
        sample_id         = id,
        population        = pop,
        froh              = NA_real_,
        n_roh             = 0L,
        total_roh_snps    = 0L,
        n_informative_snps = n_informative,
        stringsAsFactors  = FALSE
      ))
    }

    # Scan for runs of homozygosity using run-length encoding on non-missing genotypes
    # Strategy: recode as TRUE (homozygous) / FALSE (het or missing) then find runs
    hom_vec   <- is_hom   # logical vector over all loci (including NAs)

    # Walk the genome: accumulate ROH length, reset on heterozygous call.
    # Missing sites do NOT break a run (follows PLINK convention), but they
    # are NOT counted in the run length — only called homozygous SNPs are.
    # This keeps the F-ROH denominator and numerator on the same basis
    # (both count only non-missing sites) and prevents F-ROH > 1.
    n_roh          <- 0L
    total_roh_snps <- 0L
    run_len        <- 0L      # counts only non-missing homozygous SNPs in run
    run_missing    <- 0L      # missing sites bridged by the current run (not counted)

    for (k in seq_along(hom_vec)) {
      if (is.na(geno[k])) {
        # Missing site: bridge the run without counting toward its length
        run_missing <- run_missing + 1L
        next
      } else if (hom_vec[k]) {
        run_len <- run_len + 1L
      } else {
        # Heterozygous site: close current run
        if (run_len >= min_snps) {
          n_roh          <- n_roh + 1L
          total_roh_snps <- total_roh_snps + run_len  # only homozygous SNPs
        }
        run_len    <- 0L
        run_missing <- 0L
      }
    }
    # Close any run that reaches the end of the sequence
    if (run_len >= min_snps) {
      n_roh          <- n_roh + 1L
      total_roh_snps <- total_roh_snps + run_len
    }

    froh <- if (n_informative > 0) total_roh_snps / n_informative else NA_real_
    pop  <- as.character(meta_sub$population[meta_sub$sample_id == id][1])

    data.frame(
      sample_id          = id,
      population         = pop,
      froh               = round(froh, 4),
      n_roh              = n_roh,
      total_roh_snps     = total_roh_snps,
      n_informative_snps = n_informative,
      stringsAsFactors   = FALSE
    )
  })

  do.call(rbind, results)
}

#' Summarise per-individual F-ROH into per-population statistics
#'
#' @param froh_df Data frame returned by calc_froh()
#' @return Data frame with columns: population, mean_froh, sd_froh, max_froh,
#'   median_froh, n_individuals, n_with_roh
summarise_froh <- function(froh_df) {
  pops <- unique(froh_df$population)
  rows <- lapply(pops, function(pop) {
    sub   <- froh_df[froh_df$population == pop, ]
    valid <- sub$froh[!is.na(sub$froh)]
    data.frame(
      population    = pop,
      mean_froh     = if (length(valid) > 0) round(mean(valid),   4) else NA_real_,
      sd_froh       = if (length(valid) > 1) round(sd(valid),     4) else NA_real_,
      median_froh   = if (length(valid) > 0) round(median(valid), 4) else NA_real_,
      max_froh      = if (length(valid) > 0) round(max(valid),    4) else NA_real_,
      n_individuals = nrow(sub),
      n_with_roh    = sum(sub$n_roh > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Interpret an F-ROH value as a conservation risk label
#'
#' Thresholds based on Frankham (2015) and McQuillan et al. (2008).
#' Pedigree-equivalent inbreeding coefficients (diploid, random mating):
#'   F = 0.25  — full-sibling or parent-offspring mating
#'   F = 0.125 — half-sibling mating
#'   F = 0.0625 — first-cousin mating
#'
#' @param froh Numeric F-ROH value (scalar or vector)
#' @return Character risk label
froh_risk_label <- function(froh) {
  dplyr::case_when(
    is.na(froh)    ~ "Unknown",
    froh >= 0.25   ~ "Severe — equivalent to full-sib or parent-offspring",
    froh >= 0.125  ~ "High — equivalent to half-sibling mating",
    froh >= 0.0625 ~ "Moderate — equivalent to first-cousin mating",
    froh > 0       ~ "Low — background homozygosity",
    TRUE           ~ "None detected"
  )
}

# =============================================================================
# Phase 3 — Goal 2: Effective Population Size (LD-based Ne)
# =============================================================================

# ── Conservation Ne thresholds (Franklin 1980 / IUCN SSC 2023) ───────────────
NE_THRESHOLDS <- list(
  critical    = 50,    # < 50   = immediate extinction risk (50/500 rule lower bound)
  vulnerable  = 100,   # < 100  = high risk — 1% drift per generation
  concern     = 500    # < 500  = long-term viability concern (50/500 rule upper bound)
)

#' Compute LD-based effective population size (Ne) per population
#'
#' Implements the single-sample LD method of Waples & Do (2008), which
#' estimates Ne from the mean r² across all pairs of loci within a population.
#' The bias-corrected formula accounts for finite sample size.
#'
#' Formula (Waples 2006 / Waples & Do 2008):
#'   r²_corrected = r²_observed - 1/(n_samples)
#'   Ne = 1 / (3 * r²_corrected)    [for unlinked loci, Waples 2006 eq 6b]
#'
#' Confidence interval: jackknife over loci pairs.
#'
#' Caveats:
#'  - Assumes unlinked loci (no physical map used — conservative).
#'  - Assumes random mating within the population (HWE).
#'  - Very small samples (n < 10) produce unreliable estimates.
#'  - Ne is a point estimate; the jackknife CI is approximate.
#'
#' @param gt_matrix Numeric genotype matrix (samples × loci), rownames = sample IDs.
#'                  Same orientation as all other genomic functions: rows = samples, columns = loci.
#' @param metadata  Data frame with columns: sample_id, population
#' @param max_loci  Maximum loci to use per population (random subsample for speed;
#'                  default 300 — gives ~45 000 pairs, sufficient for robust estimate)
#' @param min_n     Minimum samples per population to attempt Ne estimate (default 5)
#' @return Data frame with columns:
#'   population, ne, ne_lower, ne_upper, n_samples, n_loci_used,
#'   n_pairs_used, mean_r2, mean_r2_corrected
calc_ne_ld <- function(gt_matrix, metadata, max_loci = 300L, min_n = 5L) {
  # gt_matrix is samples × loci (rows = samples, cols = loci) — canonical convention.
  # Sample IDs are in rownames.
  common_ids <- intersect(rownames(gt_matrix), as.character(metadata$sample_id))
  gt_sub     <- gt_matrix[common_ids, , drop = FALSE]
  meta_sub   <- metadata[metadata$sample_id %in% common_ids, ]

  populations <- unique(as.character(meta_sub$population))

  # Fix the random seed before any subsampling so Ne estimates are reproducible
  # for the same input data. The seed is reset at the end of each call to avoid
  # affecting the caller's RNG state.
  # Guard: .Random.seed does not exist in a fresh R session until the first RNG
  # call is made, so accessing it directly would throw "object not found".
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(
      assign(".Random.seed", old_seed, envir = globalenv()),
      add = TRUE
    )
  }
  set.seed(42L)

  results <- lapply(populations, function(pop) {
    pop_ids <- as.character(meta_sub$sample_id[meta_sub$population == pop])
    pop_gt  <- gt_sub[rownames(gt_sub) %in% pop_ids, , drop = FALSE]
    n_samp  <- nrow(pop_gt)

    if (n_samp < min_n) {
      return(data.frame(
        population      = pop,
        ne              = NA_real_,
        ne_lower        = NA_real_,
        ne_upper        = NA_real_,
        n_samples       = n_samp,
        n_loci_used     = 0L,
        n_pairs_used    = 0L,
        mean_r2         = NA_real_,
        mean_r2_corrected = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    # Keep only polymorphic loci with < 50% missing data within this population.
    # Apply over columns (margin 2) because columns = loci.
    miss_rate <- apply(pop_gt, 2, function(x) mean(is.na(x)))
    poly_mask <- apply(pop_gt, 2, function(x) {
      v <- x[!is.na(x)]; length(unique(v)) > 1
    })
    pop_gt <- pop_gt[, miss_rate <= 0.5 & poly_mask, drop = FALSE]
    n_loci <- ncol(pop_gt)

    if (n_loci < 10L) {
      return(data.frame(
        population      = pop,
        ne              = NA_real_,
        ne_lower        = NA_real_,
        ne_upper        = NA_real_,
        n_samples       = n_samp,
        n_loci_used     = n_loci,
        n_pairs_used    = 0L,
        mean_r2         = NA_real_,
        mean_r2_corrected = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    # Subsample loci for speed (seed fixed above at function entry)
    if (n_loci > max_loci) {
      keep   <- sort(sample.int(n_loci, max_loci))
      pop_gt <- pop_gt[, keep, drop = FALSE]
      n_loci <- max_loci
    }

    # Compute r² for all pairs of loci (upper triangle).
    # r² = (correlation of dosage vectors)² — valid for biallelic SNPs.
    # Mean-impute missing values before computing correlation.
    # apply(pop_gt, 2, fn) applies fn to each column (locus); fn receives a vector
    # of n_samp dosages and returns a mean-imputed vector of n_samp.
    # R's apply() transposes the result, so gt_imp is n_samp × n_loci —
    # cor(gt_imp) then gives the n_loci × n_loci correlation matrix we need.
    gt_imp <- apply(pop_gt, 2, function(col) {
      mu <- mean(col, na.rm = TRUE)
      col[is.na(col)] <- if (is.na(mu)) 0 else mu
      col
    })

    # Fast r² via cross-product of centred, scaled columns
    cor_mat    <- cor(gt_imp)                     # loci × loci correlation
    upper_idx  <- which(upper.tri(cor_mat))
    r2_vals    <- cor_mat[upper_idx]^2
    n_pairs    <- length(r2_vals)

    if (n_pairs == 0L) {
      return(data.frame(
        population      = pop,
        ne              = NA_real_,
        ne_lower        = NA_real_,
        ne_upper        = NA_real_,
        n_samples       = n_samp,
        n_loci_used     = n_loci,
        n_pairs_used    = 0L,
        mean_r2         = NA_real_,
        mean_r2_corrected = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    mean_r2   <- mean(r2_vals, na.rm = TRUE)
    # Waples (2006) bias correction for finite sample size
    r2_corr   <- mean_r2 - 1 / n_samp
    # Ne estimate — clamped to [1, 10000]
    ne_est    <- if (!is.na(r2_corr) && r2_corr > 0) {
      min(10000, max(1, 1 / (3 * r2_corr)))
    } else NA_real_

    # Jackknife CI over loci pairs (delete-one locus at a time)
    jack_ne <- vapply(seq_len(n_loci), function(drop_i) {
      keep_cols  <- seq_len(n_loci)[-drop_i]
      sub_cor    <- cor_mat[keep_cols, keep_cols]
      sub_upper  <- sub_cor[upper.tri(sub_cor)]
      sub_r2     <- mean(sub_upper^2, na.rm = TRUE)
      sub_corr   <- sub_r2 - 1 / n_samp
      if (!is.na(sub_corr) && sub_corr > 0) min(10000, max(1, 1 / (3 * sub_corr)))
      else NA_real_
    }, numeric(1))

    jack_ne  <- jack_ne[!is.na(jack_ne)]
    ne_lower <- NA_real_
    ne_upper <- NA_real_
    ne_int   <- round(ne_est)    # round here; CI is built around the same rounded value
    if (length(jack_ne) >= 3) {
      jack_se  <- sqrt((n_loci - 1) / n_loci * sum((jack_ne - mean(jack_ne))^2))
      half_w   <- round(1.96 * jack_se)     # symmetric half-width in integer Ne units
      ne_lower <- max(1,     ne_int - half_w)
      ne_upper <- min(10000, ne_int + half_w)
    }

    data.frame(
      population        = pop,
      ne                = ne_int,           # integer — LD method does not support sub-unit precision
      ne_lower          = ne_lower,
      ne_upper          = ne_upper,
      n_samples         = n_samp,
      n_loci_used       = n_loci,
      n_pairs_used      = n_pairs,
      mean_r2           = round(mean_r2, 5),
      mean_r2_corrected = round(r2_corr,  5),
      stringsAsFactors  = FALSE
    )
  })

  do.call(rbind, results)
}

#' Interpret an Ne estimate as a conservation risk label
#'
#' @param ne Numeric effective population size (scalar or vector)
#' @return Character risk label
ne_risk_label <- function(ne) {
  dplyr::case_when(
    is.na(ne)                       ~ "Unknown",
    ne < NE_THRESHOLDS$critical     ~ "Critical — immediate extinction risk",
    ne < NE_THRESHOLDS$vulnerable   ~ "High risk — rapid genetic drift",
    ne < NE_THRESHOLDS$concern      ~ "Concern — long-term viability at risk",
    TRUE                            ~ "Viable — Ne ≥ 500"
  )
}

#' Assign a colour to an Ne value for map rendering
#'
#' @param ne Numeric effective population size
#' @return Hex colour string
ne_colour <- function(ne) {
  dplyr::case_when(
    is.na(ne)                       ~ "#AAAAAA",
    ne < NE_THRESHOLDS$critical     ~ "#D32F2F",   # red    — Critical
    ne < NE_THRESHOLDS$vulnerable   ~ "#F57C00",   # orange — High risk
    ne < NE_THRESHOLDS$concern      ~ "#FBC02D",   # amber  — Concern
    TRUE                            ~ "#1565C0"    # blue   — Viable
  )
}

# =============================================================================
# (Existing) FST colour helper — retained below
# =============================================================================

#' Map FST values to hex colours for heatmap rendering
#'
#' Sequential palette: white (FST = 0, no differentiation) through teal to
#' dark teal (FST = 1, complete differentiation). Accepts a vector for
#' efficient batch colouring of entire matrices.
#'
#' @param fst_vec Numeric vector of FST values (may contain NA)
#' @return Character vector of hex colour strings, same length as fst_vec
fst_colour <- function(fst_vec) {
  pal <- scales::col_numeric(
    palette  = c("#FFFFFF", "#E0F2F1", "#26A69A", "#004D40"),
    domain   = c(0, 1),
    na.color = "#DDDDDD"
  )
  clamped        <- pmax(0, pmin(1, fst_vec))
  cols           <- pal(clamped)
  cols[is.na(fst_vec)] <- "#DDDDDD"
  cols
}
