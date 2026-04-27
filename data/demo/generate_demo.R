# =============================================================================
# generate_demo.R
# Generates a small simulated VCF file for the popgen-map demo dataset.
#
# The demo simulates a conservation genetics survey for a fictional large felid
# across five populations in a fragmented tropical landscape (Borneo-like),
# with varying levels of genetic diversity reflecting conservation status.
#
# Run from the project root:
#   source("data/demo/generate_demo.R")
#
# Output: data/demo/demo.vcf  (~150 KB plain-text VCF, not committed to git)
# =============================================================================

set.seed(42)

# ── Population parameters ────────────────────────────────────────────────────
pops <- list(
  list(name = "Central_Highland",  n = 20, he = 0.32),
  list(name = "Northern_Forest",   n = 20, he = 0.27),
  list(name = "River_Delta",       n = 18, he = 0.21),
  list(name = "Eastern_Ridge",     n = 15, he = 0.14),
  list(name = "Southern_Peatland", n = 12, he = 0.09)
)

n_loci <- 350
total_n <- sum(sapply(pops, `[[`, "n"))

# ── Sample IDs ───────────────────────────────────────────────────────────────
sample_ids <- unlist(lapply(pops, function(p) {
  paste0(p$name, "_ind", seq_len(p$n))
}))

# ── Simulate genotypes ───────────────────────────────────────────────────────
# For each population, draw per-locus allele frequencies targeting the desired He
# He = 2*p*(1-p), so p in [0.5-w, 0.5+w] where w = sqrt(0.25 - He/2)

geno_matrix <- matrix(NA_integer_, nrow = total_n, ncol = n_loci,
                      dimnames = list(sample_ids, paste0("SNP", seq_len(n_loci))))

row_offset <- 0
for (pop in pops) {
  h <- pop$he
  w <- sqrt(max(0, 0.25 - h / 2))
  p_vec <- runif(n_loci, max(0.01, 0.5 - w), min(0.99, 0.5 + w))

  rows <- seq_len(pop$n) + row_offset
  for (locus in seq_len(n_loci)) {
    geno_matrix[rows, locus] <- rbinom(pop$n, 2, p_vec[locus])
  }
  row_offset <- row_offset + pop$n
}

# Remove monomorphic loci
is_poly <- apply(geno_matrix, 2, function(x) length(unique(x[!is.na(x)])) > 1)
geno_matrix <- geno_matrix[, is_poly, drop = FALSE]
n_loci_poly <- ncol(geno_matrix)
message("Retained ", n_loci_poly, " polymorphic loci out of ", n_loci)

# ── Convert to VCF genotype strings ─────────────────────────────────────────
dosage_to_gt <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "./.",
    x == 0L  ~ "0/0",
    x == 1L  ~ "0/1",
    x == 2L  ~ "1/1",
    TRUE     ~ "./."
  )
}

gt_strings <- apply(geno_matrix, c(1, 2), dosage_to_gt)
# gt_strings: samples × loci — need to transpose to loci × samples for VCF
gt_vcf <- t(gt_strings)  # loci × samples

# ── Write VCF ────────────────────────────────────────────────────────────────
vcf_path <- "data/demo/demo.vcf"
con <- file(vcf_path, open = "w")

# VCF header
writeLines(c(
  "##fileformat=VCFv4.2",
  paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
  "##source=popgen-map generate_demo.R (simulated data)",
  "##reference=simulated",
  '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
          sample_ids), collapse = "\t")
), con)

# VCF data rows
chrom_pos <- 1000 * seq_len(n_loci_poly)   # arbitrary positions
for (i in seq_len(n_loci_poly)) {
  # Compute allele frequency for INFO field
  dosages <- geno_matrix[, i]
  af <- round(sum(dosages, na.rm = TRUE) / (2 * sum(!is.na(dosages))), 4)

  row <- paste(c(
    "1",                      # CHROM
    chrom_pos[i],             # POS
    paste0("sim_snp_", i),   # ID
    "A",                      # REF
    "T",                      # ALT
    ".",                      # QUAL
    "PASS",                   # FILTER
    paste0("AF=", af),        # INFO
    "GT",                     # FORMAT
    gt_vcf[i, ]               # sample genotypes
  ), collapse = "\t")

  writeLines(row, con)
}

close(con)
message("Demo VCF written to: ", vcf_path)
message("Samples: ", total_n, " | Loci: ", n_loci_poly)
