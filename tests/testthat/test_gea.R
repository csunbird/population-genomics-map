# =============================================================================
# test_gea.R — Unit tests for Phase 4 adaptive potential functions
#
# Tests:
#   calc_pop_allele_freqs()
#   calc_rda()
#   calc_genomic_offset()
#   offset_risk_label()
#   offset_colour()
#   validate_env_data()
#
# NOTE: calc_lfmm2_gea() requires the 'lfmm' CRAN package and a large-ish
# genotype matrix to converge; it is tested with a skip_if_not_installed()
# guard so CI without 'lfmm' still passes all other tests.
# =============================================================================

# Source genomics.R so tests run standalone without a Shiny session.
if (!exists("calc_pop_allele_freqs")) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
  })
  source(file.path(dirname(dirname(getwd())), "R", "genomics.R"),
         local = FALSE)
}

# =============================================================================
# Shared fixtures
# =============================================================================

#' Minimal test data: 3 populations × 8 samples each × 200 loci.
#' Allele frequencies are not random — Pop3 has higher alt-allele frequency
#' so there is a real signal for env var BIO1 (which is highest for Pop3).
make_gea_fixture <- function(seed = 42L) {
  set.seed(seed)
  n_per_pop <- 8L; n_loci <- 200L; n_pops <- 3L
  n_samp    <- n_per_pop * n_pops

  # Generate gt with controlled allele-freq gradient across populations
  gt <- matrix(NA_integer_, nrow = n_samp, ncol = n_loci)
  probs <- c(0.20, 0.55, 0.25,  # Pop1 — moderate
             0.15, 0.50, 0.35,  # Pop2 — slightly higher alt freq
             0.10, 0.45, 0.45)  # Pop3 — highest alt freq
  prob_mat <- matrix(probs, nrow = 3, byrow = TRUE)

  for (p in seq_len(n_pops)) {
    idx <- ((p - 1L) * n_per_pop + 1L):(p * n_per_pop)
    gt[idx, ] <- matrix(
      sample(c(0L, 1L, 2L), n_per_pop * n_loci, replace = TRUE,
             prob = prob_mat[p, ]),
      nrow = n_per_pop
    )
  }
  rownames(gt) <- paste0("s", seq_len(n_samp))

  meta <- data.frame(
    sample_id  = paste0("s", seq_len(n_samp)),
    population = rep(c("PopA", "PopB", "PopC"), each = n_per_pop),
    latitude   = c(rep(5.0, n_per_pop), rep(2.0, n_per_pop), rep(-1.0, n_per_pop)),
    longitude  = c(rep(30,  n_per_pop), rep(33,  n_per_pop), rep(36,   n_per_pop)),
    stringsAsFactors = FALSE
  )

  # Environmental data: BIO1 increases with population (correlated with AF)
  env <- data.frame(
    population = c("PopA", "PopB", "PopC"),
    BIO1       = c(220L, 240L, 260L),   # strong positive gradient
    BIO12      = c(2800L, 2600L, 2400L), # negative gradient
    stringsAsFactors = FALSE
  )

  # Future env: BIO1 +15, BIO12 -150 uniformly (warming, drying)
  fut_env <- data.frame(
    population = c("PopA", "PopB", "PopC"),
    BIO1       = c(235L, 255L, 275L),
    BIO12      = c(2650L, 2450L, 2250L),
    stringsAsFactors = FALSE
  )

  list(gt = gt, meta = meta, env = env, fut_env = fut_env)
}

# =============================================================================
# calc_pop_allele_freqs()
# =============================================================================

test_that("calc_pop_allele_freqs: returns populations × loci matrix", {
  d   <- make_gea_fixture()
  af  <- calc_pop_allele_freqs(d$gt, d$meta)
  expect_true(is.matrix(af))
  expect_equal(nrow(af), 3L)
  expect_equal(ncol(af), ncol(d$gt))
  expect_equal(sort(rownames(af)), c("PopA", "PopB", "PopC"))
})

test_that("calc_pop_allele_freqs: frequencies are in [0, 1]", {
  d   <- make_gea_fixture()
  af  <- calc_pop_allele_freqs(d$gt, d$meta)
  valid <- af[!is.na(af)]
  expect_true(all(valid >= 0))
  expect_true(all(valid <= 1))
})

test_that("calc_pop_allele_freqs: returns NA for all-missing locus", {
  d <- make_gea_fixture()
  d$gt[, 5L] <- NA_integer_   # set locus 5 entirely to NA
  af <- calc_pop_allele_freqs(d$gt, d$meta)
  expect_true(all(is.na(af[, 5L])))
})

test_that("calc_pop_allele_freqs: single population returns single-row matrix", {
  d <- make_gea_fixture()
  d$meta$population <- "One"
  af <- calc_pop_allele_freqs(d$gt, d$meta)
  expect_equal(nrow(af), 1L)
  expect_equal(rownames(af), "One")
})

# =============================================================================
# calc_rda() — requires vegan
# =============================================================================

test_that("calc_rda: returns expected list structure", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  res <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  expect_type(res, "list")
  expected_names <- c("rda", "site_scores", "biplot_scores", "species_scores",
                       "populations", "env_vars", "env_scaled", "geo_coords",
                       "axis_labels", "top_loci_idx", "pct_constrained",
                       "n_pops", "n_loci")
  expect_true(all(expected_names %in% names(res)))
})

test_that("calc_rda: populations match input", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  res <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  expect_equal(sort(res$populations), c("PopA", "PopB", "PopC"))
  expect_equal(res$n_pops, 3L)
})

test_that("calc_rda: site_scores has one row per population", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  res <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  expect_equal(nrow(res$site_scores), 3L)
})

test_that("calc_rda: biplot_scores has one row per env variable", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  res <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  expect_equal(nrow(res$biplot_scores), 2L)
  expect_equal(sort(rownames(res$biplot_scores)), c("BIO1", "BIO12"))
})

test_that("calc_rda: pct_constrained is in [0, 100]", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  res <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  expect_gte(res$pct_constrained, 0)
  expect_lte(res$pct_constrained, 100)
})

test_that("calc_rda: top_loci_idx are valid column indices", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  res <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  expect_true(all(res$top_loci_idx >= 1L))
  expect_true(all(res$top_loci_idx <= res$n_loci))
})

test_that("calc_rda: env_scaled retains center/scale attributes", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  res <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"),
                  scale_env = TRUE)
  expect_false(is.null(attr(res$env_scaled, "scaled:center")))
  expect_false(is.null(attr(res$env_scaled, "scaled:scale")))
})

test_that("calc_rda: too few populations raises an error", {
  skip_if_not_installed("vegan")
  set.seed(1L)
  gt   <- matrix(sample(0:2, 16L * 100L, replace = TRUE),
                 nrow = 16L, ncol = 100L)
  rownames(gt) <- paste0("s", seq_len(16L))
  meta <- data.frame(
    sample_id  = paste0("s", seq_len(16L)),
    population = rep(c("A", "B"), each = 8L),
    latitude   = c(rep(1, 8), rep(-1, 8)),
    longitude  = c(rep(10, 8), rep(-10, 8)),
    stringsAsFactors = FALSE
  )
  env <- data.frame(population = c("A", "B"),
                    BIO1 = c(200L, 250L), stringsAsFactors = FALSE)
  # 2 pops < 3 minimum → should error
  expect_error(calc_rda(gt, meta, env, env_vars = "BIO1"), regexp = "3 populations")
})

# =============================================================================
# calc_genomic_offset()
# =============================================================================

test_that("calc_genomic_offset: returns expected columns", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  rda <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  off <- calc_genomic_offset(rda, d$fut_env)
  expected_cols <- c("population", "offset", "offset_norm", "lat", "lon",
                      "risk_label", "risk_colour")
  expect_true(all(expected_cols %in% names(off)))
})

test_that("calc_genomic_offset: one row per matched population", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  rda <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  off <- calc_genomic_offset(rda, d$fut_env)
  expect_equal(nrow(off), 3L)
  expect_equal(sort(off$population), c("PopA", "PopB", "PopC"))
})

test_that("calc_genomic_offset: offset values are non-negative", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  rda <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  off <- calc_genomic_offset(rda, d$fut_env)
  expect_true(all(off$offset >= 0, na.rm = TRUE))
})

test_that("calc_genomic_offset: normalised offset is in [0, 1]", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  rda <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  off <- calc_genomic_offset(rda, d$fut_env)
  valid <- off$offset_norm[!is.na(off$offset_norm)]
  expect_true(all(valid >= 0))
  expect_true(all(valid <= 1))
})

test_that("calc_genomic_offset: max normalised offset equals 1 when > 1 pop", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  rda <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  off <- calc_genomic_offset(rda, d$fut_env)
  if (nrow(off) > 1L) {
    expect_equal(max(off$offset_norm, na.rm = TRUE), 1)
  }
})

test_that("calc_genomic_offset: risk_colour is a valid hex string", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  rda <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  off <- calc_genomic_offset(rda, d$fut_env)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", off$risk_colour)))
})

test_that("calc_genomic_offset: error when no populations matched", {
  skip_if_not_installed("vegan")
  d   <- make_gea_fixture()
  rda <- calc_rda(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"))
  bad_fut <- d$fut_env
  bad_fut$population <- c("X", "Y", "Z")   # no match
  expect_error(calc_genomic_offset(rda, bad_fut), regexp = "No populations matched")
})

# =============================================================================
# offset_risk_label()
# =============================================================================

test_that("offset_risk_label: correct labels at thresholds", {
  expect_equal(offset_risk_label(0.80), "Severe — highest climate exposure")
  expect_equal(offset_risk_label(0.75), "Severe — highest climate exposure")   # boundary >= 0.75
  expect_equal(offset_risk_label(0.74), "High — significant climate exposure")
  expect_equal(offset_risk_label(0.50), "High — significant climate exposure")
  expect_equal(offset_risk_label(0.49), "Moderate — some climate exposure")
  expect_equal(offset_risk_label(0.25), "Moderate — some climate exposure")
  expect_equal(offset_risk_label(0.24), "Low — least climate exposure")
  expect_equal(offset_risk_label(0.00), "Low — least climate exposure")
  expect_equal(offset_risk_label(NA),   "Unknown")
})

test_that("offset_risk_label: vectorised input", {
  vals   <- c(0.1, 0.3, 0.6, 0.9, NA)
  labels <- offset_risk_label(vals)
  expect_length(labels, 5L)
  expect_equal(labels[5], "Unknown")
})

# =============================================================================
# offset_colour()
# =============================================================================

test_that("offset_colour: returns valid hex strings", {
  vals <- c(0, 0.3, 0.6, 0.9, NA)
  cols <- offset_colour(vals)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
})

test_that("offset_colour: severe = red, low = blue", {
  expect_equal(offset_colour(0.80), "#D32F2F")
  expect_equal(offset_colour(0.10), "#1565C0")
})

# =============================================================================
# validate_env_data()
# =============================================================================

make_meta_for_validate <- function() {
  data.frame(
    sample_id  = paste0("s", 1:6),
    population = rep(c("PopA", "PopB"), each = 3L),
    latitude   = rep(0, 6),
    longitude  = rep(0, 6),
    stringsAsFactors = FALSE
  )
}

test_that("validate_env_data: valid CSV passes", {
  meta <- make_meta_for_validate()
  env  <- data.frame(population = c("PopA", "PopB"),
                     BIO1 = c(200L, 250L), BIO12 = c(3000L, 2500L),
                     stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_true(val$valid)
  expect_length(val$errors, 0L)
  expect_equal(sort(val$matched_pops), c("PopA", "PopB"))
  expect_equal(sort(val$env_vars), c("BIO1", "BIO12"))
})

test_that("validate_env_data: missing 'population' column is an error", {
  meta <- make_meta_for_validate()
  env  <- data.frame(BIO1 = c(200L, 250L), stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_false(val$valid)
  expect_gt(length(val$errors), 0L)
})

test_that("validate_env_data: no env variable columns is an error", {
  meta <- make_meta_for_validate()
  env  <- data.frame(population = c("PopA", "PopB"), stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_false(val$valid)
})

test_that("validate_env_data: non-numeric env column is an error", {
  meta <- make_meta_for_validate()
  env  <- data.frame(population = c("PopA", "PopB"),
                     BIO1 = c("hot", "cold"),
                     stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_false(val$valid)
  expect_true(any(grepl("numeric", val$errors)))
})

test_that("validate_env_data: no population names matched is an error", {
  meta <- make_meta_for_validate()
  env  <- data.frame(population = c("X", "Y"),
                     BIO1 = c(200L, 250L), stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_false(val$valid)
  expect_true(any(grepl("matched", val$errors)))
})

test_that("validate_env_data: warns about populations only in env CSV", {
  meta <- make_meta_for_validate()
  env  <- data.frame(population = c("PopA", "PopB", "PopC"),
                     BIO1 = c(200L, 250L, 300L), stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_true(val$valid)
  expect_true(any(grepl("no genomic data", val$warnings)))
})

test_that("validate_env_data: warns about populations only in genomic data", {
  meta <- make_meta_for_validate()
  env  <- data.frame(population = c("PopA"),   # PopB missing
                     BIO1 = c(200L), stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_true(val$valid)
  expect_true(any(grepl("lack environmental data", val$warnings)))
})

test_that("validate_env_data: warns about NA values in env columns", {
  meta <- make_meta_for_validate()
  env  <- data.frame(population = c("PopA", "PopB"),
                     BIO1 = c(200L, NA_integer_),
                     stringsAsFactors = FALSE)
  val  <- validate_env_data(env, meta)
  expect_true(val$valid)
  expect_true(any(grepl("NA", val$warnings)))
})

# =============================================================================
# LFMM2 — requires 'lfmm' package (skip if not installed)
# =============================================================================

test_that("calc_lfmm2_gea: returns expected columns", {
  skip_if_not_installed("lfmm")
  d   <- make_gea_fixture()
  res <- calc_lfmm2_gea(d$gt, d$meta, d$env, env_vars = "BIO1", K = 2L)
  expected_cols <- c("locus_idx", "env_var", "pvalue", "p_adj",
                      "z_score", "gif", "significant")
  expect_true(all(expected_cols %in% names(res)))
})

test_that("calc_lfmm2_gea: one row per locus per env variable", {
  skip_if_not_installed("lfmm")
  d   <- make_gea_fixture()
  res <- calc_lfmm2_gea(d$gt, d$meta, d$env, env_vars = c("BIO1", "BIO12"), K = 2L)
  expect_equal(nrow(res), ncol(d$gt) * 2L)
  expect_equal(length(unique(res$env_var)), 2L)
})

test_that("calc_lfmm2_gea: p-values are in (0, 1]", {
  skip_if_not_installed("lfmm")
  d    <- make_gea_fixture()
  res  <- calc_lfmm2_gea(d$gt, d$meta, d$env, env_vars = "BIO1", K = 2L)
  pv   <- res$pvalue[!is.na(res$pvalue)]
  expect_true(all(pv > 0))
  expect_true(all(pv <= 1))
})

test_that("calc_lfmm2_gea: adjusted p-values are >= raw p-values (BH)", {
  skip_if_not_installed("lfmm")
  d   <- make_gea_fixture()
  res <- calc_lfmm2_gea(d$gt, d$meta, d$env, env_vars = "BIO1", K = 2L)
  valid <- !is.na(res$pvalue) & !is.na(res$p_adj)
  expect_true(all(res$p_adj[valid] >= res$pvalue[valid] - 1e-9))
})
