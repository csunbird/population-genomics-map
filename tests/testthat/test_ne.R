# =============================================================================
# test_ne.R — Unit tests for Phase 3 Ne functions
# Tests: calc_ne_ld(), ne_risk_label(), ne_colour()
# =============================================================================

# Source genomics.R so tests run standalone without a Shiny session.
if (!exists("calc_ne_ld")) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
  })
  source(file.path(dirname(dirname(getwd())), "R", "genomics.R"),
         local = FALSE)
}

# =============================================================================
# Test fixtures
# =============================================================================

# gt is samples × loci (rows = samples, cols = loci); rownames = sample IDs —
# canonical app convention (same as parse_vcf output).
make_ne_gt_meta <- function(n_loci = 150, n_per_pop = 12, seed = 99) {
  set.seed(seed)
  n_samp <- n_per_pop * 2
  gt     <- matrix(
    sample(c(0L, 1L, 2L), n_samp * n_loci, replace = TRUE,
           prob = c(0.25, 0.50, 0.25)),
    nrow = n_samp, ncol = n_loci
  )
  rownames(gt) <- paste0("s", seq_len(n_samp))
  meta <- data.frame(
    sample_id  = paste0("s", seq_len(n_samp)),
    population = rep(c("Pop1", "Pop2"), each = n_per_pop),
    latitude   = c(rep(5,  n_per_pop), rep(-5,  n_per_pop)),
    longitude  = c(rep(30, n_per_pop), rep(-30, n_per_pop)),
    stringsAsFactors = FALSE
  )
  list(gt = gt, meta = meta)
}

# Population with very small n (below min_n threshold)
make_tiny_pop_meta <- function(n_loci = 100, n_samp = 3, seed = 7) {
  set.seed(seed)
  gt <- matrix(
    sample(c(0L, 1L, 2L), n_samp * n_loci, replace = TRUE),
    nrow = n_samp, ncol = n_loci
  )
  rownames(gt) <- paste0("t", seq_len(n_samp))
  meta <- data.frame(
    sample_id  = paste0("t", seq_len(n_samp)),
    population = rep("Tiny", n_samp),
    latitude   = rep(0, n_samp),
    longitude  = rep(0, n_samp),
    stringsAsFactors = FALSE
  )
  list(gt = gt, meta = meta)
}

# =============================================================================
# calc_ne_ld() — core tests
# =============================================================================

test_that("calc_ne_ld: returns expected columns", {
  d   <- make_ne_gt_meta()
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L)
  expect_s3_class(res, "data.frame")
  expected_cols <- c("population", "ne", "ne_lower", "ne_upper",
                     "n_samples", "n_loci_used", "n_pairs_used",
                     "mean_r2", "mean_r2_corrected")
  expect_true(all(expected_cols %in% names(res)))
})

test_that("calc_ne_ld: one row per population", {
  d   <- make_ne_gt_meta()
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L)
  expect_equal(nrow(res), 2L)
  expect_equal(sort(res$population), c("Pop1", "Pop2"))
})

test_that("calc_ne_ld: Ne values are positive or NA", {
  d   <- make_ne_gt_meta()
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L)
  ne_vals <- res$ne[!is.na(res$ne)]
  expect_true(all(ne_vals >= 1))
})

test_that("calc_ne_ld: Ne is capped at 10000", {
  d   <- make_ne_gt_meta()
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L)
  ne_vals <- res$ne[!is.na(res$ne)]
  expect_true(all(ne_vals <= 10000))
})

test_that("calc_ne_ld: CI lower <= Ne <= CI upper (when not NA)", {
  d   <- make_ne_gt_meta()
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L)
  valid <- !is.na(res$ne) & !is.na(res$ne_lower) & !is.na(res$ne_upper)
  if (any(valid)) {
    expect_true(all(res$ne_lower[valid] <= res$ne[valid] + 1e-6))
    expect_true(all(res$ne_upper[valid] >= res$ne[valid] - 1e-6))
  }
})

test_that("calc_ne_ld: populations below min_n return NA for Ne", {
  d   <- make_tiny_pop_meta()
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L, min_n = 5L)
  # n_samp = 3 < min_n = 5 → Ne should be NA
  expect_equal(nrow(res), 1L)
  expect_true(is.na(res$ne[1]))
})

test_that("calc_ne_ld: n_loci_used does not exceed max_loci", {
  d   <- make_ne_gt_meta(n_loci = 200)
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 60L)
  expect_true(all(res$n_loci_used <= 60L | is.na(res$n_loci_used)))
})

test_that("calc_ne_ld: reproducible with same seed", {
  d    <- make_ne_gt_meta(n_loci = 200)
  res1 <- calc_ne_ld(d$gt, d$meta, max_loci = 80L)
  res2 <- calc_ne_ld(d$gt, d$meta, max_loci = 80L)
  expect_equal(res1$ne, res2$ne)
})

test_that("calc_ne_ld: n_samples column matches actual population size", {
  d   <- make_ne_gt_meta(n_per_pop = 8)
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L)
  expect_true(all(res$n_samples == 8L | is.na(res$n_samples)))
})

test_that("calc_ne_ld: handles loci with all-NA values without error", {
  d <- make_ne_gt_meta(n_loci = 100)
  d$gt[, 1:10] <- NA_integer_   # set first 10 loci (columns) entirely to NA
  expect_no_error(calc_ne_ld(d$gt, d$meta, max_loci = 40L))
})

test_that("calc_ne_ld: single population without error", {
  d <- make_ne_gt_meta(n_per_pop = 10)
  d$meta$population <- "OnePop"
  res <- calc_ne_ld(d$gt, d$meta, max_loci = 50L)
  expect_equal(nrow(res), 1L)
  expect_equal(res$population, "OnePop")
})

# =============================================================================
# ne_risk_label()
# =============================================================================

test_that("ne_risk_label: correct label at IUCN thresholds", {
  # Thresholds use strict less-than (<), so the boundary value itself falls
  # into the next-higher (less severe) category.
  expect_equal(ne_risk_label(30),  "Critical — immediate extinction risk")   # 30 < 50
  expect_equal(ne_risk_label(49),  "Critical — immediate extinction risk")   # 49 < 50
  expect_equal(ne_risk_label(50),  "High risk — rapid genetic drift")        # 50 NOT < 50
  expect_equal(ne_risk_label(51),  "High risk — rapid genetic drift")
  expect_equal(ne_risk_label(99),  "High risk — rapid genetic drift")        # 99 < 100
  expect_equal(ne_risk_label(100), "Concern — long-term viability at risk")  # 100 NOT < 100
  expect_equal(ne_risk_label(101), "Concern — long-term viability at risk")
  expect_equal(ne_risk_label(499), "Concern — long-term viability at risk")  # 499 < 500
  expect_equal(ne_risk_label(500), "Viable — Ne ≥ 500")                 # 500 NOT < 500
  expect_equal(ne_risk_label(9999),"Viable — Ne ≥ 500")
  expect_equal(ne_risk_label(NA),  "Unknown")
})

test_that("ne_risk_label: vectorised input returns correct length", {
  vals   <- c(20, 75, 300, 600, NA)
  labels <- ne_risk_label(vals)
  expect_length(labels, 5L)
  expect_equal(labels[5], "Unknown")
})

# =============================================================================
# ne_colour()
# =============================================================================

test_that("ne_colour: returns valid hex colour string", {
  expect_match(ne_colour(30),   "^#[0-9A-Fa-f]{6}$")
  expect_match(ne_colour(75),   "^#[0-9A-Fa-f]{6}$")
  expect_match(ne_colour(300),  "^#[0-9A-Fa-f]{6}$")
  expect_match(ne_colour(1000), "^#[0-9A-Fa-f]{6}$")
  expect_match(ne_colour(NA),   "^#[0-9A-Fa-f]{6}$")
})

test_that("ne_colour: critical Ne uses red, viable Ne uses blue", {
  expect_equal(ne_colour(30),   "#D32F2F")   # critical
  expect_equal(ne_colour(1000), "#1565C0")   # viable
})

# =============================================================================
# NE_THRESHOLDS constant
# =============================================================================

test_that("NE_THRESHOLDS is defined with correct keys", {
  expect_true(exists("NE_THRESHOLDS"))
  expect_equal(NE_THRESHOLDS$critical,   50)
  expect_equal(NE_THRESHOLDS$vulnerable, 100)
  expect_equal(NE_THRESHOLDS$concern,    500)
})
