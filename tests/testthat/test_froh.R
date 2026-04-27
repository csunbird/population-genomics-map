# =============================================================================
# test_froh.R — Unit tests for Phase 3 F-ROH functions
# Tests: calc_froh(), summarise_froh(), froh_risk_label()
# =============================================================================

# Source genomics.R so tests run standalone without a Shiny session.
if (!exists("calc_froh")) {
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

# Build a tiny gt matrix and metadata for two populations (A, B)
# gt is samples × loci (rows = samples, cols = loci); rownames = sample IDs;
# values 0/1/2 (dosage), NA = missing — canonical app convention.
make_gt_meta <- function(n_loci = 100, n_per_pop = 5, seed = 42) {
  set.seed(seed)
  n_samp <- n_per_pop * 2
  gt     <- matrix(
    sample(c(0L, 1L, 2L, NA_integer_), n_samp * n_loci, replace = TRUE,
           prob = c(0.25, 0.50, 0.25, 0.0)),
    nrow = n_samp, ncol = n_loci
  )
  rownames(gt) <- paste0("s", seq_len(n_samp))
  meta <- data.frame(
    sample_id  = paste0("s", seq_len(n_samp)),
    population = rep(c("A", "B"), each = n_per_pop),
    latitude   = c(rep(10, n_per_pop), rep(-10, n_per_pop)),
    longitude  = c(rep(20, n_per_pop), rep(-20, n_per_pop)),
    stringsAsFactors = FALSE
  )
  list(gt = gt, meta = meta)
}

# Fully homozygous individuals: all SNPs 0 or 2 (never 1)
make_fully_hom_gt <- function(n_loci = 200, n_samp = 3, seed = 7) {
  set.seed(seed)
  gt <- matrix(
    sample(c(0L, 2L), n_samp * n_loci, replace = TRUE),
    nrow = n_samp, ncol = n_loci
  )
  rownames(gt) <- paste0("h", seq_len(n_samp))
  meta <- data.frame(
    sample_id  = paste0("h", seq_len(n_samp)),
    population = rep("Hom", n_samp),
    latitude   = rep(0, n_samp),
    longitude  = rep(0, n_samp),
    stringsAsFactors = FALSE
  )
  list(gt = gt, meta = meta)
}

# Fully heterozygous individuals: all SNPs 1
make_fully_het_gt <- function(n_loci = 200, n_samp = 3) {
  gt <- matrix(1L, nrow = n_samp, ncol = n_loci)
  rownames(gt) <- paste0("e", seq_len(n_samp))
  meta <- data.frame(
    sample_id  = paste0("e", seq_len(n_samp)),
    population = rep("Het", n_samp),
    latitude   = rep(0, n_samp),
    longitude  = rep(0, n_samp),
    stringsAsFactors = FALSE
  )
  list(gt = gt, meta = meta)
}

# =============================================================================
# calc_froh() — core tests
# =============================================================================

test_that("calc_froh: returns expected columns", {
  d   <- make_gt_meta()
  res <- calc_froh(d$gt, d$meta, min_snps = 10L)
  expect_s3_class(res, "data.frame")
  expected_cols <- c("sample_id", "population", "froh",
                     "n_roh", "total_roh_snps", "n_informative_snps")
  expect_true(all(expected_cols %in% names(res)))
})

test_that("calc_froh: one row per sample", {
  d   <- make_gt_meta(n_per_pop = 4)
  res <- calc_froh(d$gt, d$meta, min_snps = 10L)
  expect_equal(nrow(res), nrow(d$meta))
})

test_that("calc_froh: F-ROH in [0, 1] for all individuals", {
  d   <- make_gt_meta()
  res <- calc_froh(d$gt, d$meta, min_snps = 10L)
  expect_true(all(res$froh >= 0 & res$froh <= 1, na.rm = TRUE))
})

test_that("calc_froh: fully heterozygous individual has F-ROH = 0", {
  d   <- make_fully_het_gt(n_loci = 200, n_samp = 2)
  res <- calc_froh(d$gt, d$meta, min_snps = 10L)
  expect_true(all(res$froh == 0))
})

test_that("calc_froh: fully homozygous individual has F-ROH > 0 (with short min_snps)", {
  d   <- make_fully_hom_gt(n_loci = 200, n_samp = 2)
  res <- calc_froh(d$gt, d$meta, min_snps = 5L)
  expect_true(all(res$froh > 0))
})

test_that("calc_froh: min_snps = 1 produces same or more ROH than min_snps = 100", {
  d     <- make_gt_meta(n_loci = 200)
  res1  <- calc_froh(d$gt, d$meta, min_snps = 1L)
  res100 <- calc_froh(d$gt, d$meta, min_snps = 100L)
  # Shorter min means shorter runs qualify: total_roh_snps can only be >= with short min
  expect_true(all(res1$froh >= res100$froh - 1e-9))
})

test_that("calc_froh: single population works (no subsetting error)", {
  d <- make_gt_meta(n_per_pop = 3)
  d$meta$population <- "SinglePop"
  res <- calc_froh(d$gt, d$meta, min_snps = 10L)
  expect_equal(nrow(res), nrow(d$meta))
  expect_true(all(res$population == "SinglePop"))
})

test_that("calc_froh: handles loci with all-missing values without error", {
  d <- make_gt_meta(n_loci = 60)
  # Set first 5 loci entirely to NA
  d$gt[1:5, ] <- NA_integer_
  expect_no_error(calc_froh(d$gt, d$meta, min_snps = 5L))
})

test_that("calc_froh: population ordering matches metadata", {
  d   <- make_gt_meta(n_per_pop = 3)
  res <- calc_froh(d$gt, d$meta, min_snps = 5L)
  expect_equal(sort(unique(res$population)), c("A", "B"))
})

test_that("calc_froh: throws informative error on duplicate sample IDs", {
  d <- make_gt_meta(n_per_pop = 3)
  # Introduce a duplicate ID
  d$meta$sample_id[2] <- d$meta$sample_id[1]
  expect_error(
    calc_froh(d$gt, d$meta, min_snps = 5L),
    regexp = "Duplicate sample IDs"
  )
})

# =============================================================================
# summarise_froh()
# =============================================================================

test_that("summarise_froh: returns one row per population", {
  d    <- make_gt_meta(n_per_pop = 4)
  froh <- calc_froh(d$gt, d$meta, min_snps = 5L)
  summ <- summarise_froh(froh)
  expect_equal(nrow(summ), 2L)
  expect_true(all(c("population", "mean_froh", "max_froh",
                    "n_individuals", "n_with_roh") %in% names(summ)))
})

test_that("summarise_froh: mean_froh in [0, 1]", {
  d    <- make_gt_meta(n_per_pop = 5)
  froh <- calc_froh(d$gt, d$meta, min_snps = 5L)
  summ <- summarise_froh(froh)
  expect_true(all(summ$mean_froh >= 0 & summ$mean_froh <= 1))
})

test_that("summarise_froh: n_individuals equals population sample size", {
  d    <- make_gt_meta(n_per_pop = 6)
  froh <- calc_froh(d$gt, d$meta, min_snps = 5L)
  summ <- summarise_froh(froh)
  expect_true(all(summ$n_individuals == 6L))
})

# =============================================================================
# froh_risk_label()
# =============================================================================

test_that("froh_risk_label: correct labels at boundary values", {
  # Pedigree-equivalent thresholds (diploid):
  # F = 0.25 → full-sib or parent-offspring (NOT half-sib)
  # F = 0.125 → half-sibling (NOT first-cousin)
  # F = 0.0625 → first-cousin (NOT second-cousin)
  expect_equal(froh_risk_label(0.30),  "Severe — equivalent to full-sib or parent-offspring")
  expect_equal(froh_risk_label(0.25),  "Severe — equivalent to full-sib or parent-offspring")
  expect_equal(froh_risk_label(0.13),  "High — equivalent to half-sibling mating")
  expect_equal(froh_risk_label(0.125), "High — equivalent to half-sibling mating")
  expect_equal(froh_risk_label(0.07),  "Moderate — equivalent to first-cousin mating")
  expect_equal(froh_risk_label(0.0625),"Moderate — equivalent to first-cousin mating")
  expect_equal(froh_risk_label(0.01),  "Low — background homozygosity")
  expect_equal(froh_risk_label(0),     "None detected")
  expect_equal(froh_risk_label(NA),    "Unknown")
})

test_that("froh_risk_label: vectorised over a vector input", {
  vals   <- c(0, 0.03, 0.07, 0.13, 0.30, NA)
  labels <- froh_risk_label(vals)
  expect_length(labels, 6L)
  expect_equal(labels[6], "Unknown")
})
