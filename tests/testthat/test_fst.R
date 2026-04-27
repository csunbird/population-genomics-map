# =============================================================================
# test_fst.R — Unit tests for Phase 2 FST functions
# Tests: hudson_fst_pair(), calc_fst_matrix(), fst_label(), fst_colour()
# =============================================================================

# Source genomics.R (and its dependencies) so these tests can run standalone
# without a running Shiny session.
if (!exists("FST_THRESHOLDS")) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(scales)
  })
  source(file.path(dirname(dirname(getwd())), "R", "genomics.R"),
         local = FALSE)
}

# =============================================================================
# hudson_fst_pair()
# =============================================================================

test_that("hudson_fst_pair: identical populations → FST = 0 (clamped)", {
  set.seed(1)
  gt <- matrix(rbinom(100, 2, 0.35), nrow = 10, ncol = 10)
  fst <- hudson_fst_pair(gt, gt)
  # Identical gt matrices: numerator per locus is negative (sampling correction
  # dominates), so sum(nums) < 0, ratio < 0, clamped to 0.
  expect_equal(fst, 0)
})

test_that("hudson_fst_pair: complete differentiation → FST ≈ 1.0", {
  # Pop A: all homozygous alt (dosage 2) → pA = 1
  # Pop B: all homozygous ref (dosage 0) → pB = 0
  # Per locus: num = (1-0)^2 - 1*(0)/(n-1) - 0*(1)/(n-1) = 1; den = 1*(1) + 0*(0) = 1
  gt_A <- matrix(rep(2L, 50), nrow = 10, ncol = 5)
  gt_B <- matrix(rep(0L, 50), nrow = 10, ncol = 5)
  fst  <- hudson_fst_pair(gt_A, gt_B)
  expect_equal(fst, 1.0, tolerance = 1e-9)
})

test_that("hudson_fst_pair: result in [0, 1] for random populations", {
  set.seed(42)
  gt_A <- matrix(rbinom(200, 2, 0.25), nrow = 10, ncol = 20)
  gt_B <- matrix(rbinom(200, 2, 0.75), nrow = 10, ncol = 20)
  fst  <- hudson_fst_pair(gt_A, gt_B)
  expect_gte(fst, 0)
  expect_lte(fst, 1)
  expect_false(is.na(fst))
})

test_that("hudson_fst_pair: returns NA for single-individual population", {
  gt_A <- matrix(c(1L, 0L, 2L), nrow = 1, ncol = 3)   # 1 individual → nA = 1 < 2
  gt_B <- matrix(rbinom(30, 2, 0.5), nrow = 10, ncol = 3)
  fst  <- hudson_fst_pair(gt_A, gt_B)
  expect_true(is.na(fst))
})

test_that("hudson_fst_pair: returns NA when all loci are missing in one pop", {
  gt_A <- matrix(NA_integer_, nrow = 5, ncol = 4)
  gt_B <- matrix(rbinom(20, 2, 0.4), nrow = 5, ncol = 4)
  fst  <- hudson_fst_pair(gt_A, gt_B)
  expect_true(is.na(fst))
})

test_that("hudson_fst_pair: handles NAs within loci gracefully", {
  set.seed(7)
  gt_A <- matrix(sample(c(0L, 1L, 2L, NA_integer_), 50, replace = TRUE),
                 nrow = 10, ncol = 5)
  gt_B <- matrix(sample(c(0L, 1L, 2L, NA_integer_), 50, replace = TRUE),
                 nrow = 10, ncol = 5)
  fst <- suppressWarnings(hudson_fst_pair(gt_A, gt_B))
  expect_true(is.na(fst) || (fst >= 0 && fst <= 1))
})

test_that("hudson_fst_pair: intermediate FST between 0 and 0.5 for mildly different pops", {
  # Mild divergence: pA ≈ 0.4, pB ≈ 0.6 → expect FST well below 1
  set.seed(3)
  gt_A <- matrix(rbinom(200, 2, 0.4), nrow = 20, ncol = 10)
  gt_B <- matrix(rbinom(200, 2, 0.6), nrow = 20, ncol = 10)
  fst  <- hudson_fst_pair(gt_A, gt_B)
  expect_gte(fst, 0)
  expect_lt(fst, 0.5)
})

# =============================================================================
# calc_fst_matrix()
# =============================================================================

test_that("calc_fst_matrix: diagonal = 0, matrix is symmetric", {
  set.seed(10)
  n_ind <- 10; n_loci <- 40
  gt <- matrix(rbinom(3 * n_ind * n_loci, 2, 0.4),
               nrow = 3 * n_ind, ncol = n_loci)
  rownames(gt) <- paste0("ind", seq_len(3 * n_ind))
  meta <- data.frame(
    sample_id  = rownames(gt),
    population = rep(c("PopA", "PopB", "PopC"), each = n_ind),
    latitude   = 0,
    longitude  = 0,
    stringsAsFactors = FALSE
  )

  mat <- calc_fst_matrix(gt, meta)

  # Correct dimensions and names
  expect_equal(dim(mat), c(3L, 3L))
  expect_equal(rownames(mat), c("PopA", "PopB", "PopC"))
  expect_equal(colnames(mat), c("PopA", "PopB", "PopC"))

  # Diagonal = 0
  expect_equal(diag(mat), c(PopA = 0, PopB = 0, PopC = 0))

  # Symmetry
  expect_equal(mat["PopA", "PopB"], mat["PopB", "PopA"])
  expect_equal(mat["PopA", "PopC"], mat["PopC", "PopA"])
  expect_equal(mat["PopB", "PopC"], mat["PopC", "PopB"])
})

test_that("calc_fst_matrix: off-diagonal values in [0, 1]", {
  set.seed(99)
  n_ind <- 8; n_loci <- 30
  gt <- matrix(rbinom(3 * n_ind * n_loci, 2, 0.5),
               nrow = 3 * n_ind, ncol = n_loci)
  rownames(gt) <- paste0("s", seq_len(3 * n_ind))
  meta <- data.frame(
    sample_id  = rownames(gt),
    population = rep(c("X", "Y", "Z"), each = n_ind),
    latitude   = 0, longitude = 0,
    stringsAsFactors = FALSE
  )
  mat <- calc_fst_matrix(gt, meta)
  off <- mat[upper.tri(mat)]
  off_valid <- off[!is.na(off)]
  expect_true(all(off_valid >= 0))
  expect_true(all(off_valid <= 1))
})

test_that("calc_fst_matrix: works with two populations", {
  set.seed(55)
  gt_a <- matrix(rbinom(80, 2, 0.3), nrow = 8, ncol = 10)
  gt_b <- matrix(rbinom(80, 2, 0.7), nrow = 8, ncol = 10)
  gt   <- rbind(gt_a, gt_b)
  rownames(gt) <- paste0("i", seq_len(16))
  meta <- data.frame(
    sample_id  = rownames(gt),
    population = rep(c("A", "B"), each = 8),
    latitude   = 0, longitude = 0,
    stringsAsFactors = FALSE
  )
  mat <- calc_fst_matrix(gt, meta)
  expect_equal(dim(mat), c(2L, 2L))
  expect_gte(mat["A", "B"], 0)
  expect_lte(mat["A", "B"], 1)
})

# =============================================================================
# fst_label()
# =============================================================================

test_that("fst_label: returns correct labels at threshold boundaries", {
  expect_equal(fst_label(NA_real_),   "Unknown")
  expect_equal(fst_label(0),          "Little differentiation")
  expect_equal(fst_label(0.04),       "Little differentiation")
  expect_equal(fst_label(0.05),       "Moderate differentiation")
  expect_equal(fst_label(0.14),       "Moderate differentiation")
  expect_equal(fst_label(0.15),       "Great differentiation")
  expect_equal(fst_label(0.24),       "Great differentiation")
  expect_equal(fst_label(0.25),       "Very great differentiation")
  expect_equal(fst_label(1.00),       "Very great differentiation")
})

test_that("fst_label: vectorised output has same length as input", {
  vals   <- c(0.01, 0.10, 0.20, 0.30, NA)
  labels <- fst_label(vals)
  expect_equal(length(labels), 5L)
  expect_equal(labels[5], "Unknown")
})

# =============================================================================
# fst_colour()
# =============================================================================

test_that("fst_colour: returns valid hex strings", {
  cols <- fst_colour(c(0, 0.25, 0.5, 1, NA))
  expect_equal(length(cols), 5L)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
})

test_that("fst_colour: NA maps to grey #DDDDDD", {
  col_na <- fst_colour(NA_real_)
  expect_equal(col_na, "#DDDDDD")
})

test_that("fst_colour: FST = 0 maps to near-white", {
  col_zero <- fst_colour(0)
  # Should be exactly #FFFFFF (first palette stop)
  expect_equal(col_zero, "#FFFFFF")
})

test_that("fst_colour: higher FST produces darker colour (monotone)", {
  # Convert hex → luminance proxy: sum of RGB values (lower = darker)
  hex_to_sum <- function(h) {
    r <- strtoi(substr(h, 2, 3), 16L)
    g <- strtoi(substr(h, 4, 5), 16L)
    b <- strtoi(substr(h, 6, 7), 16L)
    r + g + b
  }
  vals <- c(0, 0.25, 0.5, 0.75, 1.0)
  cols <- fst_colour(vals)
  lums <- sapply(cols, hex_to_sum)
  # Luminance should be (weakly) decreasing as FST increases
  expect_true(all(diff(lums) <= 0))
})

test_that("fst_colour: clamping — values outside [0,1] map to palette endpoints", {
  col_neg  <- fst_colour(-0.5)
  col_over <- fst_colour(1.5)
  col_zero <- fst_colour(0)
  col_one  <- fst_colour(1)
  expect_equal(col_neg,  col_zero)
  expect_equal(col_over, col_one)
})
