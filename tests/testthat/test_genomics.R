# =============================================================================
# test_genomics.R — Unit tests for R/genomics.R and R/utils.R
# Run with: testthat::test_file("tests/testthat/test_genomics.R")
#       or: devtools::test()
# =============================================================================
library(testthat)
source("../../R/genomics.R")
source("../../R/utils.R")

# ── gt_to_numeric ─────────────────────────────────────────────────────────────

test_that("gt_to_numeric converts genotype strings correctly", {
  expect_equal(gt_to_numeric("0/0"), 0L)
  expect_equal(gt_to_numeric("0/1"), 1L)
  expect_equal(gt_to_numeric("1/0"), 1L)
  expect_equal(gt_to_numeric("1/1"), 2L)
  expect_true(is.na(gt_to_numeric("./.")))
  expect_true(is.na(gt_to_numeric(NA)))
  expect_true(is.na(gt_to_numeric(".|.")))
})

test_that("gt_to_numeric handles phased genotypes", {
  expect_equal(gt_to_numeric("0|1"), 1L)
  expect_equal(gt_to_numeric("1|1"), 2L)
  expect_equal(gt_to_numeric("0|0"), 0L)
})

# ── generate_demo_data ────────────────────────────────────────────────────────

test_that("generate_demo_data returns correct structure", {
  demo <- generate_demo_data(seed = 1)
  expect_type(demo, "list")
  expect_true(all(c("gt", "metadata", "scenario_info") %in% names(demo)))
})

test_that("demo genotype matrix has correct dimensions", {
  demo <- generate_demo_data(seed = 42)
  expect_true(nrow(demo$gt) > 0)
  expect_true(ncol(demo$gt) > 0)
  expect_true(all(demo$gt[!is.na(demo$gt)] %in% c(0L, 1L, 2L)))
})

test_that("demo metadata has required columns", {
  demo <- generate_demo_data(seed = 42)
  meta <- demo$metadata
  expect_true(all(c("sample_id", "population", "latitude", "longitude") %in% names(meta)))
  expect_equal(length(unique(meta$population)), 5)
})

test_that("demo sample IDs match genotype matrix rownames", {
  demo     <- generate_demo_data(seed = 42)
  gt_ids   <- rownames(demo$gt)
  meta_ids <- demo$metadata$sample_id
  expect_gt(length(intersect(gt_ids, meta_ids)), 0)
})

test_that("demo produces reproducible results", {
  d1 <- generate_demo_data(seed = 99)
  d2 <- generate_demo_data(seed = 99)
  expect_identical(d1$gt, d2$gt)
})

# ── calc_pop_stats ────────────────────────────────────────────────────────────

test_that("calc_pop_stats returns one row per population", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  expect_equal(nrow(stats), length(unique(demo$metadata$population)))
})

test_that("He values are in [0, 0.5]", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  valid_he <- stats$He[!is.na(stats$He)]
  expect_true(all(valid_he >= 0 & valid_he <= 0.5))
})

test_that("Ho values are in [0, 1]", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  valid_ho <- stats$Ho[!is.na(stats$Ho)]
  expect_true(all(valid_ho >= 0 & valid_ho <= 1))
})

test_that("Central Highland has higher He than Southern Peatland", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  he_ch <- stats$He[stats$population == "Central Highland"]
  he_sp <- stats$He[stats$population == "Southern Peatland"]
  expect_gt(he_ch, he_sp)
})

test_that("calc_pop_stats includes population coordinates", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  expect_true(all(c("lat", "lon") %in% names(stats)))
  expect_true(all(!is.na(stats$lat)))
})

test_that("inbreeding F is bounded [-1, 1]", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  valid_f <- stats$inbreeding_f[!is.na(stats$inbreeding_f)]
  expect_true(all(valid_f >= -1 & valid_f <= 1))
})

test_that("pi uses diploid correction (2n/(2n-1)) not haploid (n/(n-1))", {
  # For a single pop with n individuals and He = h,
  # pi should be (2n/(2n-1)) * h, not (n/(n-1)) * h
  # We verify pi >= He (since 2n/(2n-1) > 1) and pi < He * (n/(n-1))
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  valid <- stats[!is.na(stats$He) & !is.na(stats$pi) & stats$n > 1, ]
  # pi should equal (2n/(2n-1)) * He (within rounding)
  expected_pi <- round((2 * valid$n / (2 * valid$n - 1)) * valid$He, 4)
  expect_equal(valid$pi, expected_pi)
})

# ── he_risk_label ─────────────────────────────────────────────────────────────

test_that("he_risk_label returns correct labels", {
  expect_equal(he_risk_label(0.05), "Critical — very low diversity")
  expect_equal(he_risk_label(0.12), "High risk — low diversity")
  expect_equal(he_risk_label(0.18), "Moderate risk")
  expect_equal(he_risk_label(0.25), "Low risk — moderate diversity")
  expect_equal(he_risk_label(0.35), "Healthy — high diversity")
  expect_equal(he_risk_label(NA),   "Unknown")
})

test_that("he_risk_label handles boundary values", {
  expect_equal(he_risk_label(0.10), "High risk — low diversity")  # 0.10 is NOT critical
  expect_equal(he_risk_label(0.15), "Moderate risk")              # 0.15 is NOT high risk
  expect_equal(he_risk_label(0.20), "Low risk — moderate diversity")
  expect_equal(he_risk_label(0.30), "Healthy — high diversity")
})

# ── validate_gt_meta ──────────────────────────────────────────────────────────

make_test_gt <- function(sample_ids, n_loci = 10) {
  m <- matrix(
    sample(c(0L, 1L, 2L), length(sample_ids) * n_loci, replace = TRUE),
    nrow = length(sample_ids),
    ncol = n_loci,
    dimnames = list(sample_ids, paste0("SNP", seq_len(n_loci)))
  )
  m
}

make_test_meta <- function(sample_ids, pop = "Pop_A", lat = 1.5, lon = 114.5) {
  data.frame(
    sample_id  = sample_ids,
    population = rep(pop, length(sample_ids)),
    latitude   = rep(lat, length(sample_ids)),
    longitude  = rep(lon, length(sample_ids)),
    stringsAsFactors = FALSE
  )
}

test_that("validate_gt_meta succeeds when all IDs match", {
  ids <- paste0("s", 1:10)
  gt  <- make_test_gt(ids)
  meta <- make_test_meta(ids)
  result <- validate_gt_meta(gt, meta)
  expect_true(result$valid)
  expect_equal(length(result$errors), 0)
})

test_that("validate_gt_meta returns invalid when no IDs match", {
  gt   <- make_test_gt(paste0("vcf_", 1:5))
  meta <- make_test_meta(paste0("csv_", 1:5))
  result <- validate_gt_meta(gt, meta)
  expect_false(result$valid)
  expect_true(any(grepl("No sample IDs match", result$errors)))
})

test_that("validate_gt_meta warns on partial mismatch", {
  ids_gt   <- c("s1", "s2", "s3", "s4", "s5")
  ids_meta <- c("s1", "s2", "s3", "s4", "s6")  # s5 ↔ s6 mismatch
  gt   <- make_test_gt(ids_gt)
  meta <- make_test_meta(ids_meta)
  result <- validate_gt_meta(gt, meta)
  expect_true(result$valid)
  expect_gt(length(result$warnings), 0)
})

test_that("validate_gt_meta returns invalid on duplicate sample IDs", {
  ids  <- c("s1", "s1", "s2", "s3")
  gt   <- make_test_gt(c("s1", "s2", "s3"))
  meta <- make_test_meta(ids)
  result <- validate_gt_meta(gt, meta)
  expect_false(result$valid)
  expect_true(any(grepl("Duplicate sample IDs", result$errors)))
})

test_that("validate_gt_meta returns invalid on out-of-range latitude", {
  ids  <- paste0("s", 1:5)
  gt   <- make_test_gt(ids)
  meta <- make_test_meta(ids, lat = 95)
  result <- validate_gt_meta(gt, meta)
  expect_false(result$valid)
  expect_true(any(grepl("latitude values outside", result$errors)))
})

test_that("validate_gt_meta returns invalid on out-of-range longitude", {
  ids  <- paste0("s", 1:5)
  gt   <- make_test_gt(ids)
  meta <- make_test_meta(ids, lon = 200)
  result <- validate_gt_meta(gt, meta)
  expect_false(result$valid)
  expect_true(any(grepl("longitude values outside", result$errors)))
})

test_that("validate_gt_meta warns on small population size", {
  ids  <- paste0("s", 1:4)   # below default min_n_pop = 5
  gt   <- make_test_gt(ids)
  meta <- make_test_meta(ids)
  result <- validate_gt_meta(gt, meta)
  expect_true(any(grepl("fewer than 5", result$warnings)))
})

test_that("validate_gt_meta matched_ids contains only overlapping IDs", {
  ids_gt   <- paste0("s", 1:6)
  ids_meta <- paste0("s", 4:8)
  gt   <- make_test_gt(ids_gt)
  meta <- make_test_meta(ids_meta)
  result <- validate_gt_meta(gt, meta)
  expect_equal(sort(result$matched_ids), paste0("s", 4:6))
})

# ── parse_metadata_csv ────────────────────────────────────────────────────────

test_that("parse_metadata_csv reads a valid CSV", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(
    sample_id  = c("s1", "s2"),
    population = c("Pop_A", "Pop_A"),
    latitude   = c(1.5, 1.51),
    longitude  = c(114.5, 114.51)
  ), tmp, row.names = FALSE)

  meta <- parse_metadata_csv(tmp)
  expect_equal(nrow(meta), 2)
  expect_true(all(c("sample_id", "population", "latitude", "longitude") %in% names(meta)))
  unlink(tmp)
})

test_that("parse_metadata_csv accepts alternative column names", {
  df <- data.frame(ind = "s1", pop = "A", lat = 1.0, long = 114.0,
                   stringsAsFactors = FALSE)
  meta <- parse_metadata_csv(df)
  expect_true("sample_id" %in% names(meta))
  expect_true("population" %in% names(meta))
})

test_that("parse_metadata_csv errors on missing required columns", {
  expect_error(
    parse_metadata_csv(data.frame(sample_id = "s1", population = "P")),
    "missing required columns"
  )
})

test_that("parse_metadata_csv excludes rows with invalid coordinates", {
  df <- data.frame(
    sample_id  = c("s1", "s2", "s3"),
    population = c("A",  "A",  "B"),
    latitude   = c(1.5, NA, 2.0),
    longitude  = c(114.5, 114.5, 115.0),
    stringsAsFactors = FALSE
  )
  expect_warning(meta <- parse_metadata_csv(df))
  expect_equal(nrow(meta), 2)
})
