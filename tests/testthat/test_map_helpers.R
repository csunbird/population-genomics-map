# =============================================================================
# test_map_helpers.R — Unit tests for mod_map.R helper functions
# Run with: testthat::test_file("tests/testthat/test_map_helpers.R")
# =============================================================================
library(testthat)
source("../../R/genomics.R")   # he_colour, he_legend_colours
source("../../R/utils.R")       # fmt, build_popup_html
source("../../R/mod_map.R")     # colour_for_metric, marker_radius, metric_label

# ── colour_for_metric ─────────────────────────────────────────────────────────

test_that("colour_for_metric returns one colour per input for He", {
  vals <- c(0.05, 0.12, 0.18, 0.25, 0.35, NA)
  cols <- colour_for_metric("He", vals)
  expect_length(cols, 6)
  expect_match(cols[!is.na(cols)], "^#[0-9A-Fa-f]{6}$")
})

test_that("colour_for_metric critical He is red", {
  cols <- colour_for_metric("He", c(0.05))
  expect_equal(cols, "#D32F2F")
})

test_that("colour_for_metric healthy He is blue", {
  cols <- colour_for_metric("He", c(0.40))
  expect_equal(cols, "#1565C0")
})

test_that("colour_for_metric handles all-NA vector for He", {
  cols <- colour_for_metric("He", c(NA, NA))
  expect_true(all(cols == "#AAAAAA"))
})

test_that("colour_for_metric inbreeding_f returns diverging colours", {
  # High positive F should be reddish
  col_pos <- colour_for_metric("inbreeding_f", c(0.4))
  # Negative F should be teal-ish
  col_neg <- colour_for_metric("inbreeding_f", c(-0.4))
  expect_false(col_pos == col_neg)
})

test_that("colour_for_metric Ho returns valid hex strings", {
  cols <- colour_for_metric("Ho", c(0.1, 0.3, 0.5, NA))
  expect_length(cols, 4)
  expect_match(cols[!is.na(cols) & cols != "#AAAAAA"], "^#[0-9A-Fa-f]{6}$")
})

test_that("colour_for_metric n returns valid hex strings", {
  cols <- colour_for_metric("n", c(5, 10, 20, 50))
  expect_length(cols, 4)
  expect_match(cols, "^#[0-9A-Fa-f]{6}$")
})

# ── marker_radius ─────────────────────────────────────────────────────────────

test_that("marker_radius fixed returns all 14", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  radii <- marker_radius("fixed", stats, "He")
  expect_true(all(radii == 14L))
  expect_length(radii, nrow(stats))
})

test_that("marker_radius n gives larger radius to bigger populations", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  radii <- marker_radius("n", stats, "He")
  # Largest population should have largest radius
  expect_equal(which.max(radii), which.max(stats$n))
})

test_that("marker_radius stays within [8, 24]", {
  demo  <- generate_demo_data(seed = 42)
  stats <- calc_pop_stats(demo$gt, demo$metadata)
  radii <- marker_radius("n", stats, "He")
  expect_true(all(radii >= 8 & radii <= 24))
})

test_that("marker_radius handles all-NA safely", {
  df <- data.frame(n = c(NA, NA, NA), He = c(NA, NA, NA))
  radii <- marker_radius("n", df, "n")
  expect_true(all(radii == 14L))
})

test_that("marker_radius handles all-same value safely", {
  df <- data.frame(n = c(10, 10, 10), He = c(0.2, 0.2, 0.2))
  radii <- marker_radius("n", df, "n")
  expect_true(all(radii == 14L))
})

# ── metric_label ──────────────────────────────────────────────────────────────

test_that("metric_label returns expected strings", {
  expect_equal(metric_label("He"),           "He")
  expect_equal(metric_label("Ho"),           "Ho")
  expect_equal(metric_label("pi"),           "\u03c0")
  expect_equal(metric_label("inbreeding_f"), "F")
  expect_equal(metric_label("n"),            "n")
  expect_equal(metric_label("unknown"),      "unknown")  # passthrough
})

# ── null-coalescing operator ──────────────────────────────────────────────────

test_that("%||% returns left when non-NULL and non-NA", {
  expect_equal("a" %||% "b", "a")
  expect_equal(0   %||% 1,    0)
})

test_that("%||% returns right when left is NULL", {
  expect_equal(NULL %||% "default", "default")
})

test_that("%||% returns right when left is NA", {
  expect_equal(NA %||% "fallback", "fallback")
})
