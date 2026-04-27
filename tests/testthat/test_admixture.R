# =============================================================================
# test_admixture.R — Unit tests for ADMIXTURE helper functions
#
# Tests cover:
#   parse_qmatrix()  — CSV parsing and column standardisation
#   svg_pie_uri()    — SVG pie chart data URI generation
# =============================================================================

# Source the module directly (testthat runs in the package root)
source(file.path(rprojroot::find_root(rprojroot::is_testthat), "..", "R", "mod_admixture.R"),
       local = TRUE)

# ── parse_qmatrix ─────────────────────────────────────────────────────────────

test_that("parse_qmatrix: reads standard format A (sample_id + population + K cols)", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(
    c("sample_id,population,K1,K2,K3",
      "Ind1,PopA,0.9,0.05,0.05",
      "Ind2,PopA,0.85,0.1,0.05",
      "Ind3,PopB,0.1,0.8,0.1",
      "Ind4,PopB,0.05,0.9,0.05"),
    tmp
  )
  on.exit(unlink(tmp))

  result <- parse_qmatrix(tmp)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4L)
  expect_true(all(c("sample_id", "population", "K1", "K2", "K3") %in% names(result)))
  expect_equal(result$population, c("PopA", "PopA", "PopB", "PopB"))
  expect_equal(result$sample_id,  c("Ind1", "Ind2", "Ind3", "Ind4"))
})

test_that("parse_qmatrix: handles missing population column (format B)", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(
    c("sample_id,K1,K2",
      "Ind1,0.8,0.2",
      "Ind2,0.3,0.7"),
    tmp
  )
  on.exit(unlink(tmp))

  result <- parse_qmatrix(tmp)
  expect_equal(nrow(result), 2L)
  expect_true("population" %in% names(result))
  expect_equal(unique(result$population), "Unknown")
  expect_true(all(c("K1", "K2") %in% names(result)))
})

test_that("parse_qmatrix: handles pure numeric matrix (ADMIXTURE native format)", {
  tmp <- tempfile(fileext = ".Q")
  writeLines(
    c("0.9 0.05 0.05",
      "0.85 0.1 0.05",
      "0.1 0.8 0.1"),
    tmp
  )
  on.exit(unlink(tmp))

  result <- parse_qmatrix(tmp)
  expect_equal(nrow(result), 3L)
  expect_equal(ncol(result), 5L)  # sample_id + population + K1 + K2 + K3
  expect_true(all(c("K1", "K2", "K3") %in% names(result)))
  # sample IDs auto-generated
  expect_true(all(grepl("^Ind", result$sample_id)))
})

test_that("parse_qmatrix: K columns renamed to K1, K2, ... regardless of input names", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(
    c("sample_id,population,Cluster1,Cluster2",
      "S1,Pop1,0.7,0.3",
      "S2,Pop1,0.6,0.4"),
    tmp
  )
  on.exit(unlink(tmp))

  result <- parse_qmatrix(tmp)
  expect_true("K1" %in% names(result))
  expect_true("K2" %in% names(result))
  expect_false("Cluster1" %in% names(result))
})

test_that("parse_qmatrix: raises error when no numeric K columns exist", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(
    c("sample_id,population,notes",
      "S1,Pop1,some text",
      "S2,Pop2,other text"),
    tmp
  )
  on.exit(unlink(tmp))

  expect_error(parse_qmatrix(tmp), regexp = "K columns")
})

test_that("parse_qmatrix: values are numeric in output", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(
    c("sample_id,population,K1,K2",
      "S1,Pop1,0.6,0.4",
      "S2,Pop2,0.3,0.7"),
    tmp
  )
  on.exit(unlink(tmp))

  result <- parse_qmatrix(tmp)
  expect_true(is.numeric(result$K1))
  expect_true(is.numeric(result$K2))
  expect_equal(result$K1, c(0.6, 0.3))
})

# ── svg_pie_uri ───────────────────────────────────────────────────────────────

test_that("svg_pie_uri: returns a character string data URI", {
  uri <- svg_pie_uri(c(0.6, 0.4), c("red", "blue"))
  expect_type(uri, "character")
  expect_true(startsWith(uri, "data:image/svg+xml,"))
})

test_that("svg_pie_uri: contains SVG tags", {
  uri <- svg_pie_uri(c(0.5, 0.5), c("#E41A1C", "#377EB8"))
  # URL-decode %3C → < and check for svg tag
  decoded <- gsub("%20", " ", uri)
  expect_true(grepl("svg", decoded, ignore.case = TRUE))
})

test_that("svg_pie_uri: handles single-component (full circle)", {
  uri <- svg_pie_uri(1.0, "#E41A1C")
  expect_type(uri, "character")
  expect_true(startsWith(uri, "data:image/svg+xml,"))
  # Full circle → <circle> element, not <path>
  expect_true(grepl("circle", uri))
})

test_that("svg_pie_uri: normalises proportions that don't sum to 1", {
  # Proportions sum to 2 — should still produce valid SVG (normalised internally)
  uri <- svg_pie_uri(c(1.2, 0.8), c("red", "blue"), size = 36)
  expect_type(uri, "character")
  expect_true(startsWith(uri, "data:image/svg+xml,"))
})

test_that("svg_pie_uri: returns NULL for empty or zero proportions", {
  expect_null(svg_pie_uri(numeric(0), character(0)))
  expect_null(svg_pie_uri(c(0, 0, 0), c("red", "blue", "green")))
})

test_that("svg_pie_uri: respects size parameter in SVG dimensions", {
  uri50 <- svg_pie_uri(c(0.5, 0.5), c("red", "blue"), size = 50)
  uri80 <- svg_pie_uri(c(0.5, 0.5), c("red", "blue"), size = 80)
  expect_true(grepl("width='50'", uri50))
  expect_true(grepl("width='80'", uri80))
})

test_that("svg_pie_uri: encodes # characters in colours as %23", {
  uri <- svg_pie_uri(c(0.6, 0.4), c("#E41A1C", "#377EB8"))
  # Raw '#' should not appear — must be encoded as %23
  # (After the 'data:image/svg+xml,' prefix)
  svg_part <- sub("^data:image/svg[+]xml,", "", uri)
  expect_false(grepl("#", svg_part, fixed = TRUE))
})
