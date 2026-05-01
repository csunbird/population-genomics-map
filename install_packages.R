# =============================================================================
# install_packages.R — Install all dependencies for popgen-map
# Run once before launching the app: source("install_packages.R")
# =============================================================================

# Core Shiny + UI
install.packages(c(
  "shiny",
  "bslib",
  "shinyjs"       # button enable/disable, dynamic UI interactions
))

# Interactive maps
install.packages("leaflet")

# Genomics
install.packages("vcfR")   # VCF parsing (only needed for real-data uploads)

# Data manipulation + visualisation
install.packages(c(
  "dplyr",
  "scales",
  "DT",           # interactive tables
  "htmltools"
))

# Phase 2 additions
install.packages("plotly")   # interactive PCA scatter + ADMIXTURE bar (Phase 2 Goals 3–4)
# ADMIXTURE pies use URL-encoded SVG data URIs — no extra packages needed
# install.packages("tidyr")  # needed if reshaping admixture data later

# Phase 5: Export
# rmarkdown renders the HTML conservation report from report_template.Rmd
install.packages("rmarkdown")
# knitr is typically installed alongside rmarkdown; include explicitly for completeness
install.packages("knitr")
# leaflet.extras2 provides the client-side map screenshot control (camera button)
install.packages("leaflet.extras2")
# openxlsx generates the multi-sheet Excel workbook export
install.packages("openxlsx")
# pagedown converts the rendered HTML report to PDF via Chrome/Chromium
# Optional: requires Google Chrome or Chromium to be installed on the server
# If unavailable, the app falls back to HTML report automatically
install.packages("pagedown")

# Testing
install.packages("testthat")

# Package reproducibility (recommended)
install.packages("renv")

message("\n--- All packages installed ---")
message("To initialise renv for reproducible environments, run: renv::init()")
message("To launch the app, run: shiny::runApp()")
