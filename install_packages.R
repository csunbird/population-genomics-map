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

# Testing
install.packages("testthat")

# Package reproducibility (recommended)
install.packages("renv")

message("\n--- All packages installed ---")
message("To initialise renv for reproducible environments, run: renv::init()")
message("To launch the app, run: shiny::runApp()")
