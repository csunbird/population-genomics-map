
# =============================================================================
# global.R — Package loading, global constants, module sourcing
# Runs once when the app starts (shared across all sessions)
# =============================================================================

# ── Required packages ─────────────────────────────────────────────────────────
library(shiny)
library(bslib)
library(shinyjs)
library(leaflet)
library(dplyr)
library(DT)
library(scales)
library(htmltools)
library(plotly)      # interactive PCA scatter (Phase 2 Goal 3)

# Phase 5 additions
# rmarkdown is used for HTML conservation report generation (mod_export.R)
# leaflet.extras2 provides the client-side map screenshot button
# openxlsx generates the multi-sheet Excel workbook export
# pagedown converts the HTML report to PDF (requires Chrome on the server)
# All are optional — the app degrades gracefully if they are absent
if (requireNamespace("rmarkdown",       quietly = TRUE)) library(rmarkdown)
if (requireNamespace("leaflet.extras2", quietly = TRUE)) library(leaflet.extras2)
if (requireNamespace("openxlsx",        quietly = TRUE)) library(openxlsx)
# pagedown is loaded on-demand inside mod_export.R (requireNamespace check)

# vcfR is only needed when users upload real VCF files
# (demo mode works without it)

# ── Source modules and helpers ────────────────────────────────────────────────
source("R/genomics.R")
source("R/utils.R")
source("R/mod_upload.R")
source("R/mod_map.R")
source("R/mod_stats_panel.R")
source("R/mod_structure.R")
source("R/mod_pca.R")         # Phase 2 Goal 3
source("R/mod_admixture.R")   # Phase 2 Goal 4
source("R/mod_froh.R")        # Phase 3 Goal 1 — F-ROH inbreeding
source("R/mod_ne.R")          # Phase 3 Goal 2 — Effective population size (Ne)
source("R/mod_env_upload.R")  # Phase 4 Goal 1 — Environmental data upload
source("R/mod_gea.R")         # Phase 4 Goal 2 — RDA + LFMM2 GEA
source("R/mod_offset.R")      # Phase 4 Goal 3 — Genomic offset / climate vulnerability
source("R/mod_export.R")      # Phase 5     — Export: HTML report, CSV, map PNG
source("R/mod_help.R")        # Phase 5     — In-app Help tab

# ── App metadata ──────────────────────────────────────────────────────────────
APP_TITLE   <- "PopGen Map"
APP_VERSION <- "0.6.0"
APP_PHASE   <- "Phase 5 — Export & Documentation"

# HE_THRESHOLDS is defined in R/genomics.R (sourced above) so it is available
# to all modules without duplication.