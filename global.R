
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
source("R/mod_froh.R")       # Phase 3 Goal 1 — F-ROH inbreeding
source("R/mod_ne.R")         # Phase 3 Goal 2 — Effective population size (Ne)

# ── App metadata ──────────────────────────────────────────────────────────────
APP_TITLE   <- "PopGen Map"
APP_VERSION <- "0.4.0"
APP_PHASE   <- "Phase 3 — Inbreeding & Ne"

# HE_THRESHOLDS is defined in R/genomics.R (sourced above) so it is available
# to all modules without duplication.