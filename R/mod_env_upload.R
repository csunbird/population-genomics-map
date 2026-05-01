# =============================================================================
# mod_env_upload.R — Environmental data upload module (Phase 4, Goal 1)
#
# Provides two input paths for per-population environmental data:
#   1. CSV upload: user-supplied table (population + numeric env columns)
#   2. WorldClim auto-fetch: downloads BIO variables via geodata package
#
# And two input paths for future climate data (used for genomic offset):
#   1. CSV upload: future projected environmental conditions
#   2. CMIP6 auto-fetch: downloads CMIP6 delta-downscaled projections
#
# Demo mode: automatically loads env_data and future_env_data from
# generate_demo_data() so Phase 4 tabs work without requiring uploads.
#
# Returns reactive:
#   $env_data        — data.frame (population + env vars) or NULL
#   $future_env_data — data.frame (population + env vars) or NULL
#   $env_vars        — character vector of selected env variable names
#   $ready           — logical: TRUE when current env data is loaded & valid
#   $scenario_info   — descriptive text
#   $is_worldclim    — logical: TRUE when auto-fetched from WorldClim
# =============================================================================

# ── Constants ─────────────────────────────────────────────────────────────────

# Standard WorldClim BIO variable names and descriptions
WORLDCLIM_BIO_CHOICES <- c(
  "BIO1  — Annual Mean Temperature"       = "1",
  "BIO2  — Mean Diurnal Range"            = "2",
  "BIO3  — Isothermality"                 = "3",
  "BIO4  — Temperature Seasonality"       = "4",
  "BIO5  — Max Temp of Warmest Month"     = "5",
  "BIO6  — Min Temp of Coldest Month"     = "6",
  "BIO7  — Temperature Annual Range"      = "7",
  "BIO8  — Mean Temp of Wettest Quarter"  = "8",
  "BIO9  — Mean Temp of Driest Quarter"   = "9",
  "BIO10 — Mean Temp of Warmest Quarter"  = "10",
  "BIO11 — Mean Temp of Coldest Quarter"  = "11",
  "BIO12 — Annual Precipitation"          = "12",
  "BIO13 — Precip of Wettest Month"       = "13",
  "BIO14 — Precip of Driest Month"        = "14",
  "BIO15 — Precipitation Seasonality"     = "15",
  "BIO16 — Precip of Wettest Quarter"     = "16",
  "BIO17 — Precip of Driest Quarter"      = "17",
  "BIO18 — Precip of Warmest Quarter"     = "18",
  "BIO19 — Precip of Coldest Quarter"     = "19"
)

# ── UI ────────────────────────────────────────────────────────────────────────

#' Environmental upload module UI
#'
#' Renders inside the sidebar, below the main VCF/metadata upload section.
#'
#' @param id Module namespace ID
envUploadUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ── Placeholder shown before genomic data is loaded ─────────────────────
    uiOutput(ns("env_load_placeholder")),

    # ── Main panel (hidden until genomic data loaded) ────────────────────────
    shinyjs::hidden(
      div(
        id = ns("env_panel"),

        tags$hr(style = "margin: 16px 0 10px;"),

        tags$h6(
          tagList(icon("leaf"), " Environmental data (Phase 4)"),
          style = "font-weight:700; color:#1B3A5C; margin-bottom:8px; font-size:0.92em;"
        ),

        # ── Current env source selector ──────────────────────────────────────
        radioButtons(
          ns("cur_env_source"),
          label = tags$span(
            style = "font-size:0.85em; font-weight:600;",
            "Current climate / environment source"
          ),
          choices  = c(
            "Upload CSV"              = "csv",
            "Auto-fetch WorldClim"    = "worldclim"
          ),
          selected = "csv",
          inline   = TRUE
        ),

        # ── CSV upload path ──────────────────────────────────────────────────
        shinyjs::hidden(
          div(
            id = ns("cur_csv_panel"),
            tags$p(
              class = "upload-hint",
              HTML(paste0(
                "One row per population. Required column: <code>population</code>. ",
                "All other numeric columns are treated as environmental variables ",
                "(e.g. BIO1, Temperature, Elevation)."
              ))
            ),
            fileInput(
              ns("env_csv"),
              label       = NULL,
              accept      = ".csv",
              buttonLabel = "Browse…",
              placeholder = "No file selected"
            ),
            downloadButton(
              ns("dl_env_template"),
              label = tagList(icon("download"), " Env CSV template"),
              class = "btn btn-outline-secondary btn-sm w-100"
            )
          )
        ),

        # ── WorldClim auto-fetch path ────────────────────────────────────────
        shinyjs::hidden(
          div(
            id = ns("cur_worldclim_panel"),
            tags$p(
              class = "upload-hint",
              HTML(paste0(
                "Downloads WorldClim 2.1 BIO variables at population centroids. ",
                "Requires an internet connection. Rasters are cached locally. ",
                "<strong>Note:</strong> temperature variables (BIO1–BIO11) are in ",
                "<strong>°C × 10</strong> (e.g. BIO1 = 250 means 25.0 °C). ",
                "Divide by 10 before reporting or comparing against direct measurements."
              ))
            ),
            tags$div(
              style = "font-size:0.82em; font-weight:600; margin-bottom:4px;",
              "Select BIO variables:"
            ),
            checkboxGroupInput(
              ns("wc_bio_vars"),
              label    = NULL,
              choices  = WORLDCLIM_BIO_CHOICES,
              selected = c("1", "4", "12", "15"),
              inline   = FALSE
            ),
            actionButton(
              ns("fetch_worldclim"),
              label = tagList(icon("cloud-download-alt"), " Fetch from WorldClim"),
              class = "btn btn-info btn-sm w-100"
            )
          )
        ),

        # ── Current env status ───────────────────────────────────────────────
        uiOutput(ns("cur_env_status")),

        tags$hr(style = "margin: 14px 0 10px;"),

        # ── Future climate source selector ───────────────────────────────────
        tags$div(
          style = "font-size:0.85em; font-weight:600; margin-bottom:6px;",
          "Future climate source (for genomic offset)"
        ),

        radioButtons(
          ns("fut_env_source"),
          label    = NULL,
          choices  = c(
            "None (skip offset analysis)"  = "none",
            "Upload CSV"                   = "csv",
            "Auto-fetch CMIP6"             = "cmip6"
          ),
          selected = "none",
          inline   = FALSE
        ),

        # ── Future CSV upload ────────────────────────────────────────────────
        shinyjs::hidden(
          div(
            id = ns("fut_csv_panel"),
            tags$p(
              class = "upload-hint",
              HTML(paste0(
                "Same format as current env CSV. ",
                "Must contain the same environmental variable columns."
              ))
            ),
            fileInput(
              ns("fut_env_csv"),
              label       = NULL,
              accept      = ".csv",
              buttonLabel = "Browse…",
              placeholder = "No file selected"
            )
          )
        ),

        # ── CMIP6 auto-fetch options ──────────────────────────────────────────
        shinyjs::hidden(
          div(
            id = ns("fut_cmip6_panel"),
            tags$p(
              class = "upload-hint",
              HTML(paste0(
                "Fetches WorldClim CMIP6 delta-downscaled projections ",
                "(model: <strong>ACCESS-CM2</strong>). ",
                "Uses the same BIO variables as the current WorldClim data. ",
                "Requires current WorldClim data to have been fetched first."
              ))
            ),
            div(
              class = "d-flex gap-3 flex-wrap mb-2",
              div(
                style = "flex: 1 1 120px;",
                selectInput(
                  ns("cmip6_ssp"),
                  label  = tags$span(style = "font-size:0.83em;", "SSP scenario"),
                  choices = c(
                    "SSP2-4.5 (intermediate)" = "245",
                    "SSP5-8.5 (high emissions)" = "585"
                  ),
                  selected = "585",
                  width = "100%"
                )
              ),
              div(
                style = "flex: 1 1 120px;",
                selectInput(
                  ns("cmip6_period"),
                  label   = tags$span(style = "font-size:0.83em;", "Time period"),
                  choices = c(
                    "2041–2060 (mid-century)" = "2041-2060",
                    "2061–2080 (end-century)" = "2061-2080"
                  ),
                  selected = "2041-2060",
                  width = "100%"
                )
              )
            ),
            actionButton(
              ns("fetch_cmip6"),
              label = tagList(icon("cloud-download-alt"), " Fetch CMIP6 projections"),
              class = "btn btn-warning btn-sm w-100"
            )
          )
        ),

        # ── Future env status ────────────────────────────────────────────────
        uiOutput(ns("fut_env_status")),

        # ── Active env data indicator ────────────────────────────────────────
        uiOutput(ns("env_active_indicator"))
      )
    )
  )
}

# ── Server ─────────────────────────────────────────────────────────────────────

#' Environmental upload module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive returning the upload module's data list
#'                    ($ready, $gt, $metadata, $is_demo, $scenario_info)
#' @return Reactive returning list:
#'   $env_data        — current env data frame (population + env vars), or NULL
#'   $future_env_data — future env data frame, or NULL
#'   $env_vars        — character vector of env variable names
#'   $ready           — logical
#'   $scenario_info   — descriptive string
#'   $is_worldclim    — logical
envUploadServer <- function(id, upload_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Reactive state ──────────────────────────────────────────────────────
    rv <- reactiveValues(
      env_data        = NULL,
      future_env_data = NULL,
      env_vars        = character(0),
      cur_status_type = NULL,
      cur_status_msg  = NULL,
      fut_status_type = NULL,
      fut_status_msg  = NULL,
      scenario_info   = NULL,
      is_worldclim    = FALSE
    )

    # ── Show/hide main panel based on genomic data availability ─────────────
    output$env_load_placeholder <- renderUI({
      d <- upload_data()
      if (isTRUE(d$ready)) return(NULL)
      div(
        class = "alert alert-secondary",
        style = "font-size:0.82em; padding:8px 12px; margin-top:8px;",
        icon("leaf"), " ",
        "Environmental data (Phase 4) will appear here after loading genomic data."
      )
    })

    observeEvent(upload_data(), {
      shinyjs::toggle("env_panel", condition = isTRUE(upload_data()$ready))
    })

    # ── Auto-load demo env data when demo mode is active ────────────────────
    observeEvent(upload_data(), {
      d <- upload_data()
      if (!isTRUE(d$ready)) {
        # Reset when data is cleared
        rv$env_data        <- NULL
        rv$future_env_data <- NULL
        rv$env_vars        <- character(0)
        rv$scenario_info   <- NULL
        rv$is_worldclim    <- FALSE
        rv$cur_status_type <- NULL
        rv$cur_status_msg  <- NULL
        rv$fut_status_type <- NULL
        rv$fut_status_msg  <- NULL
        return()
      }
      if (isTRUE(d$is_demo)) {
        demo <- generate_demo_data()
        env_vars_demo <- setdiff(names(demo$env_data), "population")
        rv$env_data        <- demo$env_data
        rv$future_env_data <- demo$future_env_data
        rv$env_vars        <- env_vars_demo
        rv$is_worldclim    <- FALSE
        rv$scenario_info   <- paste0(
          "Demo env data: ", length(env_vars_demo), " variables — ",
          paste(env_vars_demo, collapse = ", "), ". ",
          "Future climate loaded (SSP5-8.5, 2041-2060, simulated)."
        )
        rv$cur_status_type <- "success"
        rv$cur_status_msg  <- paste0(
          "✅ Demo env data loaded: ",
          length(unique(demo$env_data$population)), " populations, ",
          length(env_vars_demo), " variables (",
          paste(env_vars_demo, collapse = ", "), ")."
        )
        rv$fut_status_type <- "success"
        rv$fut_status_msg  <- paste0(
          "✅ Demo future climate loaded ",
          "(SSP5-8.5 equivalent, 2041-2060, simulated)."
        )
      }
    })

    # ── Source panel show/hide: current env ─────────────────────────────────
    observeEvent(input$cur_env_source, {
      shinyjs::toggle("cur_csv_panel",       condition = input$cur_env_source == "csv")
      shinyjs::toggle("cur_worldclim_panel", condition = input$cur_env_source == "worldclim")
    }, ignoreInit = FALSE)

    # ── Source panel show/hide: future env ──────────────────────────────────
    observeEvent(input$fut_env_source, {
      shinyjs::toggle("fut_csv_panel",   condition = input$fut_env_source == "csv")
      shinyjs::toggle("fut_cmip6_panel", condition = input$fut_env_source == "cmip6")
    }, ignoreInit = FALSE)

    # ── Current env CSV upload ───────────────────────────────────────────────
    observeEvent(input$env_csv, {
      d <- upload_data()
      req(d$ready, input$env_csv)
      rv$cur_status_type <- "info"
      rv$cur_status_msg  <- "⏳ Parsing environmental CSV…"

      tryCatch({
        df <- utils::read.csv(input$env_csv$datapath, stringsAsFactors = FALSE)

        val <- validate_env_data(df, d$metadata)

        if (!val$valid) {
          rv$cur_status_type <- "error"
          rv$cur_status_msg  <- paste0(
            "❌ Env CSV error: ", paste(val$errors, collapse = " | ")
          )
          return()
        }

        # Keep only matched populations (drop unmatched rows)
        df_clean <- df[as.character(df$population) %in% val$matched_pops, , drop = FALSE]

        rv$env_data        <- df_clean
        rv$env_vars        <- val$env_vars
        rv$is_worldclim    <- FALSE

        warn_txt <- if (length(val$warnings) > 0)
          paste0(" ⚠️ ", paste(val$warnings, collapse = "; "))
        else ""

        rv$scenario_info <- paste0(
          "CSV: ", length(val$matched_pops), " populations, ",
          length(val$env_vars), " variables: ",
          paste(head(val$env_vars, 5), collapse = ", "),
          if (length(val$env_vars) > 5) " …" else "."
        )
        rv$cur_status_type <- if (length(val$warnings) > 0) "warning" else "success"
        rv$cur_status_msg  <- paste0(
          "✅ Env CSV loaded: ", length(val$matched_pops), " populations, ",
          length(val$env_vars), " variable(s) — ",
          paste(head(val$env_vars, 5), collapse = ", "),
          if (length(val$env_vars) > 5) " …" else ".", warn_txt
        )
        # Clear future data if variables changed
        if (!is.null(rv$future_env_data)) {
          if (!all(val$env_vars %in% names(rv$future_env_data))) {
            rv$future_env_data <- NULL
            rv$fut_status_type <- "warning"
            rv$fut_status_msg  <- paste0(
              "⚠️ Future env cleared — variable mismatch with new current env. ",
              "Please re-upload future env CSV."
            )
          }
        }

      }, error = function(e) {
        rv$cur_status_type <- "error"
        rv$cur_status_msg  <- paste0("❌ Error reading CSV: ", e$message)
      })
    })

    # ── WorldClim auto-fetch ─────────────────────────────────────────────────
    observeEvent(input$fetch_worldclim, {
      d <- upload_data()
      req(d$ready)

      bio_nums <- as.integer(input$wc_bio_vars)
      if (length(bio_nums) == 0) {
        rv$cur_status_type <- "error"
        rv$cur_status_msg  <- "❌ Please select at least one BIO variable."
        return()
      }

      shinyjs::disable("fetch_worldclim")
      on.exit(shinyjs::enable("fetch_worldclim"), add = TRUE)

      rv$cur_status_type <- "info"
      rv$cur_status_msg  <- paste0(
        "⏳ Fetching WorldClim 2.1 BIO",
        paste(bio_nums, collapse = "/"), " (2.5 arc-min) — ",
        "may take a minute on first run…"
      )

      tryCatch({
        env_df <- fetch_worldclim_env(
          metadata = d$metadata,
          vars     = bio_nums,
          path     = tempdir()
        )
        env_vars_wc <- paste0("BIO", bio_nums)

        val <- validate_env_data(env_df, d$metadata)

        rv$env_data        <- env_df
        rv$env_vars        <- env_vars_wc
        rv$is_worldclim    <- TRUE
        rv$scenario_info   <- paste0(
          "WorldClim: ", length(unique(env_df$population)), " populations, ",
          length(env_vars_wc), " BIO variables."
        )
        rv$cur_status_type <- "success"
        rv$cur_status_msg  <- paste0(
          "✅ WorldClim data fetched: ",
          length(unique(env_df$population)), " populations, variables: ",
          paste(env_vars_wc, collapse = ", "), "."
        )
        if (length(val$warnings) > 0) {
          rv$cur_status_msg <- paste0(
            rv$cur_status_msg, " ⚠️ ",
            paste(val$warnings, collapse = "; ")
          )
          rv$cur_status_type <- "warning"
        }
      }, error = function(e) {
        rv$cur_status_type <- "error"
        rv$cur_status_msg  <- paste0(
          "❌ WorldClim fetch failed: ", e$message,
          " — Check your internet connection, or upload a CSV instead."
        )
      })
    })

    # ── Future env CSV upload ────────────────────────────────────────────────
    observeEvent(input$fut_env_csv, {
      req(input$fut_env_csv)
      d <- upload_data()
      req(d$ready)

      rv$fut_status_type <- "info"
      rv$fut_status_msg  <- "⏳ Parsing future env CSV…"

      tryCatch({
        df  <- utils::read.csv(input$fut_env_csv$datapath, stringsAsFactors = FALSE)
        val <- validate_env_data(df, d$metadata)

        if (!val$valid) {
          rv$fut_status_type <- "error"
          rv$fut_status_msg  <- paste0(
            "❌ Future env CSV error: ", paste(val$errors, collapse = " | ")
          )
          return()
        }

        # Check variable alignment with current env
        if (!is.null(rv$env_vars) && length(rv$env_vars) > 0) {
          missing_vars <- setdiff(rv$env_vars, val$env_vars)
          if (length(missing_vars) > 0) {
            rv$fut_status_type <- "error"
            rv$fut_status_msg  <- paste0(
              "❌ Future env CSV is missing variables that appear in the current env: ",
              paste(missing_vars, collapse = ", "), ". ",
              "Both files must have the same environmental variable columns."
            )
            return()
          }
        }

        df_clean <- df[as.character(df$population) %in% val$matched_pops, , drop = FALSE]
        rv$future_env_data <- df_clean
        rv$fut_status_type <- "success"
        rv$fut_status_msg  <- paste0(
          "✅ Future env CSV loaded: ",
          length(val$matched_pops), " populations, ",
          length(val$env_vars), " variable(s).",
          if (length(val$warnings) > 0)
            paste0(" ⚠️ ", paste(val$warnings, collapse = "; "))
          else ""
        )

      }, error = function(e) {
        rv$fut_status_type <- "error"
        rv$fut_status_msg  <- paste0("❌ Error reading future CSV: ", e$message)
      })
    })

    # ── CMIP6 auto-fetch ─────────────────────────────────────────────────────
    observeEvent(input$fetch_cmip6, {
      d <- upload_data()
      req(d$ready)

      if (is.null(rv$env_data) || !isTRUE(rv$is_worldclim)) {
        rv$fut_status_type <- "error"
        rv$fut_status_msg  <- paste0(
          "❌ CMIP6 fetch requires current WorldClim data to be loaded first. ",
          "Please fetch WorldClim current data before fetching CMIP6 projections."
        )
        return()
      }

      bio_nums <- as.integer(
        sub("BIO", "", rv$env_vars[grepl("^BIO[0-9]+$", rv$env_vars)])
      )
      if (length(bio_nums) == 0) {
        rv$fut_status_type <- "error"
        rv$fut_status_msg  <- "❌ Could not determine BIO variable numbers from current WorldClim data."
        return()
      }

      shinyjs::disable("fetch_cmip6")
      on.exit(shinyjs::enable("fetch_cmip6"), add = TRUE)

      ssp    <- input$cmip6_ssp    %||% "585"
      period <- input$cmip6_period %||% "2041-2060"
      ssp_label <- if (ssp == "245") "SSP2-4.5" else "SSP5-8.5"

      rv$fut_status_type <- "info"
      rv$fut_status_msg  <- paste0(
        "⏳ Fetching CMIP6 (", ssp_label, ", ", period, ") — ",
        "may take a minute on first run…"
      )

      tryCatch({
        fut_df <- fetch_cmip6_env(
          metadata = d$metadata,
          vars     = bio_nums,
          ssp      = ssp,
          period   = period,
          path     = tempdir()
        )
        rv$future_env_data <- fut_df
        rv$fut_status_type <- "success"
        rv$fut_status_msg  <- paste0(
          "✅ CMIP6 projections fetched (ACCESS-CM2, ", ssp_label, ", ", period, "): ",
          length(unique(fut_df$population)), " populations."
        )

      }, error = function(e) {
        rv$fut_status_type <- "error"
        rv$fut_status_msg  <- paste0(
          "❌ CMIP6 fetch failed: ", e$message,
          " — Check your internet connection, or upload a future env CSV instead."
        )
      })
    })

    # ── Current env status display ───────────────────────────────────────────
    output$cur_env_status <- renderUI({
      req(!is.null(rv$cur_status_msg))
      status_alert(rv$cur_status_msg, rv$cur_status_type)
    })

    # ── Future env status display ────────────────────────────────────────────
    output$fut_env_status <- renderUI({
      req(!is.null(rv$fut_status_msg))
      status_alert(rv$fut_status_msg, rv$fut_status_type)
    })

    # ── Active env data indicator ────────────────────────────────────────────
    output$env_active_indicator <- renderUI({
      req(!is.null(rv$env_data))
      n_vars <- length(rv$env_vars)
      n_pops <- length(unique(rv$env_data$population))
      has_fut <- !is.null(rv$future_env_data)

      div(
        style = paste0(
          "font-size:0.78em; color:#555; background:#f8f9fa;",
          " border:1px solid #dee2e6; border-radius:4px;",
          " padding:5px 8px; margin-top:8px; line-height:1.5;"
        ),
        icon("circle", style = "color:#388E3C; font-size:8px; margin-right:4px; vertical-align:middle;"),
        tags$strong("Current env: "),
        paste0(n_pops, " pop × ", n_vars, " var"),
        tags$br(),
        icon("circle",
             style = paste0("color:", if (has_fut) "#388E3C" else "#AAAAAA",
                            "; font-size:8px; margin-right:4px; vertical-align:middle;")),
        tags$strong("Future env: "),
        if (has_fut) paste0(length(unique(rv$future_env_data$population)), " pop (offset ready)")
        else tags$em("not loaded (offset disabled)")
      )
    })

    # ── Environmental CSV template download ──────────────────────────────────
    output$dl_env_template <- downloadHandler(
      filename = "popgen_env_template.csv",
      content  = function(file) {
        template <- data.frame(
          population  = c("Pop_A", "Pop_B", "Pop_C", "Pop_D", "Pop_E"),
          BIO1        = c(248L, 252L, 262L, 244L, 268L),
          BIO12       = c(2800L, 2950L, 3400L, 2400L, 3600L),
          Elevation_m = c(420L, 320L, 35L, 380L, 15L),
          stringsAsFactors = FALSE
        )
        utils::write.csv(template, file, row.names = FALSE)
      }
    )

    # ── Return reactive ──────────────────────────────────────────────────────
    return(reactive(list(
      env_data        = rv$env_data,
      future_env_data = rv$future_env_data,
      env_vars        = rv$env_vars,
      ready           = !is.null(rv$env_data) && length(rv$env_vars) > 0,
      scenario_info   = rv$scenario_info,
      is_worldclim    = isTRUE(rv$is_worldclim)
    )))
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Small Bootstrap alert div for env module status messages
#'
#' @param msg  Character message string
#' @param type One of "success", "warning", "error", "info"
status_alert <- function(msg, type) {
  cls <- switch(
    type %||% "info",
    "success" = "alert alert-success",
    "warning" = "alert alert-warning",
    "error"   = "alert alert-danger",
    "alert alert-info"
  )
  div(
    class = cls,
    style = "font-size:0.80em; padding:7px 12px; margin-top:6px;",
    msg
  )
}
