# =============================================================================
# mod_gea.R — Genotype-Environment Association (GEA) module (Phase 4, Goal 2)
#
# Implements two complementary GEA methods:
#   1. Partial RDA (Redundancy Analysis) — multivariate, identifies sets of
#      candidate adaptive loci; geography partialled out via Condition(lat+lon)
#   2. LFMM2 — per-locus univariate GEA controlling for K latent factors
#      representing background population structure
#
# Displays:
#   1. Summary cards: % variance constrained, n significant loci, top env var
#   2. RDA biplot: sites (populations) + biplot arrows (env), Plotly
#   3. LFMM2 Manhattan-style plot: -log10(p_adj) per locus, per env variable
#   4. Results table: significant loci by env variable with export
#
# Design:
#   - RDA and LFMM2 run independently (separate action buttons); each is
#     cached in its own reactive value so changing the LFMM2 K slider does
#     not rerun RDA and vice versa.
#   - LFMM2 requires 'lfmm' CRAN package; RDA requires 'vegan'.
#   - Both methods depend on env_data from the env upload module.
#   - In demo mode, both run automatically without requiring user action.
# =============================================================================

#' GEA module UI
#'
#' @param id Module namespace ID
geaUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("load_placeholder")),
    uiOutput(ns("scenario_banner")),
    uiOutput(ns("summary_cards")),

    shinyjs::hidden(
      div(
        id = ns("gea_ready_content"),

        tags$hr(style = "margin: 12px 0;"),

        # ── Method tabs ─────────────────────────────────────────────────────
        tabsetPanel(
          id = ns("gea_tabs"),
          type = "pills",

          # ── RDA tab ───────────────────────────────────────────────────────
          tabPanel(
            title = tagList(icon("chart-scatter"), " Partial RDA"),
            value = "rda",

            div(
              class = "mt-3",

              # Controls row
              div(
                class = "d-flex align-items-end gap-3 flex-wrap mb-2",
                div(
                  style = "flex: 0 0 auto;",
                  div(
                    style = "font-size:0.83em; font-weight:600; margin-bottom:2px;",
                    "Environmental variables for RDA:"
                  ),
                  uiOutput(ns("rda_var_selector"))
                ),
                div(
                  style = "flex: 0 0 auto;",
                  actionButton(
                    ns("run_rda"),
                    label = tagList(icon("play"), " Run RDA"),
                    class = "btn btn-primary btn-sm"
                  )
                ),
                div(
                  style = "flex: 0 0 auto;",
                  shinyjs::disabled(
                    downloadButton(
                      ns("dl_rda"),
                      label = tagList(icon("download"), " Export"),
                      class = "btn btn-outline-secondary btn-sm"
                    )
                  )
                )
              ),

              uiOutput(ns("rda_status")),

              # RDA biplot
              shinyjs::hidden(
                div(
                  id = ns("rda_plot_panel"),
                  plotly::plotlyOutput(ns("rda_biplot"), height = "480px"),

                  # Variance partitioning summary
                  uiOutput(ns("varpart_summary")),

                  # Top candidate loci table
                  tags$h6(
                    "Top candidate adaptive loci (highest |RDA1 loading|)",
                    style = "font-weight:600; color:#1B3A5C; margin: 14px 0 6px;"
                  ),
                  tags$p(
                    style = "font-size:0.76em; color:#888; margin: -2px 0 6px;",
                    HTML(paste0(
                      "Ranked by absolute loading on RDA1. ",
                      "For datasets with multiple environmental variables, RDA2 and higher axes ",
                      "may capture additional independent adaptive signals — ",
                      "inspect the biplot arrows to identify which env variables drive each axis."
                    ))
                  ),
                  DT::dataTableOutput(ns("rda_loci_table"))
                )
              )
            )
          ),

          # ── LFMM2 tab ────────────────────────────────────────────────────
          tabPanel(
            title = tagList(icon("dna"), " LFMM2"),
            value = "lfmm2",

            div(
              class = "mt-3",

              # Controls row
              div(
                class = "d-flex align-items-end gap-3 flex-wrap mb-2",
                div(
                  style = "flex: 0 0 auto;",
                  div(
                    style = "font-size:0.83em; font-weight:600; margin-bottom:2px;",
                    "Variables to test:"
                  ),
                  uiOutput(ns("lfmm_var_selector"))
                ),
                div(
                  style = "flex: 0 0 auto;",
                  tags$label(
                    HTML("K latent factors:"),
                    class  = "me-1 fw-semibold",
                    style  = "font-size:0.83em;"
                  ),
                  sliderInput(
                    ns("lfmm_k"),
                    label  = NULL,
                    min = 1L, max = 10L, value = 3L, step = 1L,
                    ticks = FALSE, width = "130px"
                  ),
                  tags$p(
                    style = "font-size:0.72em; color:#888; margin:-4px 0 0;",
                    HTML("Recommended: K &le; n<sub>pop</sub> &minus; 1. Typical start: 2&ndash;3.")
                  )
                ),
                div(
                  style = "flex: 0 0 auto;",
                  sliderInput(
                    ns("lfmm_alpha"),
                    label  = tags$span(style = "font-size:0.83em; font-weight:600;",
                                       "FDR threshold (α):"),
                    min = 0.01, max = 0.20, value = 0.05, step = 0.01,
                    ticks = FALSE, width = "150px"
                  )
                ),
                div(
                  style = "flex: 0 0 auto;",
                  actionButton(
                    ns("run_lfmm2"),
                    label = tagList(icon("play"), " Run LFMM2"),
                    class = "btn btn-primary btn-sm"
                  )
                ),
                div(
                  style = "flex: 0 0 auto;",
                  shinyjs::disabled(
                    downloadButton(
                      ns("dl_lfmm2"),
                      label = tagList(icon("download"), " Export CSV"),
                      class = "btn btn-outline-secondary btn-sm"
                    )
                  )
                )
              ),

              uiOutput(ns("lfmm_status")),

              shinyjs::hidden(
                div(
                  id = ns("lfmm_plot_panel"),

                  # Manhattan plot
                  plotly::plotlyOutput(ns("manhattan_plot"), height = "380px"),

                  # Significant loci table
                  tags$h6(
                    "Significant loci (BH-adjusted p < α)",
                    style = "font-weight:600; color:#1B3A5C; margin: 14px 0 6px;"
                  ),
                  uiOutput(ns("lfmm_table_note")),
                  DT::dataTableOutput(ns("lfmm_table"))
                )
              )
            )
          )
        )
      )
    ),

    # ── Footnotes ─────────────────────────────────────────────────────────────
    tags$p(
      class = "mt-3",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "<b>Partial RDA</b> (Legendre & Legendre 2012): multivariate regression of ",
        "Hellinger-transformed population allele frequencies on environmental variables, ",
        "with latitude and longitude as covariates (Condition) to partial out spatial autocorrelation. ",
        "Candidate adaptive loci are identified by extreme species scores on the constrained axes."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "<b>LFMM2</b> (Caye et al. 2019): per-locus ridge regression controlling for ",
        "K latent factors that capture background population structure. ",
        "P-values are calibrated by the genomic inflation factor (GIF) and corrected ",
        "for multiple testing (Benjamini-Hochberg FDR). ",
        "GIF ≈ 1 = well-calibrated; ",
        "GIF 1.5–2 = mild inflation (check batch effects or unmodelled structure); ",
        "GIF > 2 = strong undercorrection (increase K or add geographic covariates); ",
        "GIF < 0.8 = over-correction (consider reducing K)."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML(paste0(
        "&#9432; Both methods identify <em>candidate</em> adaptive loci, not causative variants. ",
        "Candidates require functional annotation and validation. ",
        "RDA is best for multivariate patterns; LFMM2 is best for single-variable precision. ",
        "Results depend on the choice and completeness of environmental predictors."
      ))
    )
  )
}

#' GEA module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive returning the upload module's data list
#' @param env_data    Reactive returning the env upload module's data list
#'                    ($env_data, $future_env_data, $env_vars, $ready)
#' @return Reactive returning list ($rda_result, $lfmm_result) for downstream use
geaServer <- function(id, upload_data, env_data) {
  moduleServer(id, function(input, output, session) {

    # ── Reactive state ──────────────────────────────────────────────────────
    rda_rv  <- reactiveValues(result = NULL, status_type = NULL, status_msg = NULL)
    lfmm_rv <- reactiveValues(result = NULL, status_type = NULL, status_msg = NULL)

    # ── Helper: are all inputs ready? ───────────────────────────────────────
    all_ready <- reactive({
      d <- upload_data(); e <- env_data()
      isTRUE(d$ready) && isTRUE(e$ready)
    })

    # ── Load placeholder ─────────────────────────────────────────────────────
    output$load_placeholder <- renderUI({
      if (isTRUE(all_ready())) return(NULL)
      d <- upload_data(); e <- env_data()
      msg <- if (!isTRUE(d$ready)) {
        "Load a VCF + metadata CSV from the sidebar to enable GEA analysis."
      } else {
        "Load environmental data from the sidebar (Phase 4 section) to enable GEA analysis."
      }
      div(
        class = "alert alert-secondary",
        style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
        icon("info-circle"), " ", msg
      )
    })

    # ── Show/hide controls ───────────────────────────────────────────────────
    observe({
      shinyjs::toggle("gea_ready_content", condition = isTRUE(all_ready()))
    })

    # ── Scenario banner ──────────────────────────────────────────────────────
    output$scenario_banner <- renderUI({
      d <- upload_data(); req(d$ready, d$scenario_info)
      cls <- if (isTRUE(d$is_demo)) "alert alert-warning" else "alert alert-info"
      div(
        class = cls,
        style = "font-size:0.82em; padding:8px 12px; margin-bottom:10px;",
        icon("info-circle"), " ", d$scenario_info
      )
    })

    # ── Dynamic variable selectors (built from env_vars) ────────────────────
    output$rda_var_selector <- renderUI({
      e <- env_data(); req(e$ready)
      checkboxGroupInput(
        session$ns("rda_vars"),
        label    = NULL,
        choices  = e$env_vars,
        selected = head(e$env_vars, 4),
        inline   = TRUE
      )
    })

    output$lfmm_var_selector <- renderUI({
      e <- env_data(); req(e$ready)
      checkboxGroupInput(
        session$ns("lfmm_vars"),
        label    = NULL,
        choices  = e$env_vars,
        selected = head(e$env_vars, 2),
        inline   = TRUE
      )
    })

    # ── Update LFMM K slider max based on number of populations ─────────────
    observeEvent(upload_data(), {
      d <- upload_data()
      req(d$ready)
      n_pop     <- length(unique(d$metadata$population))
      k_max     <- max(1L, n_pop - 1L)
      k_default <- as.integer(min(3L, k_max))
      updateSliderInput(session, "lfmm_k",
                        max   = min(10L, k_max),
                        value = k_default)
    })

    # ── Auto-run on demo data ─────────────────────────────────────────────────
    observeEvent(list(upload_data(), env_data()), {
      d <- upload_data(); e <- env_data()
      if (!isTRUE(d$ready) || !isTRUE(e$ready) || !isTRUE(d$is_demo)) return()
      if (!is.null(rda_rv$result)) return()   # already run this session
      # Trigger RDA run automatically for demo
      run_rda_analysis(d, e)
    })

    # ── Run RDA ──────────────────────────────────────────────────────────────
    run_rda_analysis <- function(d, e) {
      env_vars_sel <- if (is.null(input$rda_vars) || length(input$rda_vars) == 0)
        head(e$env_vars, 4) else input$rda_vars

      if (length(env_vars_sel) == 0) {
        rda_rv$status_type <- "error"
        rda_rv$status_msg  <- "❌ Select at least one environmental variable for RDA."
        return()
      }
      # Use only populations present in BOTH genotype metadata AND env CSV —
      # that is the set calc_rda() will actually use.  Using raw metadata count
      # produces a check that passes here but then fails inside calc_rda().
      all_meta_pops <- unique(as.character(d$metadata$population))
      env_pops      <- unique(as.character(e$env_data$population))
      n_pops        <- length(intersect(all_meta_pops, env_pops))

      if (n_pops < 3L) {
        rda_rv$status_type <- "error"
        rda_rv$status_msg  <- paste0(
          "❌ RDA requires at least 3 populations with both genomic and environmental data. ",
          "Only ", n_pops, " population(s) matched."
        )
        return()
      }
      if (n_pops <= length(env_vars_sel) + 2L) {
        rda_rv$status_type <- "warning"
        rda_rv$status_msg  <- paste0(
          "⚠️ Too many environmental variables (", length(env_vars_sel),
          ") for ", n_pops, " matched populations. ",
          "Use at most ", n_pops - 3L, " variable(s)."
        )
        return()
      }

      # Note for single-variable RDA: constrained space is 1D (valid but produces
      # a flat biplot — all populations at y = 0).  We append this to the success
      # message so users understand the ordination is correct, not broken.
      single_var_note <- if (length(env_vars_sel) == 1L) {
        paste0(
          " ℹ️ Single-variable RDA: constrained space is 1D ",
          "— populations plotted on one axis only. ",
          "Select ≥ 2 variables for a 2D biplot."
        )
      } else ""

      rda_rv$result      <- NULL
      rda_rv$status_type <- "info"
      rda_rv$status_msg  <- paste0(
        "⏳ Running partial RDA (", length(env_vars_sel), " variable(s), ",
        n_pops, " populations)…"
      )

      tryCatch({
        res <- calc_rda(
          gt_matrix = d$gt,
          metadata  = d$metadata,
          env_data  = e$env_data,
          env_vars  = env_vars_sel,
          scale_env = TRUE
        )
        rda_rv$result      <- res
        rda_rv$status_type <- if (nchar(single_var_note) > 0) "warning" else "success"
        rda_rv$status_msg  <- paste0(
          "✅ RDA complete: ", res$pct_constrained, "% variance constrained by environment | ",
          res$n_loci, " loci | ", res$n_pops, " populations.",
          single_var_note
        )
        shinyjs::show("rda_plot_panel")
        shinyjs::enable("dl_rda")
      }, error = function(e2) {
        rda_rv$status_type <- "error"
        rda_rv$status_msg  <- paste0("❌ RDA error: ", e2$message)
        shinyjs::hide("rda_plot_panel")
        shinyjs::disable("dl_rda")
      })
    }

    observeEvent(input$run_rda, {
      d <- upload_data(); e <- env_data()
      req(d$ready, e$ready)
      shinyjs::disable("run_rda")
      on.exit(shinyjs::enable("run_rda"), add = TRUE)
      run_rda_analysis(d, e)
    })

    # ── RDA status ────────────────────────────────────────────────────────────
    output$rda_status <- renderUI({
      req(!is.null(rda_rv$status_msg))
      gea_status_alert(rda_rv$status_msg, rda_rv$status_type)
    })

    # ── RDA biplot ────────────────────────────────────────────────────────────
    output$rda_biplot <- plotly::renderPlotly({
      res <- rda_rv$result; req(res)
      d   <- upload_data()

      site_sc <- res$site_scores    # populations × 2
      bp_sc   <- res$biplot_scores  # env vars × 2
      n_axes  <- ncol(site_sc)

      # Guard: vegan may return site scores without rownames when n_pops == 1
      # (edge case); ensure rownames are the population names.
      if (is.null(rownames(site_sc))) rownames(site_sc) <- res$populations

      # Axis labels
      x_lab <- if (!is.null(res$axis_labels) && length(res$axis_labels) >= 1)
        res$axis_labels[1] else "RDA1"
      y_lab <- if (!is.null(res$axis_labels) && length(res$axis_labels) >= 2)
        res$axis_labels[2] else "RDA2"

      # Population colours — re-use PCA palette
      pops    <- res$populations
      palette <- stats::setNames(
        PCA_PALETTE[((seq_along(pops) - 1L) %% length(PCA_PALETTE)) + 1L],
        pops
      )
      site_cols <- palette[rownames(site_sc)]
      # Fallback colour for any unmatched population
      site_cols[is.na(site_cols)] <- "#888888"

      # Site (population) scatter
      s1 <- if (n_axes >= 1) site_sc[, 1] else rep(0, nrow(site_sc))
      s2 <- if (n_axes >= 2) site_sc[, 2] else rep(0, nrow(site_sc))

      fig <- plotly::plot_ly() |>
        plotly::add_markers(
          x         = s1,
          y         = s2,
          text      = pops,
          hoverinfo = "text",
          marker    = list(
            color  = site_cols,
            size   = 12,
            line   = list(color = "white", width = 1)
          ),
          name      = "Population"
        )

      # Biplot arrows (env vectors) — scale to 80% of max site range
      site_range <- max(abs(c(s1, s2)), na.rm = TRUE)
      bp_range   <- if (!is.null(bp_sc) && nrow(bp_sc) > 0)
        max(abs(bp_sc), na.rm = TRUE) else 1
      arrow_scale <- 0.8 * site_range / max(bp_range, 1e-6)

      bp1 <- if (n_axes >= 1) bp_sc[, 1] * arrow_scale else rep(0, nrow(bp_sc))
      bp2 <- if (n_axes >= 2) bp_sc[, 2] * arrow_scale else rep(0, nrow(bp_sc))

      if (!is.null(bp_sc) && nrow(bp_sc) > 0) {
        for (i in seq_len(nrow(bp_sc))) {
          fig <- fig |>
            plotly::add_segments(
              x = 0, y = 0,
              xend = bp1[i], yend = bp2[i],
              line = list(color = "#D32F2F", width = 2),
              showlegend = FALSE,
              hoverinfo  = "none"
            ) |>
            plotly::add_annotations(
              x = bp1[i] * 1.1, y = bp2[i] * 1.1,
              text       = rownames(bp_sc)[i],
              showarrow  = FALSE,
              font       = list(size = 11, color = "#D32F2F")
            )
        }
      }

      fig |>
        plotly::layout(
          xaxis = list(title = x_lab, zeroline = TRUE, zerolinecolor = "#DDDDDD",
                       showgrid = TRUE, gridcolor = "#F5F5F5"),
          yaxis = list(title = y_lab, zeroline = TRUE, zerolinecolor = "#DDDDDD",
                       showgrid = TRUE, gridcolor = "#F5F5F5"),
          paper_bgcolor = "white",
          plot_bgcolor  = "white",
          showlegend    = FALSE,
          font = list(family = "Helvetica Neue, Helvetica, Arial, sans-serif", size = 12)
        ) |>
        plotly::config(
          displayModeBar = TRUE,
          modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
          displaylogo = FALSE,
          toImageButtonOptions = list(
            format   = "png",
            scale    = 2,
            filename = "popgen_rda_biplot"
          )
        )
    })

    # ── Variance partitioning summary ────────────────────────────────────────
    output$varpart_summary <- renderUI({
      res <- rda_rv$result; req(res)

      vp_text <- if (!is.null(res$varpart)) {
        fractions <- tryCatch({
          ind <- res$varpart$part$indfract
          # Column "Adj.R.squared" (col 3): adjusted fractions for each partition.
          # Row 1 = environment alone, Row 2 = geography alone, Row 3 = shared,
          # Row 4 = residual.  Multiply by 100 for %.  Negative adj-R² is possible
          # and means the fraction explains less than chance; display as 0%.
          adj_r2    <- ind[, "Adj.R.squared"]
          env_pure  <- round(max(0, adj_r2[1]) * 100, 1)
          geo_pure  <- round(max(0, adj_r2[2]) * 100, 1)
          shared    <- round(adj_r2[3] * 100, 1)          # can be negative (suppression)
          resid     <- round(max(0, adj_r2[4]) * 100, 1)
          paste0(
            "Environment alone: ", env_pure, "% | ",
            "Geography alone: ", geo_pure, "% | ",
            "Shared E+G: ", shared, "% | ",
            "Residual: ", resid, "%"
          )
        }, error = function(e) NULL)
        fractions
      } else NULL

      tags$p(
        style = "font-size:0.80em; color:#555; margin-top:8px;",
        tags$strong("Variance explained: "),
        paste0(res$pct_constrained, "% constrained by environment | "),
        if (!is.null(vp_text)) vp_text else ""
      )
    })

    # ── Top candidate loci table (from RDA) ──────────────────────────────────
    output$rda_loci_table <- DT::renderDataTable({
      res <- rda_rv$result; req(res)
      sp_sc    <- res$species_scores
      top_idx  <- res$top_loci_idx

      display <- data.frame(
        "Locus index"    = top_idx,
        "RDA1 loading"   = round(sp_sc[top_idx, 1], 4),
        "|RDA1 loading|" = round(abs(sp_sc[top_idx, 1]), 4),
        check.names      = FALSE,
        stringsAsFactors = FALSE
      )
      if (ncol(sp_sc) >= 2) {
        display[["RDA2 loading"]] <- round(sp_sc[top_idx, 2], 4)
      }
      # Note: column name is "|RDA1 loading|" — no space after the opening |
      display <- display[order(display[["|RDA1 loading|"]], decreasing = TRUE,
                                na.last = TRUE), ]

      DT::datatable(
        display,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength = 10, dom = "tp", ordering = TRUE,
          columnDefs = list(list(className = "dt-center", targets = 0:2))
        )
      ) |>
        DT::formatRound(columns = c("RDA1 loading", "|RDA1 loading|"), digits = 4)
    })

    # ── RDA export ────────────────────────────────────────────────────────────
    output$dl_rda <- downloadHandler(
      filename = function() paste0("popgen_rda_", Sys.Date(), ".csv"),
      content  = function(file) {
        res <- rda_rv$result; req(res)
        sp_sc <- res$species_scores
        out <- data.frame(
          locus_idx     = seq_len(nrow(sp_sc)),
          rda1_loading  = round(sp_sc[, 1], 4),
          abs_rda1      = round(abs(sp_sc[, 1]), 4),
          stringsAsFactors = FALSE
        )
        if (ncol(sp_sc) >= 2) {
          out$rda2_loading <- round(sp_sc[, 2], 4)
        }
        out$candidate <- out$locus_idx %in% res$top_loci_idx
        utils::write.csv(out, file, row.names = FALSE)
      }
    )

    # ── Summary cards ────────────────────────────────────────────────────────
    output$summary_cards <- renderUI({
      has_rda  <- !is.null(rda_rv$result)
      has_lfmm <- !is.null(lfmm_rv$result)
      if (!has_rda && !has_lfmm) return(NULL)

      pct_const  <- if (has_rda)  paste0(rda_rv$result$pct_constrained, "%") else "—"
      n_cand_rda <- if (has_rda)  length(rda_rv$result$top_loci_idx) else NA
      n_sig_lfmm <- if (has_lfmm) sum(lfmm_rv$result$significant, na.rm = TRUE) else NA

      top_var <- if (has_lfmm) {
        sig <- lfmm_rv$result[lfmm_rv$result$significant, ]
        if (nrow(sig) > 0) {
          tab <- sort(table(sig$env_var), decreasing = TRUE)
          names(tab)[1]
        } else "—"
      } else "—"

      div(
        class = "row g-2 mb-2",
        gea_card("% Variance (env)", pct_const,     "#1B3A5C", "percent"),
        gea_card("RDA candidate loci",
                 if (has_rda) as.character(n_cand_rda) else "—",
                 "#388E3C", "dna"),
        gea_card("Significant loci (LFMM2)",
                 if (has_lfmm) as.character(n_sig_lfmm) else "—",
                 if (!is.na(n_sig_lfmm) && n_sig_lfmm > 0) "#D32F2F" else "#388E3C",
                 "exclamation-circle"),
        gea_card("Top env variable", top_var,        "#984EA3", "thermometer-half")
      )
    })

    # ── LFMM2: run ────────────────────────────────────────────────────────────
    run_lfmm2_analysis <- function(d, e) {
      env_vars_sel <- if (is.null(input$lfmm_vars) || length(input$lfmm_vars) == 0)
        head(e$env_vars, 2) else input$lfmm_vars

      if (length(env_vars_sel) == 0) {
        lfmm_rv$status_type <- "error"
        lfmm_rv$status_msg  <- "❌ Select at least one environmental variable for LFMM2."
        return()
      }

      K_val   <- input$lfmm_k    %||% 3L
      alpha   <- input$lfmm_alpha %||% 0.05
      n_ind   <- nrow(d$gt)
      n_loci  <- ncol(d$gt)

      lfmm_rv$result      <- NULL
      lfmm_rv$status_type <- "info"
      lfmm_rv$status_msg  <- paste0(
        "⏳ Running LFMM2 (K = ", K_val, ", α = ", alpha, ", ",
        length(env_vars_sel), " variable(s), ",
        n_ind, " individuals × ", n_loci, " loci) — this may take 30–120 s…"
      )

      tryCatch({
        res <- calc_lfmm2_gea(
          gt_matrix = d$gt,
          metadata  = d$metadata,
          env_data  = e$env_data,
          env_vars  = env_vars_sel,
          K         = as.integer(K_val),
          alpha     = alpha
        )
        lfmm_rv$result      <- res
        n_sig <- sum(res$significant, na.rm = TRUE)
        gif_vals <- unique(round(res$gif, 2))
        gif_warn <- dplyr::case_when(
          any(gif_vals > 2, na.rm = TRUE) ~ paste0(
            " ⚠️ GIF > 2 detected (",
            paste(gif_vals[gif_vals > 2], collapse = ", "),
            ") — strong undercorrection; increase K or add geographic covariates."
          ),
          any(gif_vals > 1.5, na.rm = TRUE) ~ paste0(
            " ℹ️ GIF 1.5–2 detected (",
            paste(gif_vals[gif_vals > 1.5 & gif_vals <= 2], collapse = ", "),
            ") — mild inflation; check for batch effects or unmodelled structure."
          ),
          any(gif_vals < 0.8, na.rm = TRUE) ~ paste0(
            " ⚠️ GIF < 0.8 detected (",
            paste(gif_vals[gif_vals < 0.8], collapse = ", "),
            ") — over-correction; consider reducing K."
          ),
          TRUE ~ ""
        )
        lfmm_rv$status_type <- if (any(gif_vals > 2 | gif_vals < 0.8, na.rm = TRUE)) "warning" else "success"
        lfmm_rv$status_msg  <- paste0(
          "✅ LFMM2 complete: ", n_sig, " significant loci (BH α = ", alpha, "). ",
          "GIF: ", paste(gif_vals, collapse = ", "), ".", gif_warn
        )
        shinyjs::show("lfmm_plot_panel")
        shinyjs::enable("dl_lfmm2")
      }, error = function(e2) {
        lfmm_rv$status_type <- "error"
        lfmm_rv$status_msg  <- paste0("❌ LFMM2 error: ", e2$message,
                                       " — Is the 'lfmm' package installed?")
        shinyjs::hide("lfmm_plot_panel")
        shinyjs::disable("dl_lfmm2")
      })
    }

    observeEvent(input$run_lfmm2, {
      d <- upload_data(); e <- env_data()
      req(d$ready, e$ready)
      shinyjs::disable("run_lfmm2")
      on.exit(shinyjs::enable("run_lfmm2"), add = TRUE)
      run_lfmm2_analysis(d, e)
    })

    # ── LFMM2 status ─────────────────────────────────────────────────────────
    output$lfmm_status <- renderUI({
      req(!is.null(lfmm_rv$status_msg))
      gea_status_alert(lfmm_rv$status_msg, lfmm_rv$status_type)
    })

    # ── Manhattan plot ────────────────────────────────────────────────────────
    output$manhattan_plot <- plotly::renderPlotly({
      res  <- lfmm_rv$result; req(res)
      d    <- upload_data()
      alpha <- input$lfmm_alpha %||% 0.05

      env_vars <- unique(res$env_var)
      n_vars   <- length(env_vars)
      palette  <- stats::setNames(
        PCA_PALETTE[((seq_along(env_vars) - 1L) %% length(PCA_PALETTE)) + 1L],
        env_vars
      )

      fig <- plotly::plot_ly()

      for (vn in env_vars) {
        sub      <- res[res$env_var == vn, ]
        log_p    <- -log10(pmax(sub$p_adj, 1e-300))
        is_sig   <- sub$significant
        col_fill <- ifelse(is_sig, palette[vn], "#CCCCCC")

        hover_txt <- paste0(
          "Locus ", sub$locus_idx, "<br/>",
          "Env: ", vn, "<br/>",
          "-log10(p_adj) = ", round(log_p, 3), "<br/>",
          "z = ", round(sub$z_score, 3), "<br/>",
          ifelse(is_sig, "<b>SIGNIFICANT</b>", "not significant")
        )

        fig <- fig |>
          plotly::add_markers(
            x         = sub$locus_idx,
            y         = log_p,
            text      = hover_txt,
            hoverinfo = "text",
            marker    = list(
              color   = col_fill,
              size    = ifelse(is_sig, 7, 4),
              opacity = ifelse(is_sig, 0.9, 0.5)
            ),
            name      = vn
          )
      }

      # FDR significance line
      sig_threshold <- -log10(alpha)

      fig |>
        plotly::layout(
          xaxis  = list(title = "Locus index", showgrid = FALSE),
          yaxis  = list(title = "-log10(adj. p-value)", showgrid = TRUE,
                        gridcolor = "#EEEEEE"),
          shapes = list(list(
            type = "line", x0 = 0, x1 = 1, y0 = sig_threshold, y1 = sig_threshold,
            xref = "paper", yref = "y",
            line = list(color = "#D32F2F", width = 1.5, dash = "dash")
          )),
          annotations = list(list(
            x = 1, y = sig_threshold,
            xref = "paper", yref = "y",
            text = paste0("FDR α = ", alpha),
            showarrow = FALSE, xanchor = "right", yanchor = "bottom",
            font = list(size = 10, color = "#D32F2F")
          )),
          showlegend    = n_vars > 1,
          paper_bgcolor = "white",
          plot_bgcolor  = "white",
          font = list(family = "Helvetica Neue, Helvetica, Arial, sans-serif", size = 12)
        ) |>
        plotly::config(
          displayModeBar = TRUE,
          modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
          displaylogo = FALSE,
          toImageButtonOptions = list(
            format   = "png",
            scale    = 2,
            filename = "popgen_manhattan"
          )
        )
    })

    # ── LFMM2 significant loci table ─────────────────────────────────────────
    output$lfmm_table <- DT::renderDataTable({
      res <- lfmm_rv$result; req(res)
      sig <- res[res$significant, , drop = FALSE]
      if (nrow(sig) == 0) {
        # Show the top 100 loci by adjusted p-value rather than all rows —
        # the full table can be millions of cells in large datasets.
        # NOTE: order() returns a plain vector — do NOT use matrix subscript [, ].
        top100 <- res[order(res$p_adj)[seq_len(min(100L, nrow(res)))], , drop = FALSE]
        sig <- top100
      }
      display <- data.frame(
        "Locus"        = sig$locus_idx,
        "Env variable" = sig$env_var,
        # Keep raw p-values — do NOT pre-round here.  Values < 5e-7 would become
        # 0.000000 if rounded to 6dp.  DT::formatSignif(digits=3) renders them
        # correctly as e.g. "2.34e-10" regardless of magnitude.
        "adj. p-value" = sig$p_adj,
        "p-value"      = sig$pvalue,
        "z-score"      = round(sig$z_score, 3),
        "GIF"          = round(sig$gif,     3),
        "Significant"  = sig$significant,
        check.names    = FALSE,
        stringsAsFactors = FALSE
      )
      DT::datatable(
        display,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength = 15,
          dom        = "ftp",
          ordering   = TRUE,
          columnDefs = list(list(className = "dt-center", targets = 0:5))
        )
      ) |>
        DT::formatSignif(columns = c("adj. p-value", "p-value"), digits = 3) |>
        DT::formatStyle(
          "Significant",
          backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("#BBDEFB", "#F5F5F5"))
        )
    })

    # ── LFMM2 table note (shown when no significant loci) ────────────────────
    output$lfmm_table_note <- renderUI({
      res <- lfmm_rv$result; req(res)
      n_sig <- sum(res$significant, na.rm = TRUE)
      if (n_sig == 0) {
        tags$p(
          style = "font-size:0.78em; color:#888; margin-bottom:4px;",
          tags$em(paste0(
            "No loci reached significance at this FDR threshold. ",
            "Showing top 100 loci by adjusted p-value."
          ))
        )
      } else NULL
    })

    # ── LFMM2 export ─────────────────────────────────────────────────────────
    output$dl_lfmm2 <- downloadHandler(
      filename = function() paste0("popgen_lfmm2_", Sys.Date(), ".csv"),
      content  = function(file) {
        res <- lfmm_rv$result; req(res)
        utils::write.csv(res, file, row.names = FALSE)
      }
    )

    # ── Return both results ───────────────────────────────────────────────────
    return(reactive(list(
      rda_result  = rda_rv$result,
      lfmm_result = lfmm_rv$result
    )))
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Summary card for GEA panel
gea_card <- function(label, value, colour, icon_name) {
  div(
    class = "col-6 col-md-3",
    div(
      class = "card border-0 shadow-sm h-100",
      style = paste0("border-left: 4px solid ", colour, " !important;"),
      div(
        class = "card-body p-2 d-flex align-items-center gap-2",
        div(
          style = paste0("color:", colour, "; font-size:1.3em;"),
          icon(icon_name)
        ),
        div(
          tags$div(style = "font-size:1.2em; font-weight:700; line-height:1.2;",
                   as.character(value)),
          tags$div(style = "font-size:0.72em; color:#777; line-height:1.2;",
                   label)
        )
      )
    )
  )
}

#' Small Bootstrap alert for GEA status messages
gea_status_alert <- function(msg, type) {
  cls <- switch(
    type %||% "info",
    "success" = "alert alert-success",
    "warning" = "alert alert-warning",
    "error"   = "alert alert-danger",
    "alert alert-info"
  )
  div(
    class = cls,
    style = "font-size:0.80em; padding:7px 12px; margin-top:6px; margin-bottom:6px;",
    msg
  )
}
