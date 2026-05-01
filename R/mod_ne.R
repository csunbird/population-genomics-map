# =============================================================================
# mod_ne.R — Effective Population Size (Ne) module (Phase 3, Goal 2)
#
# Estimates Ne using the single-sample linkage disequilibrium (LD) method
# of Waples & Do (2008) and displays:
#   1. Summary cards: Ne per population with traffic-light colouring
#   2. Plotly bar chart: Ne ± 95% CI with IUCN threshold lines at 50 and 500
#   3. Results table with export
#
# The map module receives ne_stats() so it can colour markers by Ne.
#
# Design decisions:
#  - Ne computed once per upload (reactive); max_loci slider controls speed
#    vs precision trade-off but is not shown by default (advanced users only)
#  - r² computed from dosage correlations — no external LD packages required
#  - Jackknife CI over loci; large Ne estimates are capped at 10 000
#  - Ne < 1 clamped to 1 (numerical artifact of negative bias-corrected r²)
# =============================================================================

#' Ne module UI
#'
#' @param id Module namespace ID
neUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("load_placeholder")),
    uiOutput(ns("scenario_banner")),
    uiOutput(ns("summary_cards")),

    shinyjs::hidden(
      div(
        id = ns("ne_ready_content"),

        tags$hr(style = "margin: 12px 0;"),

        # ── Controls ───────────────────────────────────────────────────────
        div(
          class = "d-flex align-items-center gap-4 flex-wrap mb-2",
          div(
            style = "flex: 0 0 auto;",
            tags$label("Max loci per population:",
                       class = "me-1 fw-semibold",
                       style = "font-size:0.9em;"),
            sliderInput(
              ns("max_loci"),
              label  = NULL,
              min    = 50L, max = 500L, value = 300L, step = 50L,
              ticks  = FALSE,
              width  = "200px"
            ),
            tags$p(
              style = "font-size:0.74em; color:#888; margin:-4px 0 0;",
              "More loci = more precision, slower computation."
            )
          ),
          div(
            style = "flex: 0 0 auto; margin-top:22px;",
            downloadButton(
              ns("dl_ne"),
              label = tagList(icon("download"), " Export CSV"),
              class = "btn btn-outline-secondary btn-sm"
            )
          )
        ),

        # ── Bar chart ──────────────────────────────────────────────────────
        plotly::plotlyOutput(ns("ne_plot"), height = "420px"),

        # ── Results table ──────────────────────────────────────────────────
        tags$h6(
          "Per-population Ne estimates",
          style = "font-weight:600; color:#1B3A5C; margin: 14px 0 6px;"
        ),
        DT::dataTableOutput(ns("ne_table"))
      )
    ),

    # ── Footnotes ─────────────────────────────────────────────────────────────
    tags$p(
      class = "mt-2",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "Ne estimated from the mean r&sup2; of all pairwise locus combinations ",
        "within each population (Waples &amp; Do 2008, single-sample LD method). ",
        "Sample-size bias is corrected using the Waples (2006) formula. ",
        "95% confidence intervals are computed by jackknife over loci. ",
        "Ne estimates are capped at 10,000 (effectively panmictic)."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML(paste0(
        "&#9432; IUCN/SSC thresholds: Ne &lt; 50 = immediate extinction risk; ",
        "Ne &lt; 100 = high risk; Ne &lt; 500 = long-term viability concern ",
        "(Franklin 1980 / IUCN SSC 2023). ",
        "The LD method assumes random mating (Hardy-Weinberg) within populations. ",
        "Violations (selfing, overlapping generations, structured subpopulations) ",
        "can bias Ne downward. ",
        "Small sample sizes (n &lt; 20) produce wide confidence intervals — ",
        "treat such estimates as indicative only. ",
        "Without a physical map, physically linked loci cannot be excluded; ",
        "this inflates r&sup2; and leads to underestimation of Ne."
      ))
    )
  )
}

#' Ne module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive returning the upload module's data list
#'                    ($ready, $gt, $metadata, $is_demo, $scenario_info)
#' @return Reactive returning per-population Ne data frame (or NULL)
#'         — consumed by mapServer to colour markers by Ne
neServer <- function(id, upload_data) {
  moduleServer(id, function(input, output, session) {

    # ── Load placeholder ─────────────────────────────────────────────────────
    output$load_placeholder <- renderUI({
      d <- upload_data()
      if (isTRUE(d$ready)) return(NULL)
      div(
        class = "alert alert-secondary",
        style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
        icon("info-circle"), " ",
        "Load a VCF file and metadata CSV from the sidebar to estimate Ne."
      )
    })

    # ── Show/hide controls ───────────────────────────────────────────────────
    observeEvent(upload_data(), {
      shinyjs::toggle("ne_ready_content", condition = isTRUE(upload_data()$ready))
    })

    # ── Compute Ne (reactive — reruns when data OR max_loci changes) ─────────
    ne_result <- reactive({
      d <- upload_data()
      req(d$ready)
      max_l <- input$max_loci %||% 300L
      tryCatch(
        calc_ne_ld(d$gt, d$metadata, max_loci = as.integer(max_l)),
        error = function(e) {
          showNotification(
            paste0("Ne estimation error: ", e$message),
            type = "error", duration = 10
          )
          NULL
        }
      )
    })

    # ── Scenario banner ──────────────────────────────────────────────────────
    output$scenario_banner <- renderUI({
      d <- upload_data()
      req(d$ready, d$scenario_info)
      cls <- if (isTRUE(d$is_demo)) "alert alert-warning" else "alert alert-info"
      div(
        class = cls,
        style = "font-size:0.82em; padding:8px 12px; margin-bottom:10px;",
        icon("info-circle"), " ", d$scenario_info
      )
    })

    # ── Summary cards ────────────────────────────────────────────────────────
    output$summary_cards <- renderUI({
      df <- ne_result()
      req(df)

      valid_ne    <- df$ne[!is.na(df$ne)]
      min_ne      <- if (length(valid_ne) > 0) round(min(valid_ne), 0) else NA
      min_pop     <- if (length(valid_ne) > 0) df$population[which.min(df$ne)] else "—"
      n_critical  <- sum(!is.na(df$ne) & df$ne < NE_THRESHOLDS$critical)
      n_concern   <- sum(!is.na(df$ne) & df$ne < NE_THRESHOLDS$concern)
      crit_col    <- if (n_critical > 0) "#D32F2F" else "#388E3C"
      conc_col    <- if (n_concern  > 0) "#F57C00" else "#388E3C"

      # Truncate long population name for card display
      min_pop_short <- if (nchar(min_pop) > 18) paste0(substr(min_pop, 1, 15), "…") else min_pop

      div(
        class = "row g-2 mb-2",
        ne_card("Lowest Ne",           if (!is.na(min_ne)) as.character(min_ne) else "—",
                "#D32F2F", "exclamation-triangle"),
        ne_card("Lowest Ne population", min_pop_short,    "#984EA3", "map-marker-alt"),
        ne_card("Populations Ne < 50",  n_critical,       crit_col,  "skull-crossbones"),
        ne_card("Populations Ne < 500", n_concern,        conc_col,  "exclamation-circle")
      )
    })

    # ── Bar chart ────────────────────────────────────────────────────────────
    output$ne_plot <- plotly::renderPlotly({
      df <- ne_result()
      req(df)

      # Colour bars by Ne conservation threshold
      bar_colours <- vapply(df$ne, ne_colour, character(1))

      # Cap display at 5000 for readability.
      # Bars that are capped receive a "↑" bar label so users know the true Ne
      # is higher than the bar height implies.
      NE_DISPLAY_CAP <- 5000
      is_capped      <- !is.na(df$ne) & df$ne > NE_DISPLAY_CAP
      ne_display     <- pmin(df$ne,       NE_DISPLAY_CAP)
      ne_lo_display  <- pmin(df$ne_lower, NE_DISPLAY_CAP)
      ne_hi_display  <- pmin(df$ne_upper, NE_DISPLAY_CAP)

      # Bar-top labels: show actual Ne for capped bars, empty for others
      bar_labels <- ifelse(
        is_capped,
        paste0("Ne = ", format(round(df$ne), big.mark = ",")),
        ""
      )

      hover_txt <- paste0(
        "<b>", htmltools::htmlEscape(df$population), "</b><br/>",
        "Ne = ", ifelse(is.na(df$ne), "—", format(round(df$ne), big.mark = ",")),
        ifelse(!is.na(df$ne_lower),
               paste0(" (95% CI: ", format(round(df$ne_lower), big.mark = ","), "–",
                      format(round(df$ne_upper), big.mark = ","), ")"),
               ""),
        ifelse(is_capped, "<br/><i>(bar capped at 5 000 for display)</i>", ""),
        "<br/>", ne_risk_label(df$ne),
        "<br/>n = ", df$n_samples, " | loci = ", df$n_loci_used
      )

      fig <- plotly::plot_ly(
        x          = df$population,
        y          = ne_display,
        type       = "bar",
        marker     = list(color = bar_colours, line = list(width = 0)),
        error_y    = if (!all(is.na(ne_lo_display))) {
          list(
            type       = "data",
            symmetric  = FALSE,
            array      = ne_hi_display - ne_display,
            arrayminus = ne_display - ne_lo_display,
            color      = "#555555",
            thickness  = 1.5,
            width      = 6
          )
        } else NULL,
        text       = bar_labels,
        textposition = "outside",
        hovertext  = hover_txt,
        hoverinfo  = "text"
      )

      # IUCN threshold lines
      threshold_shapes <- list(
        ne_hline(NE_THRESHOLDS$critical,   "#D32F2F", "dash"),
        ne_hline(NE_THRESHOLDS$vulnerable, "#F57C00", "dash"),
        ne_hline(NE_THRESHOLDS$concern,    "#FBC02D", "dot")
      )

      # Threshold annotations (xref = "paper" so x = 1 always means the right edge)
      y_max <- max(ne_display, NE_THRESHOLDS$concern * 1.1, na.rm = TRUE)

      threshold_annotations <- list(
        ne_ann(NE_THRESHOLDS$critical,   "Ne = 50 (critical)",  "#D32F2F"),
        ne_ann(NE_THRESHOLDS$vulnerable, "Ne = 100 (high risk)","#F57C00"),
        ne_ann(NE_THRESHOLDS$concern,    "Ne = 500 (concern)",  "#8B6914")
      )

      fig |>
        plotly::layout(
          xaxis  = list(title = "", showgrid = FALSE, tickangle = -30),
          yaxis  = list(
            title    = "Effective population size (Ne)",
            range    = c(0, y_max * 1.1),
            showgrid = TRUE,
            gridcolor = "#EEEEEE"
          ),
          shapes      = threshold_shapes,
          annotations = threshold_annotations,
          paper_bgcolor = "white",
          plot_bgcolor  = "white",
          showlegend  = FALSE,
          font = list(
            family = "Helvetica Neue, Helvetica, Arial, sans-serif",
            size   = 12
          )
        ) |>
        plotly::config(
          displayModeBar         = TRUE,
          modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
          displaylogo            = FALSE,
          toImageButtonOptions   = list(
            format   = "png",
            scale    = 2,
            filename = "popgen_ne"
          )
        )
    })

    # ── Results table ────────────────────────────────────────────────────────
    output$ne_table <- DT::renderDataTable({
      df <- ne_result()
      req(df)
      # Keep numeric columns as numeric so DT sorts them correctly.
      # NA values appear as blank cells in DT (standard behaviour).
      display <- data.frame(
        "Population"   = df$population,
        "Ne"           = round(df$ne),
        "95% CI lower" = round(df$ne_lower),
        "95% CI upper" = round(df$ne_upper),
        "n samples"    = df$n_samples,
        "Loci used"    = df$n_loci_used,
        "Mean r²"      = round(df$mean_r2, 5),
        "Risk"         = ne_risk_label(df$ne),
        check.names    = FALSE,
        stringsAsFactors = FALSE
      )
      DT::datatable(
        display,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength = 10,
          dom        = "t",
          ordering   = TRUE,
          columnDefs = list(list(className = "dt-center", targets = 1:6))
        )
      ) |>
        DT::formatRound(columns = c("Ne", "95% CI lower", "95% CI upper"), digits = 0) |>
        DT::formatRound(columns = "Mean r²", digits = 5) |>
        DT::formatStyle(
          "Risk",
          backgroundColor = DT::styleEqual(
            c("Critical — immediate extinction risk",
              "High risk — rapid genetic drift",
              "Concern — long-term viability at risk",
              "Viable — Ne ≥ 500",
              "Unknown"),
            c("#FFCDD2", "#FFE0B2", "#FFF9C4", "#BBDEFB", "#F5F5F5")
          )
        )
    })

    # ── CSV download ─────────────────────────────────────────────────────────
    output$dl_ne <- downloadHandler(
      filename = function() paste0("popgen_ne_", Sys.Date(), ".csv"),
      content  = function(file) {
        df <- ne_result()
        req(df)
        df$risk <- ne_risk_label(df$ne)
        write.csv(df, file, row.names = FALSE)
      }
    )

    # ── Return ne_result for map overlay ─────────────────────────────────────
    ne_result
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Individual summary card for the Ne panel
ne_card <- function(label, value, colour, icon_name) {
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
          tags$div(
            style = "font-size:1.2em; font-weight:700; line-height:1.2;",
            as.character(value)
          ),
          tags$div(
            style = "font-size:0.72em; color:#777; line-height:1.2;",
            label
          )
        )
      )
    )
  )
}

#' Plotly horizontal threshold line for Ne chart
ne_hline <- function(y, colour, dash = "dash") {
  list(
    type  = "line",
    x0    = 0, x1    = 1,
    y0    = y, y1    = y,
    xref  = "paper", yref = "y",
    line  = list(color = colour, width = 1.5, dash = dash)
  )
}

#' Plotly annotation label for Ne threshold line
ne_ann <- function(y, label, colour) {
  list(
    x         = 1,
    y         = y,
    xref      = "paper",
    yref      = "y",
    text      = label,
    showarrow = FALSE,
    xanchor   = "right",
    yanchor   = "bottom",
    font      = list(size = 10, color = colour)
  )
}
