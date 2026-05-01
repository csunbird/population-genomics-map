# =============================================================================
# mod_offset.R — Genomic Offset / Climate Vulnerability module (Phase 4, Goal 3)
#
# Computes population-level genomic offset between current and projected future
# climate conditions using the RDA-based method of Fitzpatrick & Keller (2015).
#
# Genomic offset = Euclidean distance in RDA-constrained space between a
# population's current and projected future environmental position. Larger
# values indicate a wider "adaptive gap" — the population must evolve further
# to remain locally adapted under the projected future climate.
#
# Displays:
#   1. Summary cards: mean / max offset, highest-risk population, n at risk
#   2. Leaflet map: populations coloured + sized by offset, hover details
#   3. Bar chart: normalised offset per population with risk colour bands
#   4. Results table with export
#
# Dependencies:
#   - rda_result from geaServer (must be run before offset can be computed)
#   - future_env_data from envUploadServer
# =============================================================================

#' Genomic offset module UI
#'
#' @param id Module namespace ID
offsetUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("load_placeholder")),
    uiOutput(ns("scenario_banner")),
    uiOutput(ns("summary_cards")),

    shinyjs::hidden(
      div(
        id = ns("offset_ready_content"),

        tags$hr(style = "margin: 12px 0;"),

        # ── Controls ────────────────────────────────────────────────────────
        div(
          class = "d-flex align-items-center gap-3 flex-wrap mb-2",
          div(
            actionButton(
              ns("run_offset"),
              label = tagList(icon("play"), " Compute Genomic Offset"),
              class = "btn btn-primary btn-sm"
            )
          ),
          div(
            shinyjs::disabled(
              downloadButton(
                ns("dl_offset"),
                label = tagList(icon("download"), " Export CSV"),
                class = "btn btn-outline-secondary btn-sm"
              )
            )
          )
        ),

        uiOutput(ns("offset_status")),

        # ── Map + bar chart in two columns ────────────────────────────────
        shinyjs::hidden(
          div(
            id = ns("offset_results_panel"),

            div(
              class = "row g-3 mt-1",

              # Map column
              div(
                class = "col-12 col-lg-6",
                tags$h6(
                  "Spatial distribution of genomic offset",
                  style = "font-weight:600; color:#1B3A5C; margin-bottom:6px;"
                ),
                leaflet::leafletOutput(ns("offset_map"), height = "380px")
              ),

              # Bar chart column
              div(
                class = "col-12 col-lg-6",
                tags$h6(
                  "Normalised genomic offset by population",
                  style = "font-weight:600; color:#1B3A5C; margin-bottom:6px;"
                ),
                plotly::plotlyOutput(ns("offset_bar"), height = "380px")
              )
            ),

            # Results table
            tags$h6(
              "Per-population genomic offset",
              style = "font-weight:600; color:#1B3A5C; margin: 16px 0 6px;"
            ),
            DT::dataTableOutput(ns("offset_table"))
          )
        )
      )
    ),

    # ── Footnotes ──────────────────────────────────────────────────────────────
    tags$p(
      class = "mt-3",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "Genomic offset is computed using the RDA-based method of ",
        "Fitzpatrick &amp; Keller (2015, <i>Ecology Letters</i>). ",
        "Current and future environmental conditions are projected through the ",
        "RDA biplot loadings into the constrained (genotype-environment) ordination space. ",
        "Offset = Euclidean distance between current and future projections. ",
        "Normalised offset (0–1) is scaled relative to the maximum observed offset ",
        "across all populations in this dataset."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML(paste0(
        "&#9432; Genomic offset reflects the magnitude of climate-driven selection pressure, ",
        "not the probability of extinction — populations with high adaptive potential, ",
        "high connectivity, or large Ne may tolerate larger offsets. ",
        "Results depend on the quality and representativeness of the environmental predictors ",
        "and on the RDA model. Always interpret in conjunction with Ne, F-ROH, and FST results. ",
        "Climate projections carry inherent uncertainty; ensemble models are recommended ",
        "for conservation decision-making."
      ))
    )
  )
}

#' Genomic offset module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive: upload module data list
#' @param env_data    Reactive: env upload module data list ($env_data, $future_env_data)
#' @param gea_data    Reactive: GEA module data list ($rda_result)
#' @return Reactive returning per-population offset data frame (or NULL)
offsetServer <- function(id, upload_data, env_data, gea_data) {
  moduleServer(id, function(input, output, session) {

    rv <- reactiveValues(
      result      = NULL,
      status_type = NULL,
      status_msg  = NULL
    )

    # ── Ready condition ──────────────────────────────────────────────────────
    all_ready <- reactive({
      d <- upload_data(); e <- env_data(); g <- gea_data()
      isTRUE(d$ready) && isTRUE(e$ready) && !is.null(g$rda_result) && !is.null(e$future_env_data)
    })

    inputs_partial <- reactive({
      d <- upload_data(); e <- env_data()
      isTRUE(d$ready) && isTRUE(e$ready)
    })

    # ── Clear stale results when genomic data resets ─────────────────────────
    # Without this, summary cards outside the hidden div would show results
    # from a previous dataset until the user re-runs offset computation.
    observeEvent(upload_data(), {
      if (!isTRUE(upload_data()$ready)) {
        rv$result      <- NULL
        rv$status_type <- NULL
        rv$status_msg  <- NULL
        shinyjs::hide("offset_results_panel")
        shinyjs::disable("dl_offset")
      }
    })

    # ── Load placeholder ─────────────────────────────────────────────────────
    output$load_placeholder <- renderUI({
      d <- upload_data(); e <- env_data(); g <- gea_data()

      if (!isTRUE(d$ready)) {
        return(div(
          class = "alert alert-secondary",
          style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
          icon("info-circle"), " ",
          "Load genomic data from the sidebar to enable genomic offset analysis."
        ))
      }
      if (!isTRUE(e$ready)) {
        return(div(
          class = "alert alert-info",
          style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
          icon("info-circle"), " ",
          "Load environmental data (sidebar, Phase 4 section) to enable offset analysis."
        ))
      }
      if (is.null(g$rda_result)) {
        return(div(
          class = "alert alert-info",
          style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
          icon("info-circle"), " ",
          "Run the Partial RDA on the GEA tab first — genomic offset is computed ",
          "from the RDA constrained space."
        ))
      }
      if (is.null(e$future_env_data)) {
        return(div(
          class = "alert alert-warning",
          style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
          icon("exclamation-triangle"), " ",
          "No future climate data loaded. Upload a future env CSV or fetch CMIP6 ",
          "projections from the sidebar (Phase 4 section) to compute offset."
        ))
      }
      NULL
    })

    # ── Show/hide controls ───────────────────────────────────────────────────
    observe({
      shinyjs::toggle("offset_ready_content", condition = isTRUE(inputs_partial()))
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

    # ── Auto-run on demo data ────────────────────────────────────────────────
    observeEvent(list(upload_data(), env_data(), gea_data()), {
      d <- upload_data(); g <- gea_data(); e <- env_data()
      if (!isTRUE(d$is_demo)) return()
      if (is.null(g$rda_result) || is.null(e$future_env_data)) return()
      if (!is.null(rv$result)) return()
      run_offset_calc(g$rda_result, e$future_env_data)
    })

    # ── Compute offset ───────────────────────────────────────────────────────
    run_offset_calc <- function(rda_result, future_env) {
      rv$result      <- NULL
      rv$status_type <- "info"
      rv$status_msg  <- "⏳ Computing genomic offset…"

      tryCatch({
        df <- calc_genomic_offset(rda_result, future_env)
        rv$result      <- df
        n_severe <- sum(!is.na(df$offset_norm) & df$offset_norm >= 0.75)
        # n_high: strictly the High band (0.50–0.75), not including Severe
        n_high   <- sum(!is.na(df$offset_norm) & df$offset_norm >= 0.50 & df$offset_norm < 0.75)
        rv$status_type <- if (n_severe > 0) "warning" else "success"
        rv$status_msg  <- paste0(
          "✅ Offset computed: ", nrow(df), " populations | ",
          "Mean offset = ", round(mean(df$offset_norm, na.rm = TRUE), 3), " | ",
          if (n_severe > 0) paste0(n_severe, " population(s) at severe risk. ")
          else if (n_high > 0) paste0(n_high, " population(s) at high risk. ")
          else "No populations at high or severe risk."
        )
        shinyjs::show("offset_results_panel")
        shinyjs::enable("dl_offset")
      }, error = function(e2) {
        rv$status_type <- "error"
        rv$status_msg  <- paste0("❌ Offset error: ", e2$message)
        shinyjs::hide("offset_results_panel")
        shinyjs::disable("dl_offset")
      })
    }

    observeEvent(input$run_offset, {
      g <- gea_data(); e <- env_data()
      if (is.null(g$rda_result)) {
        rv$status_type <- "error"
        rv$status_msg  <- paste0(
          "❌ Run the Partial RDA on the GEA tab first — ",
          "genomic offset requires a fitted RDA model."
        )
        return()
      }
      if (is.null(e$future_env_data)) {
        rv$status_type <- "error"
        rv$status_msg  <- paste0(
          "❌ No future climate data loaded. ",
          "Upload a future env CSV or fetch CMIP6 from the sidebar."
        )
        return()
      }
      shinyjs::disable("run_offset")
      on.exit(shinyjs::enable("run_offset"), add = TRUE)
      run_offset_calc(g$rda_result, e$future_env_data)
    })

    # ── Status display ───────────────────────────────────────────────────────
    output$offset_status <- renderUI({
      req(!is.null(rv$status_msg))
      cls <- switch(
        rv$status_type %||% "info",
        "success" = "alert alert-success",
        "warning" = "alert alert-warning",
        "error"   = "alert alert-danger",
        "alert alert-info"
      )
      div(
        class = cls,
        style = "font-size:0.80em; padding:7px 12px; margin-top:6px; margin-bottom:6px;",
        rv$status_msg
      )
    })

    # ── Summary cards ────────────────────────────────────────────────────────
    output$summary_cards <- renderUI({
      df <- rv$result; req(df)

      valid     <- df$offset_norm[!is.na(df$offset_norm)]
      mean_off  <- if (length(valid) > 0) round(mean(valid), 3) else NA
      max_off   <- if (length(valid) > 0) round(max(valid), 3) else NA
      worst_pop <- if (length(valid) > 0) df$population[which.max(df$offset_norm)] else "—"
      n_severe  <- sum(!is.na(df$offset_norm) & df$offset_norm >= 0.75)
      sev_col   <- if (n_severe > 0) "#D32F2F" else "#388E3C"

      # Truncate long pop name
      worst_short <- if (nchar(worst_pop) > 18) paste0(substr(worst_pop, 1, 15), "…")
                     else worst_pop

      div(
        class = "row g-2 mb-2",
        offset_card("Mean norm. offset",  if (!is.na(mean_off)) as.character(mean_off) else "—",
                    "#1B3A5C", "chart-line"),
        offset_card("Max norm. offset",   if (!is.na(max_off)) as.character(max_off) else "—",
                    "#D32F2F", "exclamation-triangle"),
        offset_card("Highest-risk population", worst_short,
                    "#984EA3", "map-marker-alt"),
        offset_card("Populations severe risk", n_severe,
                    sev_col, "skull-crossbones")
      )
    })

    # ── Leaflet map ──────────────────────────────────────────────────────────
    output$offset_map <- leaflet::renderLeaflet({
      df <- rv$result; req(df)

      # Circle radius: scale by offset (5–20 px).
      # Replace NA offset_norm with 0 so Leaflet doesn't crash on missing values;
      # such populations still appear on the map at minimum size.
      safe_norm <- ifelse(is.na(df$offset_norm), 0, df$offset_norm)
      r_scaled  <- 5 + 15 * safe_norm

      leaflet::leaflet(df) |>
        leaflet::addProviderTiles(
          leaflet::providers$CartoDB.Positron,
          options = leaflet::providerTileOptions(opacity = 0.85)
        ) |>
        leaflet::addCircleMarkers(
          lng         = ~lon,
          lat         = ~lat,
          radius      = r_scaled,
          color       = ~risk_colour,
          fillColor   = ~risk_colour,
          fillOpacity = 0.75,
          weight      = 2,
          opacity     = 0.9,
          popup = ~paste0(
            "<b>", htmltools::htmlEscape(population), "</b><br/>",
            "Genomic offset: ", round(offset, 4), "<br/>",
            "Normalised: ", round(offset_norm, 3), "<br/>",
            "<span style='color:", risk_colour, ";'>", risk_label, "</span>"
          ),
          label = ~htmltools::htmlEscape(population)
        ) |>
        leaflet::addLegend(
          position = "bottomright",
          colors   = c("#D32F2F", "#F57C00", "#FBC02D", "#1565C0"),
          labels   = c("Severe (≥ 0.75)", "High (≥ 0.50)",
                       "Moderate (≥ 0.25)", "Low (< 0.25)"),
          title    = "Genomic offset",
          opacity  = 0.85
        )
    })

    # ── Bar chart ────────────────────────────────────────────────────────────
    output$offset_bar <- plotly::renderPlotly({
      df <- rv$result; req(df)

      # Sort by offset descending for readability
      df <- df[order(df$offset_norm, decreasing = TRUE), ]

      hover_txt <- paste0(
        "<b>", htmltools::htmlEscape(df$population), "</b><br/>",
        "Normalised offset: ", round(df$offset_norm, 3), "<br/>",
        "Raw offset: ", round(df$offset, 4), "<br/>",
        df$risk_label
      )

      fig <- plotly::plot_ly(
        x         = df$population,
        y         = df$offset_norm,
        type      = "bar",
        marker    = list(color = df$risk_colour, line = list(width = 0)),
        text      = round(df$offset_norm, 3),
        textposition = "outside",
        hovertext = hover_txt,
        hoverinfo = "text"
      )

      # Risk band shading
      band_shapes <- list(
        offset_band(0.75, 1.05, "#FFCDD2", 0.25),   # severe
        offset_band(0.50, 0.75, "#FFE0B2", 0.20),   # high
        offset_band(0.25, 0.50, "#FFF9C4", 0.15),   # moderate
        offset_band(0.00, 0.25, "#BBDEFB", 0.10)    # low
      )

      band_annotations <- list(
        offset_band_ann(0.875, "Severe",   "#D32F2F"),
        offset_band_ann(0.625, "High",     "#F57C00"),
        offset_band_ann(0.375, "Moderate", "#8B6914"),
        offset_band_ann(0.125, "Low",      "#1565C0")
      )

      fig |>
        plotly::layout(
          xaxis  = list(title = "", showgrid = FALSE, tickangle = -30),
          yaxis  = list(
            title     = "Normalised genomic offset",
            range     = c(0, min(1.12, max(df$offset_norm, na.rm = TRUE) * 1.15 + 0.05)),
            showgrid  = TRUE,
            gridcolor = "#EEEEEE"
          ),
          shapes      = band_shapes,
          annotations = band_annotations,
          paper_bgcolor = "white",
          plot_bgcolor  = "white",
          showlegend    = FALSE,
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
            filename = "popgen_genomic_offset"
          )
        )
    })

    # ── Results table ─────────────────────────────────────────────────────────
    output$offset_table <- DT::renderDataTable({
      df <- rv$result; req(df)

      display <- data.frame(
        "Population"       = df$population,
        "Offset (raw)"     = round(df$offset,      4),
        "Offset (norm.)"   = round(df$offset_norm, 3),
        "Latitude"         = round(df$lat,          4),
        "Longitude"        = round(df$lon,          4),
        "Risk category"    = df$risk_label,
        check.names        = FALSE,
        stringsAsFactors   = FALSE
      )
      # Sort by offset descending
      display <- display[order(display[["Offset (norm.)"]], decreasing = TRUE), ]

      DT::datatable(
        display,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength = 10,
          dom        = "t",
          ordering   = TRUE,
          columnDefs = list(list(className = "dt-center", targets = 1:4))
        )
      ) |>
        DT::formatRound(columns = c("Offset (raw)", "Offset (norm.)"), digits = 4) |>
        DT::formatStyle(
          "Risk category",
          backgroundColor = DT::styleEqual(
            c("Severe — highest climate exposure",
              "High — significant climate exposure",
              "Moderate — some climate exposure",
              "Low — least climate exposure",
              "Unknown"),
            c("#FFCDD2", "#FFE0B2", "#FFF9C4", "#BBDEFB", "#F5F5F5")
          )
        )
    })

    # ── CSV download ─────────────────────────────────────────────────────────
    output$dl_offset <- downloadHandler(
      filename = function() paste0("popgen_offset_", Sys.Date(), ".csv"),
      content  = function(file) {
        df <- rv$result; req(df)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    # ── Return result for external use (e.g. map overlay) ────────────────────
    reactive(rv$result)
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Summary card for the offset panel
offset_card <- function(label, value, colour, icon_name) {
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

#' Background band shape for the offset bar chart
#'
#' @param y0,y1     Y-axis extent of the band
#' @param fill      Fill hex colour
#' @param opacity   Fill opacity
offset_band <- function(y0, y1, fill, opacity = 0.15) {
  list(
    type    = "rect",
    x0      = 0, x1      = 1,
    y0      = y0, y1     = y1,
    xref    = "paper", yref = "y",
    fillcolor = fill,
    opacity   = opacity,
    line      = list(width = 0)
  )
}

#' Annotation label for an offset risk band
offset_band_ann <- function(y, label, colour) {
  list(
    x         = 1,
    y         = y,
    xref      = "paper",
    yref      = "y",
    text      = label,
    showarrow = FALSE,
    xanchor   = "right",
    yanchor   = "middle",
    font      = list(size = 9, color = colour)
  )
}
