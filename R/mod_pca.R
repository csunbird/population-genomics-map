# =============================================================================
# mod_pca.R — PCA scatter plot module (Phase 2, Goal 3)
#
# Displays an interactive Plotly scatter plot of the first N principal
# components, coloured by population. Individuals are plotted as points;
# hovering shows sample ID and population. Variance explained per PC is
# shown on axis labels and in summary cards.
#
# Design decisions:
#  - PCA is computed once per upload (reactive), then re-plotted when PC
#    axis selectors change (no recomputation)
#  - Missing genotypes are mean-imputed per locus before prcomp()
#  - No genotype scaling (scale. = FALSE) — dosage values (0/1/2) share
#    a common unit and scaling would over-weight rare variants
#  - Population colour palette is consistent across PC pairs
# =============================================================================

# Colour palette for up to 20 populations (adapted from ColorBrewer Set1 + Set2)
PCA_PALETTE <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A"
)

#' PCA module UI
#'
#' @param id Module namespace ID
pcaUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("load_placeholder")),
    uiOutput(ns("scenario_banner")),
    uiOutput(ns("summary_cards")),

    # ── Axis selectors, plot, and divider ────────────────────────────────────
    # Hidden via shinyjs until data are ready (toggled server-side).
    # Using shinyjs::hidden rather than a renderUI gate keeps the selectInputs
    # in the DOM at all times so updateSelectInput() works regardless of timing.
    shinyjs::hidden(
      div(
        id = ns("pca_ready_content"),

        tags$hr(style = "margin: 12px 0;"),

        div(
          class = "d-flex align-items-center gap-3 flex-wrap mb-2",
          div(
            style = "flex: 0 0 auto;",
            tags$label("X axis:", class = "me-1 fw-semibold",
                       style = "font-size:0.9em;"),
            selectInput(ns("pc_x"), label = NULL,
                        choices = paste0("PC", 1:10), selected = "PC1",
                        width = "100px")
          ),
          div(
            style = "flex: 0 0 auto;",
            tags$label("Y axis:", class = "me-1 fw-semibold",
                       style = "font-size:0.9em;"),
            selectInput(ns("pc_y"), label = NULL,
                        choices = paste0("PC", 1:10), selected = "PC2",
                        width = "100px")
          ),
          div(
            style = "flex: 0 0 auto; margin-top:22px;",
            downloadButton(
              ns("dl_scores"),
              label = tagList(icon("download"), " Export scores"),
              class = "btn btn-outline-secondary btn-sm"
            )
          )
        ),

        plotly::plotlyOutput(ns("pca_plot"), height = "480px")
      )
    ),

    # ── Footnotes ───────────────────────────────────────────────────────────
    tags$p(
      class = "mt-2",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "PCA performed on the full genotype matrix (samples &times; loci). ",
        "Missing genotypes are imputed with the per-locus mean allele dosage ",
        "before computing principal components. ",
        "No scaling applied (all loci are on the 0/1/2 dosage scale). ",
        "Invariant loci are removed before decomposition."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML(paste0(
        "&#9432; PCA of genotype data reflects relatedness and population structure, ",
        "not environmental gradients. Tight clusters indicate within-population similarity; ",
        "overlap may indicate recent gene flow or shared ancestry. ",
        "In datasets with strong linkage disequilibrium (e.g. chromosomal inversions), ",
        "the first PCs may reflect LD blocks rather than population divergence — ",
        "LD-pruned data can give cleaner structure in those cases. ",
        "Individuals with high missing-data rates are pulled toward the population centroid after imputation."
      ))
    )
  )
}

#' PCA module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive returning the upload module's data list
#'                    (fields: $ready, $gt, $metadata, $is_demo, $scenario_info)
pcaServer <- function(id, upload_data) {
  moduleServer(id, function(input, output, session) {

    # ── Load placeholder (shown until data are ready) ────────────────────────
    output$load_placeholder <- renderUI({
      d <- upload_data()
      if (isTRUE(d$ready)) return(NULL)
      div(
        class = "alert alert-secondary",
        style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
        icon("info-circle"), " ",
        "Load a VCF file and metadata CSV from the sidebar to compute PCA."
      )
    })

    # ── Toggle controls + plot visibility when data loads ───────────────────
    observeEvent(upload_data(), {
      shinyjs::toggle("pca_ready_content", condition = isTRUE(upload_data()$ready))
    })

    # ── Compute PCA ──────────────────────────────────────────────────────────
    # Runs once per upload — axis changes only re-plot, not recompute.
    pca_result <- reactive({
      d <- upload_data()
      req(d$ready)
      tryCatch(
        calc_pca(d$gt, d$metadata),
        error = function(e) {
          showNotification(
            paste0("PCA error: ", e$message),
            type     = "error",
            duration = 10
          )
          NULL
        }
      )
    })

    # ── Update axis selectors based on available PCs ─────────────────────────
    observeEvent(pca_result(), {
      res <- pca_result()
      req(res)
      pc_choices <- paste0("PC", seq_len(res$n_pcs))
      pct        <- round(res$var_explained * 100, 1)
      pc_labels  <- paste0("PC", seq_len(res$n_pcs),
                            " (", pct, "%)")
      names(pc_choices) <- pc_labels

      updateSelectInput(session, "pc_x",
                        choices  = pc_choices,
                        selected = pc_choices[1])
      updateSelectInput(session, "pc_y",
                        choices  = pc_choices,
                        selected = if (length(pc_choices) >= 2) pc_choices[2]
                                   else pc_choices[1])
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
      res <- pca_result()
      req(res)

      scores  <- res$scores
      n_inds  <- nrow(scores)
      n_pops  <- length(unique(scores$population))
      pc1_pct <- round(res$var_explained[1] * 100, 1)
      pc2_pct <- if (length(res$var_explained) >= 2)
                   round(res$var_explained[2] * 100, 1) else NA

      div(
        class = "row g-2 mb-2",
        pca_card("Individuals",   n_inds,              "#1B3A5C", "users"),
        pca_card("Populations",   n_pops,              "#2E75B6", "layer-group"),
        pca_card("PC1 variance",  paste0(pc1_pct, "%"), "#388E3C", "chart-bar"),
        pca_card("PC2 variance",
                 if (!is.na(pc2_pct)) paste0(pc2_pct, "%") else "—",
                 "#2E75B6", "chart-bar")
      )
    })

    # ── PCA scatter plot ─────────────────────────────────────────────────────
    output$pca_plot <- plotly::renderPlotly({
      res <- pca_result()
      req(res)

      scores <- res$scores
      pc_x   <- input$pc_x %||% "PC1"
      pc_y   <- input$pc_y %||% "PC2"

      # Guard: selected PCs must exist in the scores data frame
      if (!pc_x %in% names(scores) || !pc_y %in% names(scores)) return(NULL)

      # Guard: same PC selected for both axes → meaningless diagonal plot
      if (pc_x == pc_y) {
        showNotification(
          "X and Y axes show the same PC. Select two different PCs.",
          type = "warning", duration = 6
        )
        return(NULL)
      }

      pct_x  <- round(res$var_explained[as.integer(sub("PC", "", pc_x))] * 100, 1)
      pct_y  <- round(res$var_explained[as.integer(sub("PC", "", pc_y))] * 100, 1)

      pops        <- unique(scores$population)
      pop_colours <- stats::setNames(
        PCA_PALETTE[((seq_along(pops) - 1L) %% length(PCA_PALETTE)) + 1L],
        pops
      )

      # Pre-compute hover text outside the formula — avoids fragile get() lookup
      hover_text <- paste0(
        "<b>", htmltools::htmlEscape(scores$sample_id), "</b><br/>",
        "Population: ", htmltools::htmlEscape(scores$population), "<br/>",
        pc_x, ": ", round(scores[[pc_x]], 3), "<br/>",
        pc_y, ": ", round(scores[[pc_y]], 3)
      )

      plotly::plot_ly(
        data   = scores,
        x      = stats::as.formula(paste0("~`", pc_x, "`")),
        y      = stats::as.formula(paste0("~`", pc_y, "`")),
        color  = ~population,
        colors = pop_colours,
        text      = hover_text,
        hoverinfo = "text",
        type      = "scatter",
        mode      = "markers",
        marker    = list(size = 9, opacity = 0.8, line = list(width = 0.5, color = "white"))
      ) |>
        plotly::layout(
          xaxis = list(
            title      = paste0(pc_x, "  (", pct_x, "% variance)"),
            zeroline   = TRUE,
            zerolinecolor = "#CCCCCC",
            showgrid   = TRUE,
            gridcolor  = "#EEEEEE"
          ),
          yaxis = list(
            title      = paste0(pc_y, "  (", pct_y, "% variance)"),
            zeroline   = TRUE,
            zerolinecolor = "#CCCCCC",
            showgrid   = TRUE,
            gridcolor  = "#EEEEEE"
          ),
          legend  = list(title = list(text = "<b>Population</b>")),
          paper_bgcolor = "white",
          plot_bgcolor  = "white",
          font          = list(family = "Helvetica Neue, Helvetica, Arial, sans-serif",
                               size   = 12)
        ) |>
        plotly::config(
          displayModeBar  = TRUE,
          modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
          displaylogo     = FALSE,
          toImageButtonOptions = list(
            format   = "png",
            scale    = 2,
            filename = "popgen_pca"
          )
        )
    })

    # ── CSV download ─────────────────────────────────────────────────────────
    output$dl_scores <- downloadHandler(
      filename = function() paste0("popgen_pca_scores_", Sys.Date(), ".csv"),
      content  = function(file) {
        res <- pca_result()
        req(res)
        write.csv(res$scores, file, row.names = FALSE)
      }
    )
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Individual summary card for the PCA panel
pca_card <- function(label, value, colour, icon_name) {
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
