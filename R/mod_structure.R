# =============================================================================
# mod_structure.R — Population structure module (Phase 2, Goal 1)
#
# Displays a pairwise FST heatmap computed with the Hudson et al. (1992)
# ratio-of-averages estimator. All computation is reactive and runs only
# when new data are loaded — metric display changes do not re-trigger it.
#
# Design decisions:
#  - FST is computed once per upload, then cached in a reactive
#  - Heatmap is a custom HTML table (no extra dependencies, full colour control)
#  - Colour palette: white (0) → teal (1) matching standard heatmap convention
#  - Conservation thresholds from Wright (1978) shown as indicative guidance
# =============================================================================

#' Structure module UI
#'
#' @param id Module namespace ID
structureUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("scenario_banner")),
    uiOutput(ns("summary_cards")),
    tags$hr(style = "margin: 12px 0;"),

    div(
      class = "d-flex justify-content-between align-items-center mb-2",
      tags$h6(
        "Pairwise Fₛₜ (Hudson et al. 1992)",
        style = "margin:0; font-weight:600; color:#1B3A5C;"
      ),
      downloadButton(
        ns("dl_fst"),
        label = tagList(icon("download"), " Export CSV"),
        class = "btn btn-outline-secondary btn-sm"
      )
    ),

    uiOutput(ns("fst_heatmap")),

    tags$p(
      class = "mt-2",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "Fₛₜ estimated using the Hudson et al. (1992) ratio-of-averages method ",
        "(sum of numerators / sum of denominators across loci). ",
        "Values range from 0 (no differentiation) to 1 (complete differentiation). ",
        "Diagonal = 0 by definition."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML(paste0(
        "&#9432; Wright (1978) guidelines: Fₛₜ &lt; 0.05 = little; ",
        "0.05&#x2013;0.15 = moderate; ",
        "0.15&#x2013;0.25 = great; ",
        "&gt; 0.25 = very great differentiation. ",
        "These thresholds were calibrated on allozyme data; ",
        "SNP-based Fₛₜ typically runs higher. Treat categories as indicative, not definitive."
      ))
    )
  )
}

#' Structure module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive returning the upload module's data list
#'                    (fields: $ready, $gt, $metadata, $is_demo, $scenario_info)
#' @param pop_stats   Reactive returning per-population statistics data frame
#' @param fst_mat     Reactive returning the pairwise FST matrix from server.R
#'                    (pre-computed at the app level to avoid redundant computation
#'                     between this module and the map FST network overlay)
structureServer <- function(id, upload_data, pop_stats, fst_mat) {
  moduleServer(id, function(input, output, session) {

    # fst_mat is now injected from server.R — no local computation needed.
    # Alias to fst_matrix for internal use (keeps the rest of the code unchanged).
    fst_matrix <- fst_mat

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
      mat <- fst_matrix()
      req(mat)

      n_pop <- nrow(mat)

      # Collect upper-triangle values (each pair counted once)
      ut_idx  <- which(upper.tri(mat), arr.ind = TRUE)
      ut_vals <- mat[ut_idx]
      valid   <- ut_vals[!is.na(ut_vals)]

      mean_fst <- if (length(valid) > 0) mean(valid) else NA_real_

      # Most-differentiated pair
      max_pair_label <- if (length(valid) > 0) {
        best     <- which.max(ut_vals)   # NA-safe: ignores NA positions
        if (!is.na(ut_vals[best])) {
          pop_r <- rownames(mat)[ut_idx[best, 1]]
          pop_c <- colnames(mat)[ut_idx[best, 2]]
          lbl   <- paste0(pop_r, " / ", pop_c)
          # Truncate very long labels for the card display
          if (nchar(lbl) > 28) paste0(substr(lbl, 1, 25), "…") else lbl
        } else "—"
      } else "—"

      max_fst_val <- if (length(valid) > 0) max(valid) else NA_real_

      # Pairs with very great differentiation (FST > 0.25)
      n_high  <- sum(valid > FST_THRESHOLDS$great, na.rm = TRUE)
      hi_col  <- if (n_high > 0) "#D32F2F" else "#388E3C"
      max_col <- if (!is.na(max_fst_val) && max_fst_val > FST_THRESHOLDS$great) "#D32F2F" else "#388E3C"

      div(
        class = "row g-2 mb-2",
        fst_card("Populations",       n_pop,               "#1B3A5C", "users"),
        fst_card("Mean Fₛₜ", fmt(mean_fst),       "#2E75B6", "chart-bar"),
        fst_card("Highest-diff. pair", max_pair_label,      max_col,   "arrows-alt-h"),
        fst_card("Very-high-diff. pairs", n_high,           hi_col,    "exclamation-triangle")
      )
    })

    # ── FST Heatmap ──────────────────────────────────────────────────────────
    output$fst_heatmap <- renderUI({
      mat <- fst_matrix()
      if (is.null(mat)) {
        return(div(
          class = "legend-box",
          style = "color:#888; font-size:0.85em;",
          "Load data to compute the pairwise Fₛₜ matrix."
        ))
      }
      render_fst_heatmap(mat)
    })

    # ── CSV download ─────────────────────────────────────────────────────────
    output$dl_fst <- downloadHandler(
      filename = function() paste0("popgen_fst_matrix_", Sys.Date(), ".csv"),
      content  = function(file) {
        mat <- fst_matrix()
        req(mat)
        df_out             <- as.data.frame(round(mat, 4))
        df_out             <- cbind(Population = rownames(df_out), df_out)
        write.csv(df_out, file, row.names = FALSE)
      }
    )
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Render the pairwise FST matrix as a colour-coded HTML table
#'
#' Pre-computes all cell colours in one vectorised call to fst_colour() for
#' efficiency, then builds the HTML table row by row.
#'
#' @param mat Named numeric matrix (populations × populations)
#' @return Shiny tagList containing the table and a colour legend
render_fst_heatmap <- function(mat) {
  pops  <- rownames(mat)
  n_pop <- length(pops)

  # Pre-compute all cell colours in one pass
  all_vals <- as.vector(mat)
  all_cols <- fst_colour(all_vals)
  col_mat  <- matrix(all_cols, nrow = n_pop, ncol = n_pop)

  # Truncate long population names in headers (title attr keeps full name)
  short_name <- function(p, max_chars = 14) {
    if (nchar(p) > max_chars) paste0(substr(p, 1, max_chars - 1), "…") else p
  }

  # Column header row — vertical text keeps the table compact
  header_cells <- lapply(pops, function(p) {
    tags$th(
      title = p,
      style = paste0(
        "padding:4px 6px; background:#1B3A5C; color:white; ",
        "font-size:0.75em; white-space:nowrap; ",
        "writing-mode:vertical-lr; transform:rotate(180deg); ",
        "min-width:36px; text-align:left; vertical-align:bottom;"
      ),
      short_name(p)
    )
  })
  header_row <- tags$tr(
    tags$th(style = "padding:4px 8px; background:#1B3A5C;", ""),
    do.call(tagList, header_cells)
  )

  # Data rows
  data_rows <- lapply(seq_len(n_pop), function(i) {
    cells <- lapply(seq_len(n_pop), function(j) {
      if (i == j) {
        # Diagonal — not informative, styled distinctly
        tags$td(
          style = paste0(
            "padding:5px 8px; background:#EEEEEE; text-align:center; ",
            "font-size:0.82em; color:#999;"
          ),
          "—"
        )
      } else {
        val     <- mat[i, j]
        bg      <- col_mat[i, j]
        # Dark text for pale cells, white text once background is sufficiently dark
        txt_col <- if (!is.na(val) && val >= 0.45) "white" else "#333"
        tags$td(
          title = if (!is.na(val)) paste0(fst_label(val), " (Fₛₜ = ", fmt(val, 3), ")") else "NA",
          style = paste0(
            "padding:5px 8px; background:", bg, "; ",
            "text-align:center; font-size:0.82em; color:", txt_col, ";"
          ),
          fmt(val, 3)
        )
      }
    })
    tags$tr(
      tags$th(
        title = pops[i],
        style = paste0(
          "padding:5px 10px; background:#f1f3f5; font-weight:600; ",
          "font-size:0.82em; white-space:nowrap; text-align:left;"
        ),
        short_name(pops[i])
      ),
      do.call(tagList, cells)
    )
  })

  # Threshold annotation below the gradient bar
  threshold_ticks <- div(
    class = "d-flex justify-content-between",
    style = "width:200px; font-size:0.72em; color:#888; margin-top:1px;",
    tags$span("0"),
    tags$span("0.05"),
    tags$span("0.15"),
    tags$span("0.25"),
    tags$span("1")
  )

  legend <- div(
    class = "d-flex align-items-start gap-3 mt-2 flex-wrap",
    style = "font-size:0.78em; color:#555;",
    div(
      tags$div(
        class = "d-flex align-items-center gap-2",
        tags$span("Low Fₛₜ (0)"),
        tags$div(
          style = paste0(
            "width:200px; height:10px; border-radius:4px; border:1px solid #ddd; ",
            "background: linear-gradient(to right, #FFFFFF, #E0F2F1, #26A69A, #004D40);"
          )
        ),
        tags$span("High (1)")
      ),
      threshold_ticks
    ),
    div(
      style = "color:#aaa; font-size:0.95em; font-style:italic;",
      HTML("little &lt; 0.05 | moderate 0.05&#x2013;0.15 | great 0.15&#x2013;0.25 | very great &gt; 0.25")
    )
  )

  tagList(
    div(
      style = "overflow-x: auto; margin-bottom: 4px;",
      tags$table(
        style = "border-collapse:collapse;",
        tags$thead(header_row),
        tags$tbody(do.call(tagList, data_rows))
      )
    ),
    legend
  )
}

#' Individual summary card for the structure panel
#'
#' Mirrors stat_card() in mod_stats_panel.R but defined here to keep modules
#' self-contained.
#'
#' @param label     Card label (shown below value)
#' @param value     Card value (shown large)
#' @param colour    Accent colour (left border + icon)
#' @param icon_name Font Awesome icon name
fst_card <- function(label, value, colour, icon_name) {
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
