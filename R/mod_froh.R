# =============================================================================
# mod_froh.R — Runs of Homozygosity (F-ROH) module (Phase 3, Goal 1)
#
# Computes per-individual F-ROH (proportion of genome in runs of homozygosity)
# and displays:
#   1. Summary cards: mean/max F-ROH per population, % individuals with ROH
#   2. Plotly boxplot: per-population F-ROH distribution with conservation
#      reference lines (equivalent inbreeding coefficients)
#   3. Per-individual data table with export
#
# Design decisions:
#  - Sliding-window ROH detection on SNP sequence (no physical map required)
#  - Missing sites do NOT break a run (follows PLINK convention)
#  - min_snps slider is module-local — does not affect other tabs
#  - F-ROH computation is cached per (upload × min_snps) combination
# =============================================================================

#' F-ROH module UI
#'
#' @param id Module namespace ID
frohUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("load_placeholder")),
    uiOutput(ns("scenario_banner")),
    uiOutput(ns("summary_cards")),

    # ── Controls ─────────────────────────────────────────────────────────────
    shinyjs::hidden(
      div(
        id = ns("froh_ready_content"),

        tags$hr(style = "margin: 12px 0;"),

        div(
          class = "d-flex align-items-center gap-4 flex-wrap mb-2",
          div(
            style = "flex: 0 0 auto;",
            tags$label(
              HTML("Min. ROH length (SNPs):"),
              class = "me-1 fw-semibold",
              style = "font-size:0.9em;"
            ),
            sliderInput(
              ns("min_snps"),
              label  = NULL,
              min    = 10, max = 200, value = 50, step = 10,
              ticks  = FALSE,
              width  = "220px"
            ),
            tags$p(
              style = "font-size:0.74em; color:#888; margin:-4px 0 0;",
              "Higher = only long, recent-inbreeding ROH counted."
            )
          ),
          div(
            style = "flex: 0 0 auto; margin-top:22px;",
            downloadButton(
              ns("dl_froh"),
              label = tagList(icon("download"), " Export CSV"),
              class = "btn btn-outline-secondary btn-sm"
            )
          )
        ),

        # ── Boxplot ─────────────────────────────────────────────────────────
        plotly::plotlyOutput(ns("froh_plot"), height = "420px"),

        # ── Per-individual table ─────────────────────────────────────────────
        tags$h6(
          "Per-individual F-ROH",
          style = "font-weight:600; color:#1B3A5C; margin: 14px 0 6px;"
        ),
        DT::dataTableOutput(ns("froh_table"))
      )
    ),

    # ── Footnotes ─────────────────────────────────────────────────────────────
    tags$p(
      class = "mt-2",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "F<sub>ROH</sub> = (SNPs in ROH) / (total informative SNPs). ",
        "A run is defined as a consecutive sequence of homozygous SNPs ",
        "of at least the specified minimum length. ",
        "Missing genotype calls bridge a run without contributing to its SNP count — ",
        "only called homozygous SNPs are counted toward run length and toward the denominator. ",
        "Values are comparable within a dataset but depend on SNP density."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML(paste0(
        "&#9432; Reference pedigree-equivalent inbreeding coefficients (diploid): ",
        "F<sub>ROH</sub> ≥ 0.25 ≈ full-sibling or parent-offspring mating (F = 1/4); ",
        "≥ 0.125 ≈ half-sibling mating (F = 1/8); ",
        "≥ 0.0625 ≈ first-cousin mating (F = 1/16). ",
        "Background F<sub>ROH</sub> in outbred populations is typically 0.01–0.05. ",
        "Unlike F = (He−Ho)/He, F<sub>ROH</sub> specifically captures ",
        "<em>recent</em> inbreeding (within the last ~10 generations) and is ",
        "not confounded by ancient population structure or allele-frequency skew. ",
        "Without a physical map, long ROH cannot be distinguished by genomic position — ",
        "interpret as a relative, not absolute, measure."
      ))
    )
  )
}

#' F-ROH module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive returning the upload module's data list
#'                    ($ready, $gt, $metadata, $is_demo, $scenario_info)
frohServer <- function(id, upload_data) {
  moduleServer(id, function(input, output, session) {

    # ── Load placeholder ─────────────────────────────────────────────────────
    output$load_placeholder <- renderUI({
      d <- upload_data()
      if (isTRUE(d$ready)) return(NULL)
      div(
        class = "alert alert-secondary",
        style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
        icon("info-circle"), " ",
        "Load a VCF file and metadata CSV from the sidebar to compute F-ROH."
      )
    })

    # ── Show/hide controls ───────────────────────────────────────────────────
    observeEvent(upload_data(), {
      shinyjs::toggle("froh_ready_content", condition = isTRUE(upload_data()$ready))
    })

    # ── Compute F-ROH (reactive — reruns when data OR min_snps changes) ──────
    froh_result <- reactive({
      d <- upload_data()
      req(d$ready)
      min_s <- input$min_snps %||% 50L
      tryCatch(
        calc_froh(d$gt, d$metadata, min_snps = as.integer(min_s)),
        error = function(e) {
          showNotification(
            paste0("F-ROH error: ", e$message),
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
      df <- froh_result()
      req(df)
      pop_sum <- summarise_froh(df)

      valid_froh   <- df$froh[!is.na(df$froh)]
      overall_mean <- if (length(valid_froh) > 0) round(mean(valid_froh), 3) else NA_real_
      max_idx      <- which.max(df$froh)   # returns integer(0) only when all values are NA
      max_indiv    <- if (length(max_idx) > 0) df$sample_id[max_idx] else "—"
      max_val      <- if (length(valid_froh) > 0) round(max(valid_froh), 3) else NA_real_
      pct_with_roh <- if (nrow(df) > 0) round(100 * mean(df$n_roh > 0, na.rm = TRUE), 1) else 0
      n_high_risk  <- sum(pop_sum$mean_froh >= 0.125, na.rm = TRUE)
      hr_col       <- if (n_high_risk > 0) "#D32F2F" else "#388E3C"

      # Guard against NA values: paste0(NA_real_) produces the literal string "NA"
      mean_str <- if (is.na(overall_mean)) "—" else as.character(overall_mean)
      max_str  <- if (is.na(max_val))      "—" else as.character(max_val)

      div(
        class = "row g-2 mb-2",
        froh_card("Mean F-ROH",
                  mean_str,                      "#1B3A5C", "dna"),
        froh_card("Max F-ROH (individual)",
                  max_str,                       "#D32F2F", "exclamation-triangle"),
        froh_card("Individuals with ROH",
                  paste0(pct_with_roh, "%"),     "#F57C00", "user"),
        froh_card("High-risk populations",
                  n_high_risk,                   hr_col,   "shield-alt")
      )
    })

    # ── Boxplot ──────────────────────────────────────────────────────────────
    output$froh_plot <- plotly::renderPlotly({
      df <- froh_result()
      req(df)

      pops    <- unique(df$population)
      n_pops  <- length(pops)
      palette <- stats::setNames(
        PCA_PALETTE[((seq_along(pops) - 1L) %% length(PCA_PALETTE)) + 1L],
        pops
      )

      fig <- plotly::plot_ly()

      for (pop in pops) {
        sub <- df[df$population == pop, ]
        fig <- plotly::add_trace(
          fig,
          type  = "box",
          y     = sub$froh,
          name  = pop,
          boxpoints = "all",
          jitter    = 0.35,
          pointpos  = 0,
          marker    = list(size = 5, opacity = 0.7),
          line      = list(color = palette[pop]),
          fillcolor = paste0(palette[pop], "40"),   # 25% opacity fill
          text      = paste0(sub$sample_id,
                             "<br/>F-ROH: ", round(sub$froh, 4),
                             "<br/>ROH runs: ", sub$n_roh),
          hoverinfo = "text"
        )
      }

      # Reference lines and annotations for inbreeding equivalents
      ref_lines <- list(
        froh_hline(0.0625, "#FBC02D"),
        froh_hline(0.125,  "#F57C00"),
        froh_hline(0.25,   "#D32F2F")
      )
      ref_annotations <- list(
        froh_ann(0.0625, "F-ROH = 0.0625 (1st-cousin)",         "#8B6914"),
        froh_ann(0.125,  "F-ROH = 0.125 (half-sibling)",        "#F57C00"),
        froh_ann(0.25,   "F-ROH = 0.25 (full-sib / parent-off)","#D32F2F")
      )

      fig |>
        plotly::layout(
          yaxis  = list(
            title    = "F-ROH",
            range    = c(-0.02, {
              y_hi <- max(df$froh, na.rm = TRUE) * 1.15 + 0.05
              if (!is.finite(y_hi)) 0.35 else min(1.0, y_hi)   # guard: all-NA → -Inf
            }),
            zeroline = TRUE,
            zerolinecolor = "#CCCCCC",
            showgrid = TRUE,
            gridcolor = "#EEEEEE"
          ),
          xaxis       = list(title = "", showgrid = FALSE),
          shapes      = ref_lines,
          annotations = ref_annotations,
          showlegend  = FALSE,
          paper_bgcolor = "white",
          plot_bgcolor  = "white",
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
            filename = "popgen_froh"
          )
        )
    })

    # ── Per-individual table ─────────────────────────────────────────────────
    output$froh_table <- DT::renderDataTable({
      df <- froh_result()
      req(df)
      display <- data.frame(
        "Sample ID"    = df$sample_id,
        "Population"   = df$population,
        "F-ROH"        = round(df$froh, 4),
        "ROH runs"     = df$n_roh,
        "ROH SNPs"     = df$total_roh_snps,
        "Informative SNPs" = df$n_informative_snps,
        "Risk"         = froh_risk_label(df$froh),
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
          columnDefs = list(list(className = "dt-center", targets = 2:5))
        )
      ) |>
        DT::formatStyle(
          "Risk",
          backgroundColor = DT::styleEqual(
            c("Severe — equivalent to full-sib or parent-offspring",
              "High — equivalent to half-sibling mating",
              "Moderate — equivalent to first-cousin mating",
              "Low — background homozygosity",
              "None detected", "Unknown"),
            c("#FFCDD2", "#FFE0B2", "#FFF9C4", "#E8F5E9", "#F5F5F5", "#F5F5F5")
          )
        )
    })

    # ── CSV download ─────────────────────────────────────────────────────────
    output$dl_froh <- downloadHandler(
      filename = function() paste0("popgen_froh_", Sys.Date(), ".csv"),
      content  = function(file) {
        df <- froh_result()
        req(df)
        df$risk <- froh_risk_label(df$froh)
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Individual summary card for the F-ROH panel
froh_card <- function(label, value, colour, icon_name) {
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

#' Plotly horizontal reference line for inbreeding equivalents
#'
#' @param y       Y-axis position (F-ROH value)
#' @param colour  Hex line colour
froh_hline <- function(y, colour) {
  list(
    type  = "line",
    x0    = 0, x1    = 1,
    y0    = y, y1    = y,
    xref  = "paper", yref = "y",
    line  = list(color = colour, width = 1.2, dash = "dot")
  )
}

#' Plotly annotation label for an F-ROH reference line
#'
#' @param y      Y-axis position (F-ROH value)
#' @param label  Short text label
#' @param colour Hex text colour
froh_ann <- function(y, label, colour) {
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
