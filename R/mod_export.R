# =============================================================================
# mod_export.R — Export panel module (Phase 5)
#
# Provides four export pathways:
#   1. HTML Conservation Report — full narrative rendered from report_template.Rmd
#   2. PDF Conservation Report  — HTML report converted via pagedown::chrome_print()
#                                  (falls back to HTML + user message if unavailable)
#   3. XLSX data workbook       — stats, FST matrix, GEA results, offset table
#   4. Statistics CSV           — per-population diversity table
#
# Plotly charts carry a built-in camera button (PNG at scale=2, set via
# toImageButtonOptions in each chart module) — no separate server handler needed.
# The Leaflet map screenshot is handled client-side by leaflet.extras2.
# =============================================================================

#' Export panel UI
#'
#' @param id Module namespace ID
exportUI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$h5(
      tagList(icon("file-export"), " Export & Reports"),
      style = "color:#1B3A5C; font-weight:700; margin-bottom:16px;"
    ),
    tags$p(
      style = "font-size:0.88em; color:#555; margin-bottom:20px;",
      "Download your results as a report or raw data files. All exports reflect ",
      "the data currently loaded in the app."
    ),

    # ── Conservation Report (HTML + PDF) ─────────────────────────────────────
    div(
      class = "export-card",
      div(
        class = "export-card-icon",
        style = "background:#E3F2FD;",
        icon("file-medical-alt", style = "color:#1565C0; font-size:1.5em;")
      ),
      div(
        class = "export-card-body",
        tags$h6("Conservation Report", style = "margin:0 0 4px; font-weight:700;"),
        tags$p(
          class = "export-card-desc",
          "Full report: risk summary, diversity table, F\\u209B\\u209C matrix, inbreeding, ",
          "Ne, GEA and genomic offset — ready to attach to a conservation management plan."
        ),
        div(
          class = "d-flex gap-2 flex-wrap",
          shinyjs::disabled(
            downloadButton(
              ns("dl_report_html"),
              label = tagList(icon("file-code"), " HTML Report"),
              class = "btn btn-primary btn-sm"
            )
          ),
          shinyjs::disabled(
            downloadButton(
              ns("dl_report_pdf"),
              label = tagList(icon("file-pdf"), " PDF Report"),
              class = "btn btn-outline-primary btn-sm"
            )
          )
        ),
        tags$p(
          style = "font-size:0.75em; color:#999; margin-top:6px;",
          icon("info-circle"), " ",
          "PDF requires the ", tags$code("pagedown"), " package and Google Chrome/Chromium. ",
          "Alternatively, open the HTML report in your browser and use File → Print → Save as PDF."
        )
      )
    ),

    tags$hr(style = "margin:18px 0;"),

    # ── XLSX Data Workbook ───────────────────────────────────────────────────
    div(
      class = "export-card",
      div(
        class = "export-card-icon",
        style = "background:#E8F5E9;",
        icon("file-excel", style = "color:#217346; font-size:1.5em;")
      ),
      div(
        class = "export-card-body",
        tags$h6("Data Workbook (Excel)", style = "margin:0 0 4px; font-weight:700;"),
        tags$p(
          class = "export-card-desc",
          "Multi-sheet Excel workbook: diversity statistics, pairwise F\\u209B\\u209C matrix, ",
          "LFMM2 results, and genomic offset table — ready for further analysis."
        ),
        div(
          class = "d-flex gap-2 flex-wrap",
          shinyjs::disabled(
            downloadButton(
              ns("dl_xlsx"),
              label = tagList(icon("file-excel"), " Download XLSX"),
              class = "btn btn-success btn-sm"
            )
          )
        ),
        tags$p(
          style = "font-size:0.75em; color:#999; margin-top:6px;",
          icon("info-circle"), " ",
          "Requires the ", tags$code("openxlsx"), " package: ",
          tags$code("install.packages('openxlsx')")
        )
      )
    ),

    tags$hr(style = "margin:18px 0;"),

    # ── Statistics CSV ───────────────────────────────────────────────────────
    div(
      class = "export-card",
      div(
        class = "export-card-icon",
        style = "background:#FFF3E0;",
        icon("table", style = "color:#E65100; font-size:1.5em;")
      ),
      div(
        class = "export-card-body",
        tags$h6("Statistics Table (CSV)", style = "margin:0 0 4px; font-weight:700;"),
        tags$p(
          class = "export-card-desc",
          "Per-population diversity metrics (He, Ho, \\u03c0, F, Ne) as a flat CSV — ",
          "import directly into Excel, R, or Python."
        ),
        shinyjs::disabled(
          downloadButton(
            ns("dl_stats_csv"),
            label = tagList(icon("file-csv"), " Download CSV"),
            class = "btn btn-warning btn-sm"
          )
        )
      )
    ),

    tags$hr(style = "margin:18px 0;"),

    # ── Map PNG ──────────────────────────────────────────────────────────────
    div(
      class = "export-card",
      div(
        class = "export-card-icon",
        style = "background:#F3E5F5;",
        icon("map", style = "color:#6A1B9A; font-size:1.5em;")
      ),
      div(
        class = "export-card-body",
        tags$h6("Map Screenshot (PNG)", style = "margin:0 0 4px; font-weight:700;"),
        tags$p(
          class = "export-card-desc",
          "Click the ",
          tags$strong("\U0001f4f7 camera icon"),
          " in the top-left corner of the Map to save the current view as PNG — ",
          "captures markers, FST network, ADMIXTURE pies, and the legend."
        ),
        actionButton(
          ns("go_to_map"),
          label = tagList(icon("map-marked-alt"), " Go to Map"),
          class = "btn btn-outline-secondary btn-sm"
        )
      )
    ),

    tags$hr(style = "margin:18px 0;"),

    # ── Chart PNG / SVG ──────────────────────────────────────────────────────
    div(
      class = "export-card",
      div(
        class = "export-card-icon",
        style = "background:#FCE4EC;",
        icon("chart-bar", style = "color:#880E4F; font-size:1.5em;")
      ),
      div(
        class = "export-card-body",
        tags$h6("Chart Images (PNG)", style = "margin:0 0 4px; font-weight:700;"),
        tags$p(
          class = "export-card-desc",
          "Every Plotly chart (PCA, ADMIXTURE bar, Ne, F-ROH, ",
          "GEA Manhattan, Offset bar) has a built-in download toolbar."
        ),
        tags$ol(
          style = "font-size:0.83em; color:#555; margin:6px 0 0 16px;",
          tags$li("Hover over the chart — a toolbar appears in the top-right corner"),
          tags$li(HTML("Click the <strong>camera icon</strong> to save as PNG at 2× (print-quality) resolution")),
          tags$li(HTML("The FST heatmap is an HTML table — use your browser's <em>Print → Save as PDF</em> or a screenshot tool to capture it"))
        )
      )
    ),

    tags$br(),

    # ── Status message ───────────────────────────────────────────────────────
    uiOutput(ns("export_status"))
  )
}

#' Export panel server
#'
#' @param id           Module namespace ID
#' @param pop_stats    Reactive: per-population diversity data frame
#' @param fst_mat      Reactive: pairwise FST matrix (or NULL)
#' @param pop_q        Reactive: per-population mean Q matrix (or NULL)
#' @param ne_stats     Reactive: Ne estimates data frame (or NULL)
#' @param gea_data     Reactive: GEA results list ($rda_result, $lfmm_result)
#' @param offset_data  Reactive: genomic offset data frame (or NULL)
#' @param upload_data  Reactive: upload module list (metadata, is_demo, etc.)
#' @param parent_session  Main session, needed so "Go to Map" can switch top-level tab
exportServer <- function(id, pop_stats, fst_mat = reactive(NULL),
                         pop_q = reactive(NULL), ne_stats = reactive(NULL),
                         gea_data = reactive(NULL), offset_data = reactive(NULL),
                         upload_data = reactive(list()),
                         parent_session = NULL) {
  moduleServer(id, function(input, output, session) {

    # ── Activate download buttons when data is ready ─────────────────────────
    observe({
      ready <- !is.null(pop_stats()) && nrow(pop_stats()) > 0
      if (ready) {
        shinyjs::enable("dl_report_html")
        shinyjs::enable("dl_report_pdf")
        shinyjs::enable("dl_stats_csv")
        shinyjs::enable("dl_xlsx")
      } else {
        shinyjs::disable("dl_report_html")
        shinyjs::disable("dl_report_pdf")
        shinyjs::disable("dl_stats_csv")
        shinyjs::disable("dl_xlsx")
      }
    })

    # ── Status banner ─────────────────────────────────────────────────────────
    output$export_status <- renderUI({
      df <- pop_stats()
      if (is.null(df) || nrow(df) == 0) {
        div(
          class = "alert alert-warning",
          style = "font-size:0.83em; padding:8px 12px;",
          icon("exclamation-triangle"), " ",
          "Load genomic data first — export buttons will activate once data is ready."
        )
      } else {
        phases_done <- character(0)
        if (!is.null(fst_mat()))               phases_done <- c(phases_done, "FST")
        if (!is.null(ne_stats()))              phases_done <- c(phases_done, "Ne")
        if (!is.null(gea_data()$rda_result))   phases_done <- c(phases_done, "RDA")
        if (!is.null(gea_data()$lfmm_result))  phases_done <- c(phases_done, "LFMM2")
        if (!is.null(offset_data()))           phases_done <- c(phases_done, "Offset")

        div(
          class = "alert alert-success",
          style = "font-size:0.83em; padding:8px 12px;",
          icon("check-circle"), " ",
          sprintf("%d population%s loaded.",
                  nrow(df), if (nrow(df) == 1) "" else "s"),
          if (length(phases_done) > 0)
            paste0(" Results available: ", paste(phases_done, collapse = ", "), ".")
          else
            " Run analyses in their tabs to include them in the report."
        )
      }
    })

    # ── HTML Conservation Report ──────────────────────────────────────────────
    output$dl_report_html <- downloadHandler(
      filename = function() {
        paste0("conservation_report_", format(Sys.Date(), "%Y%m%d"), ".html")
      },
      content = function(file) {
        df <- pop_stats()
        validate(need(!is.null(df) && nrow(df) > 0,
                      "No data loaded — please upload a VCF and coordinates CSV first."))

        report_params <- build_report_params(df, fst_mat, pop_q, ne_stats,
                                             gea_data, offset_data, upload_data)

        tmp_rmd <- file.path(tempdir(), "report_template.Rmd")
        file.copy(file.path("R", "report_template.Rmd"), tmp_rmd, overwrite = TRUE)

        tryCatch({
          rmarkdown::render(
            input       = tmp_rmd,
            output_file = file,
            params      = report_params,
            envir       = new.env(parent = globalenv()),
            quiet       = TRUE
          )
        }, error = function(e) {
          writeLines(minimal_fallback_html(df, e$message), con = file)
        })
      }
    )

    # ── PDF Conservation Report ───────────────────────────────────────────────
    output$dl_report_pdf <- downloadHandler(
      filename = function() {
        # If pagedown is unavailable the content function falls back to HTML;
        # change the extension here so the downloaded file has the right type.
        ext <- if (requireNamespace("pagedown", quietly = TRUE)) "pdf" else "html"
        paste0("conservation_report_", format(Sys.Date(), "%Y%m%d"), ".", ext)
      },
      content = function(file) {
        df <- pop_stats()
        validate(need(!is.null(df) && nrow(df) > 0,
                      "No data loaded."))

        # Check for pagedown availability (needs Chrome/Chromium too)
        has_pagedown <- requireNamespace("pagedown", quietly = TRUE)

        if (!has_pagedown) {
          # Fallback: render HTML to file — filename() already returns .html in this case
          report_params <- build_report_params(df, fst_mat, pop_q, ne_stats,
                                               gea_data, offset_data, upload_data)
          tmp_rmd <- file.path(tempdir(), "report_template.Rmd")
          file.copy(file.path("R", "report_template.Rmd"), tmp_rmd, overwrite = TRUE)
          tryCatch({
            rmarkdown::render(tmp_rmd, output_file = file,
                              params = report_params,
                              envir  = new.env(parent = globalenv()),
                              quiet  = TRUE)
          }, error = function(e) {
            writeLines(minimal_fallback_html(df,
              paste0("pagedown not installed; rmarkdown render also failed: ", e$message)),
              con = file)
          })
          showNotification(
            paste0("pagedown package not found. Saved as HTML instead. ",
                   "Install pagedown for true PDF: install.packages('pagedown')"),
            type = "warning", duration = 12
          )
          return(invisible(NULL))
        }

        # Render HTML first, then convert to PDF
        report_params <- build_report_params(df, fst_mat, pop_q, ne_stats,
                                             gea_data, offset_data, upload_data)
        tmp_rmd  <- file.path(tempdir(), "report_template.Rmd")
        tmp_html <- file.path(tempdir(),
                               paste0("conservation_report_", Sys.getpid(), ".html"))
        file.copy(file.path("R", "report_template.Rmd"), tmp_rmd, overwrite = TRUE)

        tryCatch({
          rmarkdown::render(
            input       = tmp_rmd,
            output_file = tmp_html,
            params      = report_params,
            envir       = new.env(parent = globalenv()),
            quiet       = TRUE
          )
          pagedown::chrome_print(
            input  = tmp_html,
            output = file,
            wait   = 5,
            timeout = 120,
            verbose = FALSE
          )
        }, error = function(e) {
          showNotification(
            paste0("PDF conversion failed: ", e$message,
                   ". Falling back to HTML. Open in browser and use File → Print → Save as PDF."),
            type = "error", duration = 15
          )
          # Copy the HTML as the fallback file
          if (file.exists(tmp_html)) {
            file.copy(tmp_html, file, overwrite = TRUE)
          } else {
            writeLines(minimal_fallback_html(df, e$message), con = file)
          }
        })
      }
    )

    # ── XLSX Data Workbook ────────────────────────────────────────────────────
    output$dl_xlsx <- downloadHandler(
      filename = function() {
        paste0("popgen_data_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
      },
      content = function(file) {
        df <- pop_stats()
        validate(need(!is.null(df) && nrow(df) > 0, "No data loaded."))

        if (!requireNamespace("openxlsx", quietly = TRUE)) {
          showNotification(
            paste0("openxlsx package required for Excel export. ",
                   "Install with: install.packages('openxlsx')"),
            type = "error", duration = 10
          )
          req(FALSE)   # cancels the download cleanly; avoids writing a 0-byte file
        }

        wb <- openxlsx::createWorkbook()

        # ── Sheet 1: Diversity statistics ────────────────────────────────────
        stats_out <- df
        stats_out$risk_label <- he_risk_label(df$He)
        ne_df <- ne_stats()
        if (!is.null(ne_df) && "ne" %in% names(ne_df)) {
          stats_out$Ne      <- ne_df$ne[match(stats_out$population, ne_df$population)]
          stats_out$Ne_risk <- ne_risk_label(stats_out$Ne)
        }
        openxlsx::addWorksheet(wb, "Diversity Statistics")
        openxlsx::writeData(wb, "Diversity Statistics", stats_out)
        openxlsx::setColWidths(wb, "Diversity Statistics", cols = seq_len(ncol(stats_out)),
                               widths = "auto")

        # ── Sheet 2: Pairwise FST ────────────────────────────────────────────
        fst <- fst_mat()
        if (!is.null(fst) && nrow(fst) >= 2) {
          fst_tbl <- as.data.frame(round(fst, 4))
          fst_tbl <- cbind(Population = rownames(fst_tbl), fst_tbl)
          openxlsx::addWorksheet(wb, "Pairwise FST")
          openxlsx::writeData(wb, "Pairwise FST", fst_tbl)
          openxlsx::setColWidths(wb, "Pairwise FST",
                                  cols = seq_len(ncol(fst_tbl)), widths = "auto")
        }

        # ── Sheet 3: LFMM2 results (all significant loci) ────────────────────
        gea <- gea_data()
        if (!is.null(gea$lfmm_result)) {
          lfmm_sig <- gea$lfmm_result[gea$lfmm_result$significant, , drop = FALSE]
          if (nrow(lfmm_sig) == 0) lfmm_sig <- gea$lfmm_result  # all loci if none sig
          openxlsx::addWorksheet(wb, "LFMM2 Results")
          openxlsx::writeData(wb, "LFMM2 Results", lfmm_sig)
          openxlsx::setColWidths(wb, "LFMM2 Results",
                                  cols = seq_len(ncol(lfmm_sig)), widths = "auto")
        }

        # ── Sheet 4: Genomic offset ──────────────────────────────────────────
        off_df <- offset_data()
        if (!is.null(off_df) && nrow(off_df) > 0) {
          # Drop colour column — not useful in a spreadsheet
          off_export <- off_df[, setdiff(names(off_df), "risk_colour"), drop = FALSE]
          off_export  <- off_export[order(off_export$offset_norm, decreasing = TRUE), ]
          openxlsx::addWorksheet(wb, "Genomic Offset")
          openxlsx::writeData(wb, "Genomic Offset", off_export)
          openxlsx::setColWidths(wb, "Genomic Offset",
                                  cols = seq_len(ncol(off_export)), widths = "auto")
        }

        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )

    # ── Statistics CSV ────────────────────────────────────────────────────────
    output$dl_stats_csv <- downloadHandler(
      filename = function() {
        paste0("popgen_stats_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        df <- pop_stats()
        validate(need(!is.null(df) && nrow(df) > 0, "No data loaded."))

        out            <- df
        out$risk_label <- he_risk_label(df$He)

        ne_df <- ne_stats()
        if (!is.null(ne_df) && "ne" %in% names(ne_df)) {
          out$Ne      <- ne_df$ne[match(out$population, ne_df$population)]
          out$Ne_risk <- ne_risk_label(out$Ne)
        }
        write.csv(out, file, row.names = FALSE)
      }
    )

    # ── "Go to Map" button ────────────────────────────────────────────────────
    observeEvent(input$go_to_map, {
      sess <- if (!is.null(parent_session)) parent_session else session
      updateTabsetPanel(session = sess, inputId = "main_tabs", selected = "tab_map")
    })
  })
}

# =============================================================================
# Internal helpers
# =============================================================================

#' Collect all reactive results into the flat params list the Rmd template expects
build_report_params <- function(df, fst_mat, pop_q, ne_stats,
                                gea_data, offset_data, upload_data) {
  list(
    pop_stats   = df,
    fst_mat     = fst_mat(),
    pop_q       = pop_q(),
    ne_stats    = ne_stats(),
    gea_data    = gea_data(),
    offset_data = offset_data(),
    is_demo     = isTRUE(upload_data()$is_demo),
    report_date = format(Sys.Date(), "%d %B %Y"),
    app_version = APP_VERSION
  )
}

#' Minimal plain HTML fallback — used when rmarkdown is unavailable or fails
#'
#' @param pop_stats_df Per-population statistics data frame
#' @param err_msg      Optional error message string to embed
minimal_fallback_html <- function(pop_stats_df, err_msg = NULL) {
  # Local helper — rounds numeric to string, NA → "—"
  fmt <- function(x, digits = 3) {
    ifelse(is.na(x), "&mdash;",
           formatC(round(as.numeric(x), digits), format = "f", digits = digits))
  }

  rows <- vapply(seq_len(nrow(pop_stats_df)), function(i) {
    r <- pop_stats_df[i, ]
    sprintf(
      "<tr><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
      htmltools::htmlEscape(as.character(r$population)),
      as.integer(r$n_samples %||% r$n),   # n_samples is the canonical name; r$n as fallback
      fmt(r$He), fmt(r$Ho), fmt(r$pi), fmt(r$inbreeding_f, 3),
      he_risk_label(r$He)
    )
  }, character(1))

  err_note <- if (!is.null(err_msg)) {
    sprintf(
      "<p style='color:#B00;font-size:0.85em;'>Note: rmarkdown render failed — %s. ",
      htmltools::htmlEscape(err_msg)
    )
  } else ""

  paste0(
    "<!DOCTYPE html><html><head><meta charset='utf-8'>",
    "<title>Conservation Genomics Report</title>",
    "<style>",
    "body{font-family:'Segoe UI',Arial,sans-serif;max-width:900px;margin:40px auto;padding:0 20px;color:#333;}",
    "h1{color:#1B3A5C;border-bottom:3px solid #2E75B6;padding-bottom:8px;}",
    "h2{color:#2E75B6;margin-top:1.8em;}",
    "table{border-collapse:collapse;width:100%;font-size:0.9em;margin-top:12px;}",
    "thead th{background:#1B3A5C;color:white;padding:9px 12px;text-align:left;}",
    "tbody td{padding:7px 12px;border-bottom:1px solid #eee;}",
    "tbody tr:nth-child(even){background:#fafafa;}",
    ".warn{background:#FFF3CD;border:1px solid #ffc107;border-radius:6px;padding:10px 14px;",
    "       margin-bottom:16px;font-size:0.88em;}",
    "</style></head><body>",
    "<h1>Conservation Genomics Report</h1>",
    "<p>Generated: ", format(Sys.Date(), "%d %B %Y"),
    " &nbsp;|&nbsp; <strong>PopGen Map v", APP_VERSION, "</strong></p>",
    if (!is.null(err_msg)) paste0(
      '<div class="warn">⚠️ Full report rendering failed: ',
      htmltools::htmlEscape(err_msg),
      '<br/>To generate the full report, ensure <code>rmarkdown</code> and <code>knitr</code> ',
      'are installed: <code>install.packages(c("rmarkdown","knitr"))</code></div>'
    ) else "",
    "<h2>Population Diversity Statistics</h2>",
    "<table><thead><tr>",
    "<th>Population</th><th>n</th><th>He</th><th>Ho</th>",
    "<th>&pi;</th><th>F</th><th>Risk (He)</th>",
    "</tr></thead><tbody>",
    paste(rows, collapse = "\n"),
    "</tbody></table>",
    "<p style='font-size:0.8em;color:#999;margin-top:40px;'>",
    "PopGen Map v", APP_VERSION, " &nbsp;|&nbsp; MIT Licence &nbsp;|&nbsp;",
    "<a href='https://github.com/csunbird/population-genomics-map'>GitHub</a></p>",
    "</body></html>"
  )
}
