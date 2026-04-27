# =============================================================================
# mod_stats_panel.R — Population statistics panel module
# Displays per-population diversity table and summary cards
# =============================================================================

#' Stats panel UI
#'
#' @param id Module namespace ID
statsPanelUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("scenario_banner")),
    uiOutput(ns("summary_cards")),
    tags$hr(style = "margin: 12px 0;"),
    div(
      class = "d-flex justify-content-between align-items-center mb-2",
      tags$h6("Population Statistics", style = "margin:0; font-weight:600; color:#1B3A5C;"),
      downloadButton(
        ns("dl_stats"),
        label = tagList(icon("download"), " Export CSV"),
        class = "btn btn-outline-secondary btn-sm"
      )
    ),
    DT::dataTableOutput(ns("stats_table")),
    tags$p(
      class = "mt-2",
      style = "font-size:0.78em; color:#777;",
      HTML("He = expected heterozygosity &nbsp;|&nbsp; Ho = observed heterozygosity &nbsp;|&nbsp;
           &pi; = nucleotide diversity (Nei 1978, diploid-corrected) &nbsp;|&nbsp; F = inbreeding coefficient")
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML("&#9432; Conservation thresholds (He) were calibrated on microsatellite data.
           SNP-based He at MAF &ge; 0.05 typically runs lower — treat risk categories as indicative, not definitive.
           SNPs count may differ between populations because each population independently excludes loci
           with &gt;50&thinsp;% missing data within that group.")
    )
  )
}

#' Stats panel server
#'
#' @param id        Module namespace ID
#' @param pop_stats Reactive returning a data frame of per-population statistics
#' @param is_demo   Reactive logical — is current data the built-in demo?
#' @param scenario_info Reactive character — description of current dataset
statsPanelServer <- function(id, pop_stats, is_demo, scenario_info) {
  moduleServer(id, function(input, output, session) {

    # ── Scenario banner ─────────────────────────────────────────────────
    output$scenario_banner <- renderUI({
      req(scenario_info())
      cls <- if (isTRUE(is_demo())) "alert alert-warning" else "alert alert-info"
      div(
        class = cls,
        style = "font-size:0.82em; padding:8px 12px; margin-bottom:10px;",
        icon("info-circle"), " ", scenario_info()
      )
    })

    # ── Summary cards ───────────────────────────────────────────────────
    output$summary_cards <- renderUI({
      df <- pop_stats()
      req(df, nrow(df) > 0)

      n_pops   <- nrow(df)
      n_inds   <- sum(df$n, na.rm = TRUE)
      mean_he  <- mean(df$He, na.rm = TRUE)
      # At-risk = He < 0.20 (moderate risk or worse per IUCN thresholds)
      at_risk  <- sum(df$He < HE_THRESHOLDS$vulnerable, na.rm = TRUE)

      div(
        class = "row g-2 mb-2",
        stat_card("Populations",    n_pops,             "#1B3A5C", "users"),
        stat_card("Individuals",    n_inds,             "#2E75B6", "person"),
        stat_card("Mean He",        fmt(mean_he),       "#388E3C", "dna"),
        stat_card("At-risk pops",   at_risk,
                  if (at_risk > 0) "#D32F2F" else "#388E3C", "exclamation-triangle")
      )
    })

    # ── Statistics table ────────────────────────────────────────────────
    output$stats_table <- DT::renderDataTable({
      df <- pop_stats()
      req(df, nrow(df) > 0)

      display <- data.frame(
        Population = df$population,
        n          = df$n,
        SNPs       = if ("n_loci" %in% names(df)) df$n_loci else NA_integer_,
        He         = fmt(df$He),
        Ho         = fmt(df$Ho),
        pi         = fmt(df$pi),
        F          = fmt(df$inbreeding_f, 3),
        Risk       = he_risk_label(df$He),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      # Rename pi column to \u03c0 after construction \u2014 backtick Unicode escapes
      # are not supported in R < 4.4, but string assignment works in all versions.
      names(display)[names(display) == "pi"] <- "\u03c0"

      # Enable search + pagination only when there are many populations
      n_rows <- nrow(display)
      tbl_dom <- if (n_rows > 10) "ftp" else "t"

      DT::datatable(
        display,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength = 10,
          dom        = tbl_dom,
          ordering   = TRUE,
          columnDefs = list(
            list(className = "dt-center", targets = 1:6)
          )
        )
      ) |>
        DT::formatStyle(
          "Risk",
          backgroundColor = DT::styleEqual(
            c("Critical — very low diversity", "High risk — low diversity",
              "Moderate risk", "Low risk — moderate diversity", "Healthy — high diversity",
              "Unknown"),
            c("#FFCDD2", "#FFE0B2", "#FFF9C4", "#C8E6C9", "#BBDEFB", "#F5F5F5")
          )
        )
    })

    # ── CSV download ─────────────────────────────────────────────────
    output$dl_stats <- downloadHandler(
      filename = function() paste0("popgen_map_stats_", Sys.Date(), ".csv"),
      content  = function(file) {
        df <- pop_stats()
        req(df)
        df$risk_label <- he_risk_label(df$He)
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}

# ── Helper: individual summary card ──────────────────────────────────────────
stat_card <- function(label, value, colour, icon_name) {
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
