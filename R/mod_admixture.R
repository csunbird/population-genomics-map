# =============================================================================
# mod_admixture.R — ADMIXTURE Q-matrix upload + visualisation (Phase 2 Goal 4)
#
# Accepts a Q-matrix CSV produced by ADMIXTURE (run outside the app) and
# displays:
#   1. A STRUCTURE-style stacked bar chart (Plotly) of individual ancestry
#   2. SVG pie chart markers overlaid on the Leaflet map (one per population,
#      proportions averaged within each population)
#
# Design decisions:
#  - Q-matrix is uploaded separately from VCF/metadata (module-local fileInput)
#  - Coordinates for pie markers come from pop_stats (already computed in
#    server.R; no duplication of coordinate logic)
#  - SVG pies are URL-encoded data URIs — zero extra R package dependencies
#  - The server function returns a reactive (pop_q) so mapServer can consume it
# =============================================================================

# Colour palette for up to 12 K ancestry components
ADMIX_K_PALETTE <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E"
)

# =============================================================================
# UI
# =============================================================================

#' ADMIXTURE module UI
#'
#' @param id Module namespace ID
admixtureUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("load_placeholder")),
    uiOutput(ns("scenario_banner")),

    # ── Upload panel ─────────────────────────────────────────────────────────
    div(
      class = "card border-0 shadow-sm mb-3",
      div(
        class = "card-body p-3",
        tags$h6("Upload Q-matrix",
                style = "font-weight:600; color:#1B3A5C; margin-bottom:8px;"),
        tags$p(
          style = "font-size:0.82em; color:#666; margin-bottom:8px;",
          HTML(paste0(
            "Upload the Q-matrix output from ADMIXTURE or STRUCTURE. ",
            "Expected CSV columns: <code>sample_id, population, K1, K2, ..., Kn</code>.<br/>",
            "Run ADMIXTURE externally; this module visualises the results. ",
            "The <em>population</em> column must match population names in your metadata."
          ))
        ),
        fileInput(
          ns("q_file"),
          label       = NULL,
          accept      = c(".csv", ".txt", ".Q"),
          width       = "100%",
          placeholder = "Select Q-matrix CSV…"
        ),
        uiOutput(ns("upload_status"))
      )
    ),

    # ── Summary cards ────────────────────────────────────────────────────────
    uiOutput(ns("summary_cards")),
    tags$hr(style = "margin: 12px 0;"),

    # ── Controls ─────────────────────────────────────────────────────────────
    div(
      class = "d-flex align-items-center gap-3 flex-wrap mb-2",
      div(
        style = "flex: 0 0 auto;",
        tags$label("Sort individuals by:",
                   class = "me-1 fw-semibold",
                   style = "font-size:0.9em;"),
        selectInput(
          ns("sort_by"),
          label    = NULL,
          choices  = c("Population" = "population", "K1 ancestry" = "k1"),
          selected = "population",
          width    = "160px"
        ),
        tags$p(
          style = "font-size:0.74em; color:#888; margin:2px 0 0;",
          "Ancestry sort is within each population group."
        )
      ),
      div(
        style = "flex: 0 0 auto; margin-top:22px;",
        downloadButton(
          ns("dl_qmat"),
          label = tagList(icon("download"), " Export CSV"),
          class = "btn btn-outline-secondary btn-sm"
        )
      )
    ),

    # ── STRUCTURE-style stacked bar chart ────────────────────────────────────
    plotly::plotlyOutput(ns("structure_bar"), height = "280px"),

    # ── K colour legend ──────────────────────────────────────────────────────
    uiOutput(ns("k_legend")),

    # ── Footnotes ────────────────────────────────────────────────────────────
    tags$p(
      class = "mt-2",
      style = "font-size:0.78em; color:#777;",
      HTML(paste0(
        "STRUCTURE-style bar chart: each vertical bar is one individual; ",
        "colour segments show ancestry proportions per K component. ",
        "Pie chart map overlays (toggle on Map tab) show <em>per-population mean</em> ancestry — ",
        "consult the bar chart to assess within-population variation, which the mean can mask ",
        "(e.g. a population with 50% K1 individuals and 50% K2 individuals looks identical ",
        "to a uniformly admixed population in the pie)."
      ))
    ),
    tags$p(
      class = "mt-1",
      style = "font-size:0.76em; color:#aaa; font-style:italic;",
      HTML(paste0(
        "&#9432; Choose K by comparing cross-validation (CV) error across runs — the optimal K ",
        "typically shows the lowest CV error. K should not exceed the number of sampled populations. ",
        "High admixture between populations may indicate recent or ongoing gene flow, or ",
        "a K value that is too high (over-parameterised model). ",
        "ADMIXTURE assumes Hardy-Weinberg equilibrium within ancestral populations; ",
        "inbred or structured source populations can bias ancestry estimates."
      ))
    )
  )
}

# =============================================================================
# Server
# =============================================================================

#' ADMIXTURE module server
#'
#' @param id          Module namespace ID
#' @param upload_data Reactive returning the upload module's data list
#'                    (fields: $ready, $gt, $metadata, $is_demo, $scenario_info)
#' @param pop_stats   Reactive returning per-population statistics (for coordinates)
#' @return A reactive returning the per-population mean Q data frame (or NULL)
#'         — consumed by mapServer for pie chart map overlay
admixtureServer <- function(id, upload_data, pop_stats) {
  moduleServer(id, function(input, output, session) {

    # ── Load placeholder (shown until VCF data are ready) ───────────────────
    output$load_placeholder <- renderUI({
      d <- upload_data()
      if (isTRUE(d$ready)) return(NULL)
      div(
        class = "alert alert-secondary",
        style = "font-size:0.85em; padding:10px 14px; margin-bottom:10px;",
        icon("info-circle"), " ",
        "Load a VCF file and metadata CSV from the sidebar first, then upload a Q-matrix below."
      )
    })

    # ── Demo Q-matrix: auto-load when demo data is active ───────────────────
    # This gives first-time users a working ADMIXTURE visualisation immediately.
    # The reactive merges the user-uploaded file with the demo fallback so that
    # q_data() is always defined when demo data is loaded.
    demo_q_active <- reactive({
      d <- upload_data()
      if (isTRUE(d$is_demo) && !is.null(d$demo_qmatrix) && is.null(input$q_file)) {
        d$demo_qmatrix
      } else {
        NULL
      }
    })

    # ── Parse uploaded Q-matrix (or use demo fallback) ──────────────────────
    q_data <- reactive({
      # Demo fallback: return the pre-built K=2 Q-matrix when in demo mode
      demo <- demo_q_active()
      if (!is.null(demo)) return(demo)
      # User-uploaded file
      f <- input$q_file
      req(f)
      result <- tryCatch(
        parse_qmatrix(f$datapath),
        error = function(e) {
          showNotification(
            paste0("Q-matrix error: ", e$message),
            type     = "error",
            duration = 12
          )
          NULL
        }
      )
      # Surface any row-sum warning that parse_qmatrix() attached as an attribute
      if (!is.null(result)) {
        warn_msg <- attr(result, "rowsum_warning")
        if (!is.null(warn_msg)) {
          showNotification(
            HTML(paste0(
              warn_msg,
              "Verify that all K columns are included and none are missing from the file."
            )),
            type     = "warning",
            duration = 15
          )
        }
      }
      result
    })

    # ── Upload status ────────────────────────────────────────────────────────
    output$upload_status <- renderUI({
      f    <- input$q_file
      demo <- demo_q_active()
      if (is.null(f) && is.null(demo)) {
        return(div(
          style = "font-size:0.8em; color:#888; margin-top:4px;",
          "No Q-matrix loaded — upload a CSV above."
        ))
      }
      qd <- q_data()
      if (is.null(qd)) return(NULL)
      k_cols  <- grep("^K[0-9]+$", names(qd), value = TRUE)
      q_pops  <- unique(qd$population)
      n_k     <- length(k_cols)

      # Check population name overlap with loaded metadata
      d       <- upload_data()
      overlap_msg <- NULL
      k_warn_msg  <- NULL
      if (isTRUE(d$ready)) {
        meta_pops   <- unique(as.character(d$metadata$population))
        matched     <- intersect(q_pops, meta_pops)
        unmatched   <- setdiff(q_pops, meta_pops)
        n_matched   <- length(matched)
        n_meta_pops <- length(meta_pops)

        if (n_matched == 0) {
          overlap_msg <- div(
            class = "alert alert-danger",
            style = "font-size:0.8em; padding:6px 10px; margin:6px 0 0;",
            icon("exclamation-triangle"), " ",
            HTML(paste0(
              "<strong>No population names match</strong> between Q-matrix and metadata. ",
              "Q-matrix populations: [", paste(head(q_pops, 4), collapse = ", "),
              if (length(q_pops) > 4) "…" else "", "]. ",
              "Metadata populations: [", paste(head(meta_pops, 4), collapse = ", "),
              if (length(meta_pops) > 4) "…" else "", "]. ",
              "Check for spelling differences. Pie chart overlay will not display."
            ))
          )
        } else if (length(unmatched) > 0) {
          overlap_msg <- div(
            class = "alert alert-warning",
            style = "font-size:0.8em; padding:6px 10px; margin:6px 0 0;",
            icon("exclamation-circle"), " ",
            HTML(paste0(
              n_matched, " of ", length(q_pops), " Q-matrix populations matched metadata. ",
              "Unmatched (no map pin): ", paste(head(unmatched, 4), collapse = ", "),
              if (length(unmatched) > 4) "…" else "", "."
            ))
          )
        }

        # Warn when K exceeds the number of populations (over-parameterised)
        if (n_k > n_meta_pops) {
          k_warn_msg <- div(
            class = "alert alert-warning",
            style = "font-size:0.8em; padding:6px 10px; margin:6px 0 0;",
            icon("exclamation-triangle"), " ",
            HTML(paste0(
              "<strong>K = ", n_k, " exceeds the number of populations (", n_meta_pops, ").</strong> ",
              "Over-parameterised ADMIXTURE runs can produce unreliable ancestry assignments. ",
              "Consider comparing CV error across K values and using K ≤ number of populations."
            ))
          )
        }
      }

      is_demo_q <- !is.null(demo_q_active())
      tagList(
        div(
          class = if (is_demo_q) "alert alert-warning" else "alert alert-success",
          style = "font-size:0.8em; padding:6px 10px; margin:6px 0 0;",
          icon(if (is_demo_q) "flask" else "check-circle"), " ",
          if (is_demo_q)
            HTML(paste0(
              "<strong>Demo Q-matrix loaded (K=2, simulated).</strong> ",
              sprintf("%d individuals · %d populations. ", nrow(qd), length(q_pops)),
              "Upload your own Q-matrix above to replace it."
            ))
          else
            sprintf(
              "%d individuals · %d populations · K = %d loaded.",
              nrow(qd), length(q_pops), n_k
            )
        ),
        overlap_msg,
        k_warn_msg
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
      qd <- q_data()
      req(qd)
      k_cols <- grep("^K[0-9]+$", names(qd), value = TRUE)
      n_inds <- nrow(qd)
      n_pops <- length(unique(qd$population))
      k_val  <- length(k_cols)

      # Modal ancestry component: most common dominant K across individuals
      if (k_val > 0) {
        dominant_k <- apply(qd[, k_cols, drop = FALSE], 1,
                            function(r) k_cols[which.max(r)])
        modal_k    <- names(sort(table(dominant_k), decreasing = TRUE))[1]
      } else {
        modal_k <- "—"
      }

      div(
        class = "row g-2 mb-2",
        admix_card("Individuals",  n_inds,  "#1B3A5C", "users"),
        admix_card("Populations",  n_pops,  "#2E75B6", "layer-group"),
        admix_card("K components", k_val,   "#388E3C", "palette"),
        admix_card("Modal ancestry", modal_k, "#F57C00", "chart-pie")
      )
    })

    # ── STRUCTURE-style stacked bar chart ────────────────────────────────────
    output$structure_bar <- plotly::renderPlotly({
      qd     <- q_data()
      req(qd)
      k_cols <- grep("^K[0-9]+$", names(qd), value = TRUE)
      req(length(k_cols) > 0)

      # Sort individuals — input$sort_by is either "population" or "k1"/"k2"/...
      sort_col <- toupper(input$sort_by %||% "population")   # "k1" → "K1"
      if (sort_col %in% k_cols) {
        qd <- qd[order(qd$population, -qd[[sort_col]]), ]
      } else {
        qd <- qd[order(qd$population, qd$sample_id), ]
      }
      qd$x_idx <- seq_len(nrow(qd))

      # Build figure with one trace per K component
      fig <- plotly::plot_ly()
      for (ki in seq_along(k_cols)) {
        kcol    <- k_cols[ki]
        pal_idx <- ((ki - 1L) %% length(ADMIX_K_PALETTE)) + 1L
        fig <- plotly::add_trace(
          fig,
          x           = qd$x_idx,
          y           = qd[[kcol]],
          name        = kcol,
          type        = "bar",
          marker      = list(
            color = ADMIX_K_PALETTE[pal_idx],
            line  = list(width = 0)
          ),
          text        = paste0(
            qd$sample_id, "<br/>",
            "Pop: ", qd$population, "<br/>",
            kcol, ": ", round(qd[[kcol]], 3)
          ),
          hoverinfo   = "text"
        )
      }

      # Population boundary shapes (white vertical dividers)
      # In plotly with integer x-positions, bar k occupies [k-0.5, k+0.5].
      # idx is the first row index of a new population group, so the boundary
      # between groups falls at x = idx - 0.5 (between bar idx-1 and bar idx).
      pop_changes <- which(c(TRUE, diff(as.integer(factor(qd$population))) != 0))
      shapes <- lapply(pop_changes[-1], function(idx) {
        list(
          type  = "line",
          x0    = idx - 0.5, x1 = idx - 0.5,
          y0    = 0,         y1 = 1,
          xref  = "x",       yref = "y",
          line  = list(color = "white", width = 2)
        )
      })

      # Population label annotations (centred below each group)
      pop_rle    <- rle(as.character(qd$population))
      pop_ends   <- cumsum(pop_rle$lengths)
      pop_starts <- c(1L, pop_ends[-length(pop_ends)] + 1L)
      pop_mids   <- (pop_starts + pop_ends) / 2

      annotations <- lapply(seq_along(pop_rle$values), function(i) {
        lbl <- pop_rle$values[i]
        if (nchar(lbl) > 12) lbl <- paste0(substr(lbl, 1, 11), "…")
        list(
          x         = pop_mids[i],
          y         = -0.08,
          xref      = "x",
          yref      = "paper",
          text      = lbl,
          showarrow = FALSE,
          font      = list(size = 9, color = "#555"),
          xanchor   = "center",
          yanchor   = "top"
        )
      })

      fig |>
        plotly::layout(
          barmode     = "stack",
          bargap      = 0,
          xaxis       = list(
            showticklabels = FALSE,
            showgrid       = FALSE,
            zeroline       = FALSE,
            title          = "",
            range          = c(0.5, nrow(qd) + 0.5)
          ),
          yaxis       = list(
            title    = "Ancestry proportion",
            range    = c(0, 1),
            showgrid = FALSE,
            tickvals = c(0, 0.25, 0.5, 0.75, 1)
          ),
          shapes      = shapes,
          annotations = annotations,
          showlegend  = FALSE,
          paper_bgcolor = "white",
          plot_bgcolor  = "white",
          font          = list(
            family = "Helvetica Neue, Helvetica, Arial, sans-serif",
            size   = 11
          ),
          margin = list(l = 50, r = 10, t = 10, b = 60)
        ) |>
        plotly::config(
          displayModeBar         = TRUE,
          modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
          displaylogo            = FALSE,
          toImageButtonOptions   = list(
            format   = "png",
            scale    = 2,
            filename = "popgen_admixture_bar"
          )
        )
    })

    # ── K colour legend ──────────────────────────────────────────────────────
    output$k_legend <- renderUI({
      qd <- q_data()
      req(qd)
      k_cols <- grep("^K[0-9]+$", names(qd), value = TRUE)
      req(length(k_cols) > 0)

      swatches <- lapply(seq_along(k_cols), function(ki) {
        pal_idx <- ((ki - 1L) %% length(ADMIX_K_PALETTE)) + 1L
        div(
          class = "d-flex align-items-center gap-1 me-2",
          tags$div(style = sprintf(
            "width:12px;height:12px;border-radius:2px;background:%s;flex-shrink:0;",
            ADMIX_K_PALETTE[pal_idx]
          )),
          tags$span(k_cols[ki], style = "font-size:0.8em;")
        )
      })

      div(
        class = "d-flex flex-wrap mt-2",
        style = "gap:4px;",
        do.call(tagList, swatches)
      )
    })

    # ── CSV download ─────────────────────────────────────────────────────────
    output$dl_qmat <- downloadHandler(
      filename = function() paste0("popgen_qmatrix_", Sys.Date(), ".csv"),
      content  = function(file) {
        qd <- q_data()
        req(qd)
        write.csv(qd, file, row.names = FALSE)
      }
    )

    # ── Reactive: per-population mean Q (consumed by mapServer) ─────────────
    pop_q <- reactive({
      qd <- q_data()
      if (is.null(qd)) return(NULL)
      k_cols <- grep("^K[0-9]+$", names(qd), value = TRUE)
      if (length(k_cols) == 0) return(NULL)
      # Average ancestry proportions within each population
      agg <- stats::aggregate(
        qd[, k_cols, drop = FALSE],
        by    = list(population = qd$population),
        FUN   = mean,
        na.rm = TRUE
      )
      agg
    })

    # ── Update sort options when Q-matrix loads ───────────────────────────────
    observeEvent(q_data(), {
      qd <- q_data()
      req(qd)
      k_cols <- grep("^K[0-9]+$", names(qd), value = TRUE)
      k_choices <- stats::setNames(
        tolower(k_cols),                          # "k1", "k2", …
        paste0(k_cols, " ancestry proportion")    # display label
      )
      updateSelectInput(session, "sort_by",
                        choices  = c("Population" = "population", k_choices),
                        selected = "population")
    })

    # Return pop_q directly — it is already a reactive, no wrapping needed
    pop_q
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Parse a Q-matrix CSV into a standardised data frame
#'
#' Accepts several input formats:
#'   Format A: sample_id, population, K1, K2, ..., Kn  (preferred)
#'   Format B: sample_id, K1, K2, ..., Kn              (no population column)
#'   Format C: pure numeric matrix (ADMIXTURE native .Q, space-separated)
#'
#' @param filepath Path to the uploaded file
#' @return data.frame with columns: sample_id, population, K1, ..., Kn
parse_qmatrix <- function(filepath) {

  # Try CSV first, then fall back to space-separated
  df <- tryCatch(
    utils::read.csv(filepath, stringsAsFactors = FALSE,
                    check.names = FALSE, strip.white = TRUE),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
    df <- tryCatch(
      utils::read.table(filepath, header = FALSE,
                        stringsAsFactors = FALSE, fill = TRUE),
      error = function(e) stop("Cannot parse file as CSV or whitespace-separated table.")
    )
  }

  if (nrow(df) == 0 || ncol(df) == 0) stop("File is empty.")

  col_lc <- tolower(trimws(colnames(df)))

  # Locate meta columns
  sample_idx <- which(col_lc %in%
    c("sample_id", "sample", "id", "individual", "ind", "indv", "name"))[1]
  pop_idx    <- which(col_lc %in%
    c("population", "pop", "group", "cluster", "locality", "site"))[1]

  has_sample <- !is.na(sample_idx)
  has_pop    <- !is.na(pop_idx)

  # K columns = numeric columns excluding meta
  meta_idx <- c(sample_idx, pop_idx)
  meta_idx <- meta_idx[!is.na(meta_idx)]
  k_idx    <- setdiff(seq_len(ncol(df)), meta_idx)

  # Keep only columns that are numeric or can be coerced to numeric
  k_idx <- k_idx[vapply(k_idx, function(i) {
    v <- df[[i]]
    is.numeric(v) || !any(is.na(suppressWarnings(as.numeric(v))))
  }, logical(1))]

  if (length(k_idx) == 0) {
    stop("No numeric K columns found. Check that ancestry proportion columns contain numbers.")
  }

  # Coerce K columns to numeric
  for (i in k_idx) df[[i]] <- as.numeric(df[[i]])

  # Build standardised output
  out <- data.frame(
    sample_id  = if (has_sample) as.character(df[[sample_idx]])
                 else paste0("Ind", seq_len(nrow(df))),
    population = if (has_pop)    as.character(df[[pop_idx]])
                 else "Unknown",
    stringsAsFactors = FALSE
  )

  k_data        <- df[, k_idx, drop = FALSE]
  names(k_data) <- paste0("K", seq_along(k_idx))
  out           <- cbind(out, k_data)

  # Attach a row-sum warning as an attribute rather than calling warning().
  # In Shiny, warning() goes to the R console and is invisible to users.
  # The caller (q_data reactive) reads this attribute and fires showNotification.
  new_k_cols <- paste0("K", seq_along(k_idx))
  row_sums   <- rowSums(out[, new_k_cols, drop = FALSE], na.rm = TRUE)
  if (any(abs(row_sums - 1) > 0.1, na.rm = TRUE)) {
    attr(out, "rowsum_warning") <- sprintf(
      "Q-matrix row sums range %.3f–%.3f (expected ≈1.0). ",
      min(row_sums, na.rm = TRUE), max(row_sums, na.rm = TRUE)
    )
  }

  out
}

#' Generate an SVG pie chart as a URL-encoded data URI
#'
#' Produces a compact SVG string with one path per slice, then URL-encodes
#' it for safe embedding as a Leaflet icon URL (no base64 dependency needed).
#'
#' @param props   Numeric vector of proportions (will be normalised to sum = 1)
#' @param colours Character vector of hex/named colours, one per element of props
#' @param size    Integer pixel diameter of the SVG canvas (default 40)
#' @return data URI string, or NULL if props is empty/zero
svg_pie_uri <- function(props, colours, size = 40L) {
  n     <- length(props)
  total <- sum(props, na.rm = TRUE)
  if (n == 0L || total == 0) return(NULL)

  # Defensive: recycle colours to match props length so colours[i] never returns NA
  if (length(colours) < n) colours <- rep_len(colours, n)

  props <- props / total  # normalise
  r     <- size / 2
  cx    <- r
  cy    <- r

  # Build one SVG element per slice
  paths <- vapply(seq_len(n), function(i) {
    start_ang <- sum(props[seq_len(i - 1L)]) * 2 * pi - pi / 2
    end_ang   <- start_ang + props[i] * 2 * pi

    # Full circle: use <circle> to avoid degenerate M-L-A-Z
    if (props[i] >= 0.9999) {
      return(sprintf(
        "<circle cx='%.2f' cy='%.2f' r='%.2f' fill='%s'/>",
        cx, cy, r, colours[i]
      ))
    }

    x1 <- cx + r * cos(start_ang)
    y1 <- cy + r * sin(start_ang)
    x2 <- cx + r * cos(end_ang)
    y2 <- cy + r * sin(end_ang)
    large_arc <- if (props[i] > 0.5) 1L else 0L

    # SVG arc: sweep-flag=1 → clockwise (screen coordinates)
    sprintf(
      "<path d='M %.2f %.2f L %.2f %.2f A %.2f %.2f 0 %d 1 %.2f %.2f Z' fill='%s'/>",
      cx, cy, x1, y1, r, r, large_arc, x2, y2, colours[i]
    )
  }, character(1))

  border <- sprintf(
    "<circle cx='%.2f' cy='%.2f' r='%.2f' fill='none' stroke='white' stroke-width='1.5'/>",
    cx, cy, r
  )

  svg <- sprintf(
    "<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d'>%s%s</svg>",
    as.integer(size), as.integer(size),
    paste(paths, collapse = ""),
    border
  )

  # Minimal URL encoding — only characters that break data URIs in browsers
  svg <- gsub(" ",  "%20", svg, fixed = TRUE)
  svg <- gsub("#",  "%23", svg, fixed = TRUE)
  svg <- gsub("+",  "%2B", svg, fixed = TRUE)

  paste0("data:image/svg+xml,", svg)
}

#' Draw per-population ADMIXTURE pie chart markers on a Leaflet proxy
#'
#' Clears the "admixture_pies" group and redraws one marker per population
#' that has matching coordinates in df_coords.
#'
#' @param proxy     leafletProxy object (already bound to session)
#' @param pop_q     Data frame with columns: population, K1, K2, ..., Kn
#'                  (per-population mean ancestry from admixtureServer)
#' @param df_coords pop_stats data frame (columns: population, lat, lon)
#' @param k_colours Character vector of hex colours, one per K component
#' @param size      Integer pixel diameter of each pie marker (default 40)
draw_admixture_pies <- function(proxy, pop_q, df_coords, k_colours, size = 40L) {
  proxy <- proxy |> clearGroup("admixture_pies")

  k_cols <- grep("^K[0-9]+$", names(pop_q), value = TRUE)
  if (length(k_cols) == 0) return(invisible(proxy))

  # Merge Q proportions with coordinates (inner join — only matched populations)
  coord_sub <- df_coords[, c("population", "lat", "lon")]
  merged    <- merge(pop_q, coord_sub, by = "population", all.x = FALSE)
  if (nrow(merged) == 0) return(invisible(proxy))

  half <- size / 2L

  for (i in seq_len(nrow(merged))) {
    props  <- as.numeric(merged[i, k_cols])
    uri    <- svg_pie_uri(props, k_colours[seq_along(k_cols)], size = size)
    if (is.null(uri)) next

    pop_nm <- merged$population[i]

    # Popup: population name + K breakdown table
    k_rows <- paste(
      sprintf(
        "<tr><td style='padding:2px 6px;'>%s</td><td style='padding:2px 6px; text-align:right;'>%.1f%%</td></tr>",
        k_cols, props * 100
      ),
      collapse = ""
    )
    popup_html <- sprintf(
      "<strong>%s</strong><br/><table style='font-size:0.85em;margin-top:4px;'>%s</table>",
      htmltools::htmlEscape(pop_nm), k_rows
    )

    proxy <- proxy |> addMarkers(
      lng   = merged$lon[i],
      lat   = merged$lat[i],
      icon  = leaflet::makeIcon(
        iconUrl     = uri,
        iconWidth   = size,
        iconHeight  = size,
        iconAnchorX = half,
        iconAnchorY = half
      ),
      popup = popup_html,
      label = htmltools::HTML(sprintf(
        "<strong>%s</strong><br/>ADMIXTURE pies (K=%d)",
        htmltools::htmlEscape(pop_nm), length(k_cols)
      )),
      labelOptions = labelOptions(
        style     = list("font-weight" = "normal", padding = "3px 6px"),
        textsize  = "12px",
        direction = "auto"
      ),
      group = "admixture_pies"
    )
  }

  invisible(proxy)
}

#' Individual summary card for the ADMIXTURE panel
admix_card <- function(label, value, colour, icon_name) {
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
