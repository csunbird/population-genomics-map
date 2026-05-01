# =============================================================================
# mod_map.R — Interactive Leaflet map module
# Renders population markers coloured by selected genetic diversity metric,
# plus optional overlays:
#   - Pairwise FST edge network (Phase 2 Goal 2)
#   - ADMIXTURE pie chart markers (Phase 2 Goal 4)
#
# Design decisions:
#  - Data-load trigger: updates markers AND fits map bounds to the new dataset
#  - Metric/size change trigger: updates marker colours/radii ONLY (no re-fit)
#    so the user's manual pan/zoom is preserved
#  - He uses the IUCN-calibrated categorical colour scale from genomics.R
#  - Other metrics use continuous diverging/sequential palettes
#  - FST network: categorical Wright-threshold colours; one addPolylines call per
#    edge (acceptable for typical conservation datasets of 5–30 populations)
#  - ADMIXTURE pies: SVG data-URI icons via leaflet::makeIcon; pie toggle is
#    hidden until a Q-matrix is loaded
#  - Legend shows a colour-coded key for the active metric, plus overlay legends
#    when FST network and/or ADMIXTURE pies are enabled
# =============================================================================

#' Map module UI
#'
#' @param id Module namespace ID
mapUI <- function(id) {
  ns <- NS(id)

  tagList(
    # ── Controls bar above the map ────────────────────────────────────────
    div(
      class = "map-controls",
      div(
        class = "d-flex align-items-center gap-3 flex-wrap",

        div(
          style = "flex: 0 0 auto;",
          tags$label("Colour by:", class = "me-1 fw-semibold",
                     style = "font-size:0.9em;"),
          selectInput(
            ns("colour_metric"),
            label    = NULL,
            choices  = c(
              "Expected heterozygosity (He)"     = "He",
              "Observed heterozygosity (Ho)"     = "Ho",
              "Nucleotide diversity (\u03c0)"        = "pi",
              "Inbreeding coefficient (F)"       = "inbreeding_f",
              "Effective population size (Ne)"   = "Ne",
              "Sample size (n)"                  = "n"
            ),
            selected = "He",
            width    = "260px"
          )
        ),

        div(
          style = "flex: 0 0 auto;",
          tags$label("Marker size:", class = "me-1 fw-semibold",
                     style = "font-size:0.9em;"),
          selectInput(
            ns("size_metric"),
            label    = NULL,
            choices  = c(
              "Fixed size"         = "fixed",
              "Proportional to n"  = "n",
              "Proportional to He" = "He"
            ),
            selected = "fixed",
            width    = "200px"
          )
        ),

        div(
          style = "flex: 0 0 auto; margin-top: 22px;",
          actionButton(
            ns("fit_bounds"),
            label  = tagList(icon("compress-arrows-alt"), " Fit to data"),
            class  = "btn btn-outline-secondary btn-sm",
            title  = "Re-centre map on current populations"
          )
        ),

        # ── FST network toggle ───────────────────────────────────────────────
        div(
          style = "flex: 0 0 auto; margin-top: 22px;",
          checkboxInput(
            ns("show_fst_network"),
            label = HTML("Show F<sub>ST</sub> network"),
            value = FALSE
          )
        ),

        # ── FST threshold slider (hidden until checkbox is ticked) ───────────
        shinyjs::hidden(
          div(
            id    = ns("fst_threshold_wrap"),
            style = "flex: 0 0 auto;",
            tags$label(
              HTML("Min F<sub>ST</sub>:"),
              class = "me-1 fw-semibold",
              style = "font-size:0.9em;"
            ),
            sliderInput(
              ns("fst_threshold"),
              label = NULL,
              min   = 0,
              max   = 0.5,
              value = 0.05,
              step  = 0.01,
              width = "160px"
            )
          )
        ),

        # ── ADMIXTURE pie toggle (hidden until Q-matrix is loaded) ───────────
        shinyjs::hidden(
          div(
            id    = ns("admix_pie_wrap"),
            style = "flex: 0 0 auto; margin-top: 22px;",
            checkboxInput(
              ns("show_admix_pies"),
              label = HTML("Show ADMIXTURE pies"),
              value = FALSE
            )
          )
        )
      )
    ),

    # ── Leaflet map ────────────────────────────────────────────────────────
    leafletOutput(ns("map"), height = "520px"),

    # ── Legend ─────────────────────────────────────────────────────────────
    uiOutput(ns("legend_ui"))
  )
}

#' Map module server
#'
#' @param id        Module namespace ID
#' @param pop_stats Reactive returning a data frame of per-population statistics
#'                  (columns: population, n, n_loci, Ho, He, pi, inbreeding_f, lat, lon)
#' @param fst_mat   Reactive returning the pairwise FST matrix (populations × populations)
#'                  computed by calc_fst_matrix(); may be NULL before data are loaded
#' @param pop_q     Optional reactive returning per-population mean Q data frame
#'                  (columns: population, K1, ..., Kn) from admixtureServer;
#'                  defaults to reactive(NULL) — pie overlay hidden when NULL
#' @param ne_stats  Optional reactive returning per-population Ne data frame
#'                  (columns: population, ne, ...) from neServer;
#'                  defaults to reactive(NULL) — Ne colour option greyed when NULL
mapServer <- function(id, pop_stats, fst_mat, pop_q = reactive(NULL),
                      ne_stats = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {

    # Track the population bounding box of the most recently loaded dataset
    # so the "Fit to data" button works independently of metric changes
    data_bounds   <- reactiveVal(NULL)

    # Pre-computed popup HTML — rebuilt only when the dataset changes, not on
    # every colour/size metric switch. Popup content (all stats) is independent
    # of which metric is currently displayed.
    cached_popups <- reactiveVal(NULL)

    # ── Base map (rendered once; subsequent updates use leafletProxy) ─────
    output$map <- renderLeaflet({
      m <- leaflet(options = leafletOptions(
        zoomControl   = TRUE,
        preferCanvas  = TRUE   # Canvas renderer: far fewer DOM operations during
                               # zoom/pan than the default SVG renderer
      )) |>
        addProviderTiles(
          "CartoDB.Positron",
          group   = "Light",
          options = providerTileOptions(
            updateWhenIdle   = TRUE,   # fetch new tiles only after panning stops
            updateWhenZooming = FALSE  # skip intermediate zoom frames
          )
        ) |>
        addProviderTiles(
          "Esri.WorldImagery",
          group   = "Satellite",
          options = providerTileOptions(
            opacity          = 0.75,
            updateWhenIdle   = TRUE,
            updateWhenZooming = FALSE,
            detectRetina     = FALSE   # halves tile requests on HiDPI screens
          )
        ) |>
        addLayersControl(
          baseGroups = c("Light", "Satellite"),
          options    = layersControlOptions(collapsed = FALSE)
        ) |>
        setView(lng = 0, lat = 20, zoom = 2)

      # ── Phase 5: screenshot button (PNG export) ──────────────────────────
      # leaflet.extras2::addScreenshotControl() adds a camera-icon button that
      # uses the leaflet-simple-map-screenshoter JS plugin to capture the current
      # map viewport as PNG — pure client-side, no server round-trip.
      if (requireNamespace("leaflet.extras2", quietly = TRUE)) {
        m <- m |> leaflet.extras2::addScreenshotControl(
          options = leaflet.extras2::screenshotControlOptions(
            position    = "topleft",
            title       = "Save map as PNG",
            iconUrl     = NULL,       # default camera icon
            hideElementsWithSelectors = ".leaflet-control-zoom",
            screenType  = "Canvas"    # matches preferCanvas = TRUE above
          )
        )
      }

      m
    })

    # ── Observe: new dataset loaded ───────────────────────────────────────
    # Fires when pop_stats() changes. Pre-builds popups, redraws markers, fits bounds.
    observeEvent(pop_stats(), {
      df <- pop_stats()
      req(df, nrow(df) > 0, "lat" %in% names(df), "lon" %in% names(df))

      # Build popup HTML once per dataset load — all metrics are always shown
      # in the popup, so this never needs to change when the user switches metric.
      cached_popups(vapply(seq_len(nrow(df)), function(i) {
        build_popup_html(df[i, ])
      }, character(1)))

      # Store bounds for re-use by the Fit button
      data_bounds(list(
        lng1 = min(df$lon, na.rm = TRUE) - 0.5,
        lat1 = min(df$lat, na.rm = TRUE) - 0.5,
        lng2 = max(df$lon, na.rm = TRUE) + 0.5,
        lat2 = max(df$lat, na.rm = TRUE) + 0.5
      ))

      redraw_markers(df, input$colour_metric %||% "He",
                     input$size_metric %||% "fixed", fit = TRUE)
    })

    # ── Observe: metric or size selector changed ──────────────────────────
    # Only recolours/resizes markers; does NOT move the map.
    observeEvent(
      list(input$colour_metric, input$size_metric),
      {
        df <- pop_stats()
        req(df, nrow(df) > 0)
        redraw_markers(df, input$colour_metric %||% "He",
                       input$size_metric %||% "fixed", fit = FALSE)
      },
      ignoreInit = TRUE   # don't double-fire on first load alongside the data observer
    )

    # ── Fit-to-data button ────────────────────────────────────────────────
    observeEvent(input$fit_bounds, {
      b <- data_bounds()
      req(b)
      leafletProxy("map", session) |>
        fitBounds(b$lng1, b$lat1, b$lng2, b$lat2)
    })

    # ── FST network: toggle + threshold + data changes ────────────────────
    # Single observer handles all three triggers so logic is in one place.
    observeEvent(
      list(input$show_fst_network, input$fst_threshold, fst_mat()),
      {
        show <- isTRUE(input$show_fst_network)
        shinyjs::toggle(id = "fst_threshold_wrap", condition = show)

        proxy <- leafletProxy("map", session)

        if (!show) {
          proxy |> clearGroup("fst_network")
          return()
        }

        df  <- pop_stats()
        mat <- fst_mat()
        req(df, mat, nrow(df) >= 2, nrow(mat) >= 2)
        draw_fst_network(proxy, df, mat, input$fst_threshold %||% 0.05)
      },
      ignoreInit = TRUE
    )

    # ── ADMIXTURE pies: reveal toggle when Q-matrix loads / unloads ─────────
    observeEvent(pop_q(), {
      pq <- pop_q()
      shinyjs::toggle(id = "admix_pie_wrap", condition = !is.null(pq))
      # Auto-clear pies if the Q-matrix is removed
      if (is.null(pq)) {
        leafletProxy("map", session) |> clearGroup("admixture_pies")
      }
    })

    # ── ADMIXTURE pies: draw or clear when toggle / Q-data changes ───────────
    observeEvent(
      list(input$show_admix_pies, pop_q()),
      {
        show <- isTRUE(input$show_admix_pies)
        proxy <- leafletProxy("map", session)

        if (!show) {
          proxy |> clearGroup("admixture_pies")
          return()
        }

        pq <- pop_q()
        df <- pop_stats()
        req(pq, df)

        k_cols    <- grep("^K[0-9]+$", names(pq), value = TRUE)
        k_colours <- ADMIX_K_PALETTE[
          ((seq_along(k_cols) - 1L) %% length(ADMIX_K_PALETTE)) + 1L
        ]
        draw_admixture_pies(proxy, pq, df, k_colours)
      },
      ignoreInit = TRUE
    )

    # ── Ne data arrives / changes: refresh markers if Ne is the active metric ──
    observeEvent(ne_stats(), {
      df <- pop_stats()
      req(df, nrow(df) > 0, isTRUE(input$colour_metric == "Ne"))
      redraw_markers(df, "Ne", input$size_metric %||% "fixed", fit = FALSE)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # ── Legend ─────────────────────────────────────────────────────────────
    output$legend_ui <- renderUI({
      df     <- pop_stats()
      metric <- input$colour_metric %||% "He"

      main_legend <- if (is.null(df)) {
        div(
          class = "legend-box",
          style = "color:#888; font-size:0.85em;",
          "Load data to see the map legend."
        )
      } else {
        build_legend_ui(metric, df[[metric]])
      }

      # Collect optional overlay legends
      overlay_legends <- tagList()

      # FST edge legend strip
      if (isTRUE(input$show_fst_network) && !is.null(df)) {
        overlay_legends <- tagList(
          overlay_legends,
          div(
            class = "legend-box mt-2",
            tags$strong(HTML("F<sub>ST</sub> network"), style = "font-size:0.85em;"),
            tags$br(),
            div(
              class = "d-flex gap-3 flex-wrap mt-1",
              fst_edge_swatch("#90CAF9", "< 0.05"),
              fst_edge_swatch("#1565C0", "0.05–0.15"),
              fst_edge_swatch("#E65100", "0.15–0.25"),
              fst_edge_swatch("#B71C1C", "> 0.25")
            ),
            tags$p(
              style = "font-size:0.75em; color:#777; margin: 3px 0 0;",
              "little | moderate | great | very great differentiation"
            )
          )
        )
      }

      # ADMIXTURE K colour legend
      pq <- pop_q()
      if (isTRUE(input$show_admix_pies) && !is.null(pq)) {
        k_cols    <- grep("^K[0-9]+$", names(pq), value = TRUE)
        k_swatches <- lapply(seq_along(k_cols), function(ki) {
          pal_idx <- ((ki - 1L) %% length(ADMIX_K_PALETTE)) + 1L
          div(
            class = "d-flex align-items-center gap-1 me-2",
            tags$div(style = sprintf(
              "width:12px;height:12px;border-radius:50%%;background:%s;flex-shrink:0;",
              ADMIX_K_PALETTE[pal_idx]
            )),
            tags$span(k_cols[ki], style = "font-size:0.78em;")
          )
        })
        overlay_legends <- tagList(
          overlay_legends,
          div(
            class = "legend-box mt-2",
            tags$strong("ADMIXTURE ancestry (per-population mean)", style = "font-size:0.85em;"),
            tags$br(),
            div(
              class = "d-flex flex-wrap mt-1",
              style = "gap:4px;",
              do.call(tagList, k_swatches)
            ),
            tags$p(
              style = "font-size:0.74em; color:#777; margin:3px 0 0;",
              "Pie markers show mean Q; diversity circles are beneath them. Click a pie for K breakdown."
            )
          )
        )
      }

      tagList(main_legend, overlay_legends)
    })

    # ── Internal helper: redraw all circle markers via leafletProxy ───────
    redraw_markers <- function(df, metric, sizing, fit = FALSE) {

      # For Ne, join Ne estimates onto df (Ne is not in pop_stats; it comes from
      # the Ne module reactive).  The column is stored as "Ne" (capital) to match
      # the metric key used throughout.
      if (metric == "Ne") {
        ne_df <- ne_stats()
        if (!is.null(ne_df) && "ne" %in% names(ne_df)) {
          ne_lookup <- ne_df[, c("population", "ne"), drop = FALSE]
          df$Ne <- ne_lookup$ne[match(df$population, ne_lookup$population)]
        } else {
          df$Ne <- NA_real_
        }
      }

      vals   <- df[[metric]]
      cols   <- colour_for_metric(metric, vals)
      radii  <- marker_radius(sizing, df, metric)

      # For the Ne metric the df now has a "Ne" column that build_popup_html uses
      # to add an Ne row — so we must bypass the cache and rebuild fresh.
      # For all other metrics, re-use the pre-computed cache (fast path).
      popups <- if (metric == "Ne" && "Ne" %in% names(df)) {
        vapply(seq_len(nrow(df)), function(i) build_popup_html(df[i, ]), character(1))
      } else {
        cached_popups() %||% vapply(seq_len(nrow(df)), function(i) {
          build_popup_html(df[i, ])
        }, character(1))
      }

      # Hover labels DO depend on the active metric — kept lightweight (one
      # sprintf per population, no HTML table).
      # For Ne: "—" is too ambiguous; distinguish "not estimated" from data.
      labels <- lapply(seq_len(nrow(df)), function(i) {
        val_str <- if (metric == "Ne" && is.na(vals[i])) {
          "not estimated"
        } else {
          fmt(vals[i])
        }
        htmltools::HTML(sprintf(
          "<strong>%s</strong><br/>%s: %s",
          htmltools::htmlEscape(df$population[i]),
          metric_label(metric),
          val_str
        ))
      })

      # clearGroup("pop_circles") removes only the diversity CircleMarker layer.
      # clearMarkers() would incorrectly remove addMarkers() layers (ADMIXTURE pies)
      # and does NOT remove addCircleMarkers() layers (which are Shapes/Paths in
      # Leaflet.js), causing circles to stack on every redraw.
      proxy <- leafletProxy("map", session) |>
        clearGroup("pop_circles") |>
        addCircleMarkers(
          data         = df,
          lng          = ~lon,
          lat          = ~lat,
          radius       = radii,
          fillColor    = cols,
          fillOpacity  = 0.88,
          color        = "white",
          weight       = 2,
          opacity      = 1,
          popup        = popups,
          label        = labels,
          labelOptions = labelOptions(
            style     = list("font-weight" = "normal", padding = "4px 8px"),
            textsize  = "13px",
            direction = "auto"
          ),
          group        = "pop_circles"
        )

      if (fit) {
        b <- data_bounds()
        if (!is.null(b)) {
          proxy <- proxy |> fitBounds(b$lng1, b$lat1, b$lng2, b$lat2)
        }
      }

      proxy
    }
  })
}

# =============================================================================
# Helpers (module-private)
# =============================================================================

#' Map metric values to hex colours
#'
#' @param metric Column name string
#' @param vals   Numeric vector of values
#' @return Character vector of hex colour strings, same length as vals
colour_for_metric <- function(metric, vals) {

  # Replace NA with median for palette scaling (but keep NA labels in legend)
  safe_vals <- vals
  if (any(is.na(safe_vals))) {
    med <- median(safe_vals, na.rm = TRUE)
    safe_vals[is.na(safe_vals)] <- if (is.na(med)) 0 else med
  }

  if (metric == "He") {
    # IUCN-calibrated categorical scale — defined in genomics.R
    return(vapply(vals, he_colour, character(1)))
  }

  if (metric == "Ne") {
    # Categorical Ne conservation scale — defined in genomics.R
    return(vapply(vals, ne_colour, character(1)))
  }

  if (metric == "Ho") {
    # Sequential: low Ho = amber, high Ho = blue
    # Domain extends above 0.5 if needed — Ho > 0.5 is valid in outbred populations
    max_ho <- max(c(0.5, safe_vals), na.rm = TRUE)
    pal <- scales::col_numeric(
      palette  = c("#FBC02D", "#43A047", "#1565C0"),
      domain   = c(0, max_ho),
      na.color = "#AAAAAA"
    )
    return(pal(safe_vals))
  }

  if (metric == "pi") {
    # Sequential: pi ≈ He for large n; for small n the Nei correction pushes pi slightly > He
    # Domain extends above 0.5 to avoid clipping (e.g. n=5 correction factor = 10/9)
    max_pi <- max(c(0.5, safe_vals), na.rm = TRUE)
    pal <- scales::col_numeric(
      palette  = c("#FBC02D", "#43A047", "#1565C0"),
      domain   = c(0, max_pi),
      na.color = "#AAAAAA"
    )
    return(pal(safe_vals))
  }

  if (metric == "inbreeding_f") {
    # Diverging: negative F (outbreeding) = teal, 0 = white/light, positive = red
    # F is bounded [-1, 1] theoretically; practically [-0.5, 0.5] in real data
    pal <- scales::col_numeric(
      palette  = c("#00796B", "#F5F5F5", "#D32F2F"),
      domain   = c(-0.5, 0.5),
      na.color = "#AAAAAA"
    )
    return(pal(pmax(pmin(safe_vals, 0.5), -0.5)))  # clamp to display range
  }

  if (metric == "n") {
    # Sequential: low n = grey-blue, high n = deep blue
    rng <- range(safe_vals, na.rm = TRUE)
    # Guard against degenerate domain (all populations same size)
    if (rng[1] == rng[2]) return(rep("#1565C0", length(vals)))
    pal <- scales::col_numeric(
      palette  = c("#90A4AE", "#1565C0"),
      domain   = rng,
      na.color = "#AAAAAA"
    )
    return(pal(safe_vals))
  }

  # Fallback
  rep("#1565C0", length(vals))
}

#' Compute circle marker radii
#'
#' @param sizing "fixed" | "n" | "He"
#' @param df     pop_stats data frame
#' @param metric Active colour metric (unused here but kept for signature symmetry)
#' @return Numeric vector of radii in pixels
marker_radius <- function(sizing, df, metric) {
  n_pop <- nrow(df)
  if (sizing == "fixed") return(rep(14L, n_pop))

  vals <- df[[sizing]]
  if (all(is.na(vals))) return(rep(14L, n_pop))

  # Replace NA with min so those populations get the smallest marker
  vals[is.na(vals)] <- min(vals, na.rm = TRUE)

  rng <- range(vals, na.rm = TRUE)
  if (rng[1] == rng[2]) return(rep(14L, n_pop))   # all same → fixed

  as.integer(round(scales::rescale(vals, to = c(8, 24))))
}

#' Build the HTML legend panel below the map
#'
#' @param metric Active colour metric
#' @param vals   Numeric vector (used for range display on non-He metrics)
build_legend_ui <- function(metric, vals) {

  if (metric == "He") {
    pal_info <- he_legend_colours()
    swatches <- mapply(function(col, lbl) {
      tags$div(
        class = "d-flex align-items-center gap-2 mb-1",
        tags$div(style = sprintf(
          "width:14px;height:14px;border-radius:50%%;background:%s;flex-shrink:0;", col
        )),
        tags$span(lbl, style = "font-size:0.82em;")
      )
    }, pal_info$colours, pal_info$labels, SIMPLIFY = FALSE)

    return(div(
      class = "legend-box",
      tags$strong("Expected He (IUCN thresholds)", style = "font-size:0.85em;"),
      tags$br(),
      do.call(tagList, swatches)
    ))
  }

  if (metric == "Ne") {
    ne_swatch <- function(colour, label) {
      tags$div(
        class = "d-flex align-items-center gap-2 mb-1",
        tags$div(style = sprintf(
          "width:14px;height:14px;border-radius:50%%;background:%s;flex-shrink:0;", colour
        )),
        tags$span(label, style = "font-size:0.82em;")
      )
    }
    return(div(
      class = "legend-box",
      tags$strong(HTML("Effective size Ne (IUCN/SSC)"), style = "font-size:0.85em;"),
      tags$br(),
      ne_swatch("#D32F2F", "< 50 — Critical: immediate extinction risk"),
      ne_swatch("#F57C00", "50–99 — High risk: rapid genetic drift"),
      ne_swatch("#FBC02D", "100–499 — Concern: long-term viability at risk"),
      ne_swatch("#1565C0", "≥ 500 — Viable"),
      ne_swatch("#AAAAAA", "Not estimated / insufficient data")
    ))
  }

  # For continuous metrics: render a gradient bar + min/max labels
  valid <- vals[!is.na(vals)]
  if (length(valid) == 0) return(NULL)

  gradient_css <- switch(metric,
    Ho           = "linear-gradient(to right, #FBC02D, #43A047, #1565C0)",
    pi           = "linear-gradient(to right, #FBC02D, #43A047, #1565C0)",
    inbreeding_f = "linear-gradient(to right, #00796B, #F5F5F5, #D32F2F)",
    n            = "linear-gradient(to right, #90A4AE, #1565C0)",
    "linear-gradient(to right, #FBC02D, #1565C0)"
  )

  div(
    class = "legend-box",
    tags$strong(metric_label(metric), style = "font-size:0.85em;"),
    tags$br(),
    tags$div(
      style = sprintf(
        "width:100%%;height:12px;border-radius:6px;background:%s;margin:4px 0 2px;",
        gradient_css
      )
    ),
    tags$div(
      class = "d-flex justify-content-between",
      tags$span(fmt(min(valid)), style = "font-size:0.78em; color:#555;"),
      tags$span(fmt(max(valid)), style = "font-size:0.78em; color:#555;")
    ),
    if (metric == "inbreeding_f") {
      tags$p(
        style = "font-size:0.78em; color:#777; margin-top:4px;",
        "Teal = outbreeding (F < 0) | Red = inbreeding (F > 0)"
      )
    }
  )
}

#' Human-readable metric label
metric_label <- function(metric) {
  switch(metric,
    He           = "He",
    Ho           = "Ho",
    pi           = "\u03c0",
    inbreeding_f = "F",
    Ne           = "Ne",
    n            = "n",
    metric
  )
}

#' Draw FST edge network on a Leaflet proxy map
#'
#' Adds one polyline per population pair where FST \u2265 threshold.
#' Edges are coloured by Wright (1978) FST category (matching the Structure tab),
#' with weight and visibility increasing with FST value.
#'
#' @param proxy     leafletProxy object (already bound to session)
#' @param df        pop_stats data frame (needs columns: population, lat, lon)
#' @param mat       Named pairwise FST matrix from calc_fst_matrix()
#' @param threshold Minimum FST to display (default 0.05 = little differentiation)
draw_fst_network <- function(proxy, df, mat, threshold = 0.05) {
  proxy <- proxy |> clearGroup("fst_network")

  pops  <- rownames(mat)
  pairs <- which(upper.tri(mat), arr.ind = TRUE)
  if (nrow(pairs) == 0) return(invisible(proxy))

  MAX_EDGES <- 300L   # guard against very large datasets
  n_drawn   <- 0L

  for (k in seq_len(nrow(pairs))) {
    fst_val <- mat[pairs[k, 1], pairs[k, 2]]
    if (is.na(fst_val) || fst_val < threshold) next

    pop_i <- pops[pairs[k, 1]]
    pop_j <- pops[pairs[k, 2]]
    row_i <- df[df$population == pop_i, ]
    row_j <- df[df$population == pop_j, ]
    if (nrow(row_i) == 0 || nrow(row_j) == 0) next

    # Categorical colour and weight by Wright (1978) FST threshold
    edge_col <- dplyr::case_when(
      fst_val < FST_THRESHOLDS$little   ~ "#90CAF9",   # light blue  \u2014 little diff
      fst_val < FST_THRESHOLDS$moderate ~ "#1565C0",   # dark blue   \u2014 moderate
      fst_val < FST_THRESHOLDS$great    ~ "#E65100",   # orange      \u2014 great
      TRUE                              ~ "#B71C1C"    # dark red    \u2014 very great
    )
    edge_wt <- dplyr::case_when(
      fst_val < FST_THRESHOLDS$little   ~ 1.5,
      fst_val < FST_THRESHOLDS$moderate ~ 2.5,
      fst_val < FST_THRESHOLDS$great    ~ 3.5,
      TRUE                              ~ 4.5
    )

    popup_txt <- sprintf(
      "<strong>%s &#x2194; %s</strong><br/>F<sub>ST</sub> = %s &mdash; %s",
      htmltools::htmlEscape(pop_i), htmltools::htmlEscape(pop_j),
      fmt(fst_val, 3), fst_label(fst_val)
    )

    proxy <- proxy |> addPolylines(
      lng          = c(row_i$lon[1], row_j$lon[1]),
      lat          = c(row_i$lat[1], row_j$lat[1]),
      weight       = edge_wt,
      color        = edge_col,
      opacity      = 0.75,
      group        = "fst_network",
      popup        = popup_txt,
      label        = htmltools::HTML(popup_txt),
      labelOptions = labelOptions(
        style     = list("font-weight" = "normal", padding = "3px 6px"),
        textsize  = "12px",
        direction = "auto"
      )
    )

    n_drawn <- n_drawn + 1L
    if (n_drawn >= MAX_EDGES) {
      showNotification(
        paste0(
          "FST network: showing first ", MAX_EDGES, " edges. ",
          "Increase the min FST threshold to reduce edge count."
        ),
        type = "warning", duration = 8
      )
      break
    }
  }

  invisible(proxy)
}

#' Small coloured line swatch + label for the FST edge legend
#'
#' @param colour  Hex colour string
#' @param label   Short text label shown beside the swatch
fst_edge_swatch <- function(colour, label) {
  div(
    class = "d-flex align-items-center gap-1",
    tags$div(style = sprintf(
      "width:22px; height:3px; background:%s; border-radius:2px; flex-shrink:0;",
      colour
    )),
    tags$span(label, style = "font-size:0.78em;")
  )
}
