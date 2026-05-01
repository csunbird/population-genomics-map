# =============================================================================
# server.R — Main server logic
# popgen-map Phase 5
# =============================================================================

server <- function(input, output, session) {

  # ── Upload module → genomic data ──────────────────────────────────────────
  upload_data <- uploadServer("upload")

  # ── Compute population statistics (reactive, re-runs when data changes) ──
  pop_stats <- reactive({
    d <- upload_data()
    req(d$ready)
    tryCatch(
      calc_pop_stats(d$gt, d$metadata),
      error = function(e) {
        showNotification(
          paste0("Statistics error: ", e$message),
          type     = "error",
          duration = 10
        )
        NULL
      }
    )
  })

  # ── Pairwise FST matrix (Phase 2) ────────────────────────────────────────
  # Computed once at the app level so both the map overlay AND the structure
  # heatmap share the same result without redundant computation.
  fst_mat <- reactive({
    d <- upload_data()
    req(d$ready)
    tryCatch(
      calc_fst_matrix(d$gt, d$metadata),
      error = function(e) {
        showNotification(
          paste0("FST error: ", e$message),
          type     = "error",
          duration = 10
        )
        NULL
      }
    )
  })

  # ── ADMIXTURE module (Phase 2 Goal 4) — returns pop_q for map overlay ───
  pop_q <- admixtureServer("admixture", upload_data, pop_stats)

  # ── Phase 3: F-ROH and Ne modules ───────────────────────────────────────
  frohServer("froh", upload_data)
  ne_stats <- neServer("ne", upload_data)   # returns reactive for map Ne overlay

  # ── Phase 4 Goal 1: Environmental data upload ────────────────────────────
  env_data <- envUploadServer("env_upload", upload_data)

  # ── Phase 4 Goal 2: RDA + LFMM2 GEA ─────────────────────────────────────
  gea_data <- geaServer("gea", upload_data, env_data)

  # ── Phase 4 Goal 3: Genomic offset / climate vulnerability ───────────────
  offset_data <- offsetServer("offset", upload_data, env_data, gea_data)

  # ── Map module ──────────────────────────────────────────────────────────
  mapServer("map", pop_stats, fst_mat, pop_q, ne_stats)

  # ── Stats panel module ──────────────────────────────────────────────────
  statsPanelServer(
    "stats",
    pop_stats    = pop_stats,
    is_demo      = reactive(isTRUE(upload_data()$is_demo)),
    scenario_info = reactive(upload_data()$scenario_info)
  )

  # ── Structure module (Phase 2 Goal 1: pairwise FST heatmap) ─────────────
  structureServer("structure", upload_data, pop_stats, fst_mat)

  # ── PCA module (Phase 2 Goal 3) ──────────────────────────────────────────
  pcaServer("pca", upload_data)

  # ── Phase 5: Export module ─────────────────────────────────────────────
  exportServer(
    "export",
    pop_stats      = pop_stats,
    fst_mat        = fst_mat,
    pop_q          = pop_q,
    ne_stats       = ne_stats,
    gea_data       = gea_data,
    offset_data    = offset_data,
    upload_data    = upload_data,
    parent_session = session   # needed so "Go to Map" can switch the top-level tab
  )

  # ── Help module (Phase 5 — static UI, no server logic) ───────────────────
  helpServer("help")

  # ── Switch to map tab automatically when data loads ──────────────────────
  observeEvent(upload_data(), {
    if (isTRUE(upload_data()$ready)) {
      updateTabsetPanel(session, "main_tabs", selected = "tab_map")
    }
  })
}
