# =============================================================================
# server.R — Main server logic
# popgen-map Phase 1
# =============================================================================

server <- function(input, output, session) {

  # ── Upload module → data ──────────────────────────────────────────────────
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

  # ── Switch to map tab automatically when data loads ──────────────────────
  observeEvent(upload_data(), {
    if (isTRUE(upload_data()$ready)) {
      updateTabsetPanel(session, "main_tabs", selected = "tab_map")
    }
  })
}
