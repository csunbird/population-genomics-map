# =============================================================================
# ui.R — Main user interface definition
# popgen-map Phase 4
# =============================================================================

# ── UI helpers ────────────────────────────────────────────────────────────────

#' Build one row of the He threshold table in the About tab
threshold_row <- function(colour, range, status, guidance) {
  tags$tr(
    tags$td(style = paste0("padding:5px 8px; border-bottom:1px solid #eee;
                            border-left: 4px solid ", colour, ";"),
            range),
    tags$td(style = "padding:5px 8px; border-bottom:1px solid #eee;", status),
    tags$td(style = "padding:5px 8px; border-bottom:1px solid #eee; color:#555;", guidance)
  )
}

# ── Main UI ───────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = tagList(
    tags$span(
      style = "font-weight: 700; letter-spacing: 0.02em;",
      icon("globe-africa", style = "margin-right: 6px; color: #2E75B6;"),
      APP_TITLE
    ),
    tags$span(
      style = "font-size: 0.65em; font-weight: 400; color: #777; margin-left: 12px;",
      APP_PHASE, " — v", APP_VERSION
    )
  ),
  theme = bs_theme(
    version    = 5,
    bootswatch = "flatly",
    primary    = "#2E75B6",
    # Aptos (Office 365 / Windows 11) → Calibri (Office) → Candara → Segoe UI
    # All are local system fonts — no network request, works offline.
    base_font  = font_collection("Aptos", "Calibri", "Candara",
                                 "'Segoe UI'", "Arial", "sans-serif"),
    # Consolas is the standard Windows monospace; Courier New as final fallback
    code_font  = font_collection("Consolas", "'Courier New'", "monospace")
  ),
  fillable = TRUE,

  # ── Custom CSS + shinyjs ────────────────────────────────────────────────────
  tags$head(
    tags$link(rel = "stylesheet", href = "styles.css")
  ),
  shinyjs::useShinyjs(),

  # ── Sidebar ─────────────────────────────────────────────────────────────────
  sidebar = sidebar(
    width = 310,
    open  = TRUE,

    # ── Genomic data input section ────────────────────────────────────────────
    tags$h6(
      "Data Input",
      style = "font-weight: 600; color: #1B3A5C; margin-bottom: 8px;"
    ),

    uploadUI("upload"),

    # ── Phase 4: Environmental data section ──────────────────────────────────
    envUploadUI("env_upload"),

    tags$hr(style = "margin: 16px 0 12px;"),
    tags$p(
      style = "font-size: 0.75em; color: #999; line-height: 1.5;",
      HTML(paste0(
        "<strong>Phase 1</strong>: nucleotide diversity (&pi;), ",
        "expected &amp; observed heterozygosity (He / Ho), inbreeding (F).<br/>",
        "<strong>Phase 2</strong>: pairwise F&#x209B;&#x209C; heatmap, FST network, PCA, ADMIXTURE pies.<br/>",
        "<strong>Phase 3</strong>: inbreeding (F<sub>ROH</sub>) and effective population size (Ne).<br/>",
        "<strong>Phase 4</strong>: adaptive potential — RDA, LFMM2 GEA, genomic offset.<br/>",
        "<strong>Phase 5</strong>: export HTML report, CSV, map PNG &amp; chart PNG (2&times; resolution)."
      ))
    )
  ),

  # ── Main panel ──────────────────────────────────────────────────────────────
  navset_tab(
    id = "main_tabs",

    # ── Tab 1: Map ──────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("map-marked-alt"), " Map"),
      value = "tab_map",
      div(
        class = "p-3",
        mapUI("map")
      )
    ),

    # ── Tab 2: Statistics table ──────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("table"), " Statistics"),
      value = "tab_stats",
      div(
        class = "p-3",
        statsPanelUI("stats")
      )
    ),

    # ── Tab 3: Structure (Phase 2) ───────────────────────────────────────────
    nav_panel(
      title = tagList(icon("project-diagram"), " Structure"),
      value = "tab_structure",
      div(
        class = "p-3",
        structureUI("structure")
      )
    ),

    # ── Tab 4: PCA (Phase 2 Goal 3) ──────────────────────────────────────────
    nav_panel(
      title = tagList(icon("project-diagram"), " PCA"),
      value = "tab_pca",
      div(
        class = "p-3",
        pcaUI("pca")
      )
    ),

    # ── Tab 5: ADMIXTURE (Phase 2 Goal 4) ────────────────────────────────────
    nav_panel(
      title = tagList(icon("chart-pie"), " ADMIXTURE"),
      value = "tab_admixture",
      div(
        class = "p-3",
        admixtureUI("admixture")
      )
    ),

    # ── Tab 6: F-ROH (Phase 3 Goal 1) ────────────────────────────────────────
    nav_panel(
      title = tagList(icon("dna"), " F-ROH"),
      value = "tab_froh",
      div(
        class = "p-3",
        frohUI("froh")
      )
    ),

    # ── Tab 7: Effective Population Size (Phase 3 Goal 2) ─────────────────────
    nav_panel(
      title = tagList(icon("users"), " Ne"),
      value = "tab_ne",
      div(
        class = "p-3",
        neUI("ne")
      )
    ),

    # ── Tab 8: GEA — Genotype-Environment Association (Phase 4 Goal 2) ────────
    nav_panel(
      title = tagList(icon("leaf"), " GEA"),
      value = "tab_gea",
      div(
        class = "p-3",
        geaUI("gea")
      )
    ),

    # ── Tab 9: Genomic Offset / Climate Vulnerability (Phase 4 Goal 3) ────────
    nav_panel(
      title = tagList(icon("thermometer-half"), " Offset"),
      value = "tab_offset",
      div(
        class = "p-3",
        offsetUI("offset")
      )
    ),

    # ── Tab 10: Export (Phase 5) ──────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("file-export"), " Export"),
      value = "tab_export",
      div(
        class = "p-3",
        style = "max-width: 760px;",
        exportUI("export")
      )
    ),

    # ── Tab 11: Help (Phase 5) ────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("question-circle"), " Help"),
      value = "tab_help",
      div(
        class = "p-3",
        helpUI("help")
      )
    ),

    # ── Tab 12: About ─────────────────────────────────────────────────────────
    nav_panel(
      title = tagList(icon("info-circle"), " About"),
      value = "tab_about",
      div(
        class = "p-4",
        style = "max-width: 720px;",

        tags$h4("About PopGen Map", style = "color: #1B3A5C;"),
        tags$p(
          "PopGen Map is an open-source, browser-accessible platform that turns ",
          "population genomic data (VCF files + sample coordinates) into interactive ",
          "geographic visualisations of genetic health metrics — designed for ",
          "conservation practitioners, not bioinformaticians."
        ),

        tags$h5("Phase 1 metrics", style = "color: #2E75B6; margin-top: 1.2em;"),
        tags$ul(
          tags$li(HTML("<strong>He</strong> — Expected heterozygosity: the probability that two alleles drawn at random from a population are different. Higher = more diverse.")),
          tags$li(HTML("<strong>Ho</strong> — Observed heterozygosity: the proportion of individuals actually carrying two different alleles at a locus.")),
          tags$li(HTML("<strong>&pi;</strong> — Nucleotide diversity: average number of pairwise nucleotide differences per site.")),
          tags$li(HTML("<strong>F</strong> — Inbreeding coefficient: (He &minus; Ho) / He. Positive values indicate excess homozygosity relative to Hardy-Weinberg expectations (inbreeding). This may be associated with inbreeding depression, though heavily bottlenecked populations that have purged deleterious alleles can show high F without overt fitness reduction."))
        ),

        tags$h5("Conservation thresholds (He)", style = "color: #2E75B6; margin-top: 1.2em;"),
        tags$p(style = "font-size:0.9em; color:#555;",
          HTML("Thresholds follow Frankham et al. (2014) and are consistent with IUCN genetic diversity criteria.
               <strong>Important:</strong> these thresholds were calibrated on microsatellite and allozyme data.
               SNP-based He filtered at MAF &ge; 0.05 will typically yield lower values than microsatellite He for the same population.
               Use these thresholds as indicative guidance, not hard cut-offs, when interpreting SNP data.")),
        tags$table(
          class = "about-threshold-table",
          style = "width:100%; border-collapse:collapse; font-size:0.88em; margin-bottom:0.5em;",
          tags$thead(tags$tr(
            tags$th(style="padding:5px 8px; background:#1B3A5C; color:white; border-radius:4px 0 0 0;", "He range"),
            tags$th(style="padding:5px 8px; background:#1B3A5C; color:white;", "Status"),
            tags$th(style="padding:5px 8px; background:#1B3A5C; color:white; border-radius:0 4px 0 0;", "Guidance")
          )),
          tags$tbody(
            threshold_row("#D32F2F", "< 0.10",       "\U0001f534 Critical",          "High extinction risk; urgent action required"),
            threshold_row("#F57C00", "0.10 – 0.15",  "\U0001f7e0 High risk",          "Low diversity; intervention likely needed"),
            threshold_row("#FBC02D", "0.15 – 0.20",  "\U0001f7e1 Moderate risk",      "Monitoring and genetic rescue assessment recommended"),
            threshold_row("#388E3C", "0.20 – 0.30",  "\U0001f7e2 Lower risk",         "Moderate diversity; routine monitoring"),
            threshold_row("#1565C0", "> 0.30",        "\U0001f535 Healthy",            "High diversity; reference population candidate")
          )
        ),

        tags$h5("Phase 3 metrics", style = "color: #2E75B6; margin-top: 1.2em;"),
        tags$ul(
          tags$li(HTML(paste0(
            "<strong>F<sub>ROH</sub></strong> — Runs of Homozygosity inbreeding coefficient: ",
            "proportion of the genome in consecutive homozygous SNP runs. ",
            "Captures <em>recent</em> inbreeding (last ~10 generations). ",
            "F<sub>ROH</sub> ≥ 0.25 ≈ full-sibling or parent-offspring mating; ",
            "≥ 0.125 ≈ half-sibling mating; ",
            "≥ 0.0625 ≈ first-cousin mating."
          ))),
          tags$li(HTML(paste0(
            "<strong>Ne</strong> — Effective population size: the size of an idealised random-mating ",
            "population that would experience the same rate of genetic drift as the observed population. ",
            "Estimated from the linkage disequilibrium (r²) among pairs of loci ",
            "(Waples &amp; Do 2008 single-sample method). ",
            "IUCN/SSC thresholds: Ne &lt; 50 = immediate extinction risk; ",
            "Ne &lt; 100 = high risk; Ne &lt; 500 = long-term viability concern."
          )))
        ),

        tags$h5("Phase 4 metrics", style = "color: #2E75B6; margin-top: 1.2em;"),
        tags$ul(
          tags$li(HTML(paste0(
            "<strong>Partial RDA</strong> — Redundancy Analysis: multivariate regression of ",
            "population allele frequencies on environmental variables, with geography partialled out. ",
            "Identifies which environmental gradients explain genetic variation beyond spatial patterns. ",
            "Top-loading loci are candidate adaptive variants."
          ))),
          tags$li(HTML(paste0(
            "<strong>LFMM2</strong> — Latent Factor Mixed Model (Caye et al. 2019): ",
            "per-locus test for genotype-environment association, controlling for ",
            "population structure via K latent factors. ",
            "P-values are calibrated by the genomic inflation factor (GIF) and corrected by ",
            "Benjamini-Hochberg FDR. Identifies single-locus signals of local adaptation."
          ))),
          tags$li(HTML(paste0(
            "<strong>Genomic offset</strong> — RDA-based climate vulnerability metric ",
            "(Fitzpatrick &amp; Keller 2015): Euclidean distance in constrained ordination space ",
            "between a population's current and projected future environmental position. ",
            "Larger values indicate a larger adaptive gap under climate change. ",
            "Always interpret alongside Ne and connectivity data."
          )))
        ),

        tags$h5("Input format", style = "color: #2E75B6; margin-top: 1.4em;"),
        tags$p("Upload a standard biallelic SNP VCF file (any variant caller: GATK, STACKS, FreeBayes, DeepVariant) and a sample metadata CSV with columns:"),
        tags$code("sample_id, population, latitude, longitude"),
        tags$p(
          style = "margin-top: 0.6em; font-size:0.85em; color:#555;",
          "For Phase 4 (GEA / offset), also upload a per-population environmental CSV with columns: ",
          tags$code("population, BIO1, BIO12, …"),
          " or use the WorldClim / CMIP6 auto-fetch option in the sidebar."
        ),

        tags$h5("Roadmap", style = "color: #2E75B6; margin-top: 1.4em;"),
        tags$table(
          style = "width:100%; border-collapse:collapse; font-size:0.85em;",
          tags$thead(tags$tr(
            tags$th(style="padding:4px 8px; background:#f1f3f5; color:#333;", "Phase"),
            tags$th(style="padding:4px 8px; background:#f1f3f5; color:#333;", "Features"),
            tags$th(style="padding:4px 8px; background:#f1f3f5; color:#333;", "Status")
          )),
          tags$tbody(
            tags$tr(
              tags$td(style="padding:4px 8px;","1"),
              tags$td(style="padding:4px 8px;","Diversity mapping (He, Ho, π, F)"),
              tags$td(style="padding:4px 8px; color:#388E3C;","✅ Done")
            ),
            tags$tr(
              tags$td(style="padding:4px 8px; background:#fafafa;","2"),
              tags$td(style="padding:4px 8px; background:#fafafa;","Population structure: pairwise Fₛₜ heatmap, FST network, PCA, ADMIXTURE pies"),
              tags$td(style="padding:4px 8px; background:#fafafa; color:#388E3C;","✅ Done")
            ),
            tags$tr(
              tags$td(style="padding:4px 8px;","3"),
              tags$td(style="padding:4px 8px;","Inbreeding (F-ROH) and effective population size (Ne)"),
              tags$td(style="padding:4px 8px; color:#388E3C;","✅ Done")
            ),
            tags$tr(
              tags$td(style="padding:4px 8px; background:#fafafa;","4"),
              tags$td(style="padding:4px 8px; background:#fafafa;","Adaptive potential: RDA, LFMM2 GEA, genomic offset"),
              tags$td(style="padding:4px 8px; background:#fafafa; color:#2E75B6;","✅ Done")
            ),
            tags$tr(
              tags$td(style="padding:4px 8px;","5"),
              tags$td(style="padding:4px 8px;","Export: HTML conservation report, statistics CSV, map PNG, chart PNG (2× res); in-app Help tab"),
              tags$td(style="padding:4px 8px; color:#388E3C;","✅ Done")
            )
          )
        ),

        tags$h5("Export & Documentation", style = "color: #2E75B6; margin-top: 1.4em;"),
        tags$ul(
          tags$li(HTML(paste0(
            "<strong>Conservation Report (HTML)</strong>: full narrative report with all metrics, ",
            "risk summary, FST matrix, Ne, GEA, and genomic offset — generated from the Export tab."
          ))),
          tags$li(HTML(paste0(
            "<strong>Statistics CSV</strong>: per-population diversity table for downstream analysis ",
            "in R, Python, or Excel — from the Export tab or the Statistics tab."
          ))),
          tags$li(HTML(paste0(
            "<strong>Map PNG</strong>: click the camera icon inside the Map tab to save the current ",
            "map viewport as a PNG image."
          ))),
          tags$li(HTML(paste0(
            "<strong>Chart PNG</strong>: hover over any Plotly chart (PCA, ADMIXTURE bar, ",
            "F-ROH, Ne, GEA biplot/Manhattan, Offset bar) and click the camera icon ",
            "in its toolbar to download as PNG at 2&times; (print-quality) resolution. ",
            "Note: the FST heatmap is an HTML table — use browser Print or screenshot to capture it."
          ))),
          tags$li(HTML(paste0(
            "<strong>User Guide</strong>: full documentation available at ",
            "<code>docs/USER_GUIDE.md</code> in the repository."
          )))
        ),

        tags$p(
          style = "margin-top: 1.4em; font-size:0.82em; color:#999;",
          "Open-source — MIT licence — v", APP_VERSION,
          " | ",
          tags$a("github.com/csunbird/population-genomics-map",
                 href   = "https://github.com/csunbird/population-genomics-map",
                 target = "_blank"),
          " | Contact: ",
          tags$a("seahcheyanne@gmail.com",
                 href = "mailto:seahcheyanne@gmail.com")
        )
      )
    )
  )
)
