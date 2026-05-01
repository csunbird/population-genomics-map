# =============================================================================
# mod_help.R — In-app Help tab (Phase 5)
# Accordion-style reference covering Quick Start, all phases, export, and
# troubleshooting.  No server logic — pure UI.
# =============================================================================

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Styled section heading inside an accordion panel
help_h <- function(text, color = "#2E75B6") {
  tags$h6(text,
    style = paste0("color:", color, "; font-weight:600;",
                   " margin-top:1.1em; margin-bottom:.4em;"))
}

#' Inline code span (prefixed to avoid global namespace collision)
help_ic <- function(...) tags$code(...)

#' Styled note / callout box
help_note <- function(..., icon_chr = "ℹ️", color = "#e8f4fd",
                       border = "#2E75B6") {
  tags$div(
    style = paste0(
      "background:", color, "; border-left:4px solid ", border, ";",
      " padding:8px 12px; border-radius:0 4px 4px 0;",
      " margin:.6em 0; font-size:.88em; line-height:1.55;"
    ),
    HTML(paste0(icon_chr, " ", ...))
  )
}

#' Warning callout
help_warn <- function(...) {
  help_note(..., icon_chr = "⚠️", color = "#fff8e1", border = "#F9A825")
}

#' Success callout
help_ok <- function(...) {
  help_note(..., icon_chr = "✅", color = "#e8f5e9", border = "#388E3C")
}

#' Small two-column definition list
def_list <- function(...) {
  rows <- list(...)   # each element is c(term, def)
  items <- lapply(rows, function(r) {
    tagList(
      tags$dt(style = "font-weight:600; margin-top:.5em;", r[[1]]),
      tags$dd(style = "margin-left:1.4em; color:#444;",   HTML(r[[2]]))
    )
  })
  tags$dl(style = "font-size:.88em;", .list = items)
}

# ── UI function ───────────────────────────────────────────────────────────────

helpUI <- function(id) {
  ns <- NS(id)

  tagList(

    tags$div(
      style = "max-width:800px; padding-bottom:2em;",

      tags$h4(
        icon("question-circle", style = "margin-right:6px; color:#2E75B6;"),
        "User Guide",
        style = "color:#1B3A5C; margin-bottom:.2em;"
      ),
      tags$p(
        style = "color:#666; margin-bottom:1.2em; font-size:.92em;",
        "Quick answers for every part of the app.  ",
        "A full markdown reference is also available at ",
        tags$code("docs/USER_GUIDE.md"), " in the repository."
      ),

      # ── Accordion ──────────────────────────────────────────────────────────
      accordion(
        id = ns("help_accordion"),
        open = NULL,           # all panels collapsed on first load

        # ── 1. Quick Start ─────────────────────────────────────────────────
        accordion_panel(
          title = tagList(icon("rocket"), " Quick Start"),
          value = "qs",

          help_h("What you need"),
          tags$ol(
            tags$li(
              tags$strong("VCF file"), " — biallelic SNPs, any standard caller ",
              "(GATK, STACKS, FreeBayes, DeepVariant).  Must end in ",
              help_ic(".vcf"), " or ", help_ic(".vcf.gz"), "."
            ),
            tags$li(
              tags$strong("Metadata CSV"), " — one row per sample with columns: ",
              help_ic("sample_id, population, latitude, longitude"), "."
            )
          ),

          help_h("Three steps to your first map"),
          tags$ol(
            style = "line-height:1.8;",
            tags$li("Open the sidebar → upload your VCF and metadata CSV."),
            tags$li("The app switches automatically to the ", tags$strong("Map"), " tab."),
            tags$li("Use the metric selector above the map to explore He, Ho, π, F, Ne, and Offset overlays.")
          ),

          help_ok(
            "Want to explore without real data?  ",
            "Click <strong>Load demo dataset</strong> in the sidebar to load a ",
            "synthetic 4-population dataset covering all analyses."
          ),

          help_h("Minimum data requirements"),
          def_list(
            c("Samples",  "≥ 2 samples per population; ≥ 2 populations for FST / ADMIXTURE."),
            c("Loci",     "≥ 100 biallelic SNPs recommended; ≥ 500 for stable LFMM2 estimates."),
            c("VCF",      "Must be diploid.  Multi-allelic sites are silently filtered."),
            c("Metadata", paste0("sample_id values must match the VCF sample column exactly ",
                                 "(case-sensitive)."))
          )
        ),

        # ── 2. Input formats ────────────────────────────────────────────────
        accordion_panel(
          title = tagList(icon("file-upload"), " Input Formats"),
          value = "inputs",

          help_h("Genomic data (VCF)"),
          tags$p(style = "font-size:.9em;",
            "Standard VCF 4.1+ format.  The app keeps only biallelic SNPs (",
            help_ic("ALT"), " field has a single non-reference allele) and samples ",
            "present in both the VCF and the metadata CSV."
          ),
          help_note("Indels, structural variants, and monomorphic sites are silently dropped during import."),

          help_h("Sample metadata CSV"),
          tags$p(style = "font-size:.9em;",
            "Required columns (exact names, any order):"
          ),
          tags$table(
            class = "table table-sm table-bordered",
            style = "font-size:.85em; width:auto;",
            tags$thead(tags$tr(
              tags$th("Column"),  tags$th("Type"),    tags$th("Example")
            )),
            tags$tbody(
              tags$tr(tags$td(help_ic("sample_id")),  tags$td("character"), tags$td("Sample_001")),
              tags$tr(tags$td(help_ic("population")), tags$td("character"), tags$td("PopA")),
              tags$tr(tags$td(help_ic("latitude")),   tags$td("numeric"),   tags$td("-34.5")),
              tags$tr(tags$td(help_ic("longitude")),  tags$td("numeric"),   tags$td("150.2"))
            )
          ),

          help_h("Environmental data CSV (Phase 4)"),
          tags$p(style = "font-size:.9em;",
            "Required for GEA and genomic offset if you are not using the ",
            "WorldClim auto-fetch.  One row per population:"
          ),
          tags$ul(style = "font-size:.88em;",
            tags$li("First column must be ", help_ic("population"),
                    " matching names in your metadata."),
            tags$li("Remaining columns are numeric environmental variables ",
                    "(e.g. ", help_ic("BIO1"), ", ", help_ic("BIO12"), ")."),
            tags$li("No missing values allowed — every cell must be filled.")
          ),
          help_warn(
            "WorldClim BIO temperature variables (BIO1–BIO11) are stored in ",
            "<strong>°C × 10</strong>.  BIO1 = 250 means 25.0 °C, not 250 °C.  ",
            "This is the standard WorldClim convention and is handled correctly ",
            "by the auto-fetch; if you upload your own CSV, use the same convention."
          )
        ),

        # ── 3. Phase 1 – Diversity Metrics ──────────────────────────────────
        accordion_panel(
          title = tagList(icon("dna"), " Phase 1 — Diversity Metrics"),
          value = "p1",

          help_h("What is computed"),
          def_list(
            c("He — Expected heterozygosity",
              "Probability that two alleles drawn at random are different.  ",
              "Ranges 0 – 1; higher is more diverse."),
            c("Ho — Observed heterozygosity",
              "Proportion of heterozygous genotypes in the sample.  ",
              "Should be close to He in a randomly mating population."),
            c("π — Nucleotide diversity",
              "Average pairwise nucleotide differences per site.  ",
              "Similar scale to He for biallelic SNPs."),
            c("F — Inbreeding coefficient",
              "(He − Ho) / He.  Positive values indicate excess homozygosity.  ",
              "Negative values can arise from technical artefacts or outbreeding.")
          ),

          help_h("Conservation thresholds (He)"),
          tags$p(style = "font-size:.85em; color:#555;",
            "Thresholds follow Frankham et al. (2014).  ",
            tags$strong("These were calibrated on microsatellite data."),
            "  SNP-based He filtered at MAF ≥ 0.05 runs systematically lower.  ",
            "Treat them as indicative, not hard cut-offs."
          ),
          tags$table(
            class = "table table-sm",
            style = "font-size:.85em; width:auto;",
            tags$thead(tags$tr(
              tags$th("He"),         tags$th("Status"),
              tags$th("Guidance")
            )),
            tags$tbody(
              tags$tr(tags$td("< 0.10"),      tags$td("🔴 Critical"),
                      tags$td("Urgent action")),
              tags$tr(tags$td("0.10 – 0.15"), tags$td("🟠 High risk"),
                      tags$td("Intervention likely needed")),
              tags$tr(tags$td("0.15 – 0.20"), tags$td("🟡 Moderate risk"),
                      tags$td("Monitor; consider genetic rescue")),
              tags$tr(tags$td("0.20 – 0.30"), tags$td("🟢 Lower risk"),
                      tags$td("Routine monitoring")),
              tags$tr(tags$td("> 0.30"),       tags$td("🔵 Healthy"),
                      tags$td("High diversity; reference candidate"))
            )
          ),

          help_h("Reading the map"),
          tags$ul(style = "font-size:.88em;",
            tags$li("Circle size scales with sample size (number of samples in that population)."),
            tags$li("Circle colour encodes the selected metric — see the legend on the right."),
            tags$li("Click any circle for a detailed tooltip with all Phase 1 statistics."),
            tags$li("Use the metric selector buttons above the map to switch between He, Ho, π, F, Ne, and Offset.")
          )
        ),

        # ── 4. Phase 2 – Population Structure ───────────────────────────────
        accordion_panel(
          title = tagList(icon("project-diagram"), " Phase 2 — Population Structure"),
          value = "p2",

          help_h("FST heatmap (Structure tab)"),
          tags$p(style = "font-size:.88em;",
            "Pairwise FST is computed with the Hudson et al. (1992) ratio-of-averages estimator.  ",
            "Values range 0–1: 0 = identical allele frequencies, 1 = fixed for ",
            "different alleles.  The heatmap uses a blue–red diverging colour scale ",
            "centred on the dataset median.  Populations are reordered by hierarchical ",
            "clustering (Ward linkage) to bring similar populations together."
          ),
          def_list(
            c("FST < 0.05",  "Negligible differentiation — gene flow likely."),
            c("0.05 – 0.15", "Moderate differentiation."),
            c("0.15 – 0.25", "Large differentiation — limited gene flow."),
            c("> 0.25",      "Very large differentiation — genetic isolation.")
          ),

          help_h("PCA (PCA tab)"),
          tags$p(style = "font-size:.88em;",
            "Principal Component Analysis on the sample genotype matrix (mean imputation ",
            "for missing data).  Each dot is one sample, coloured by population.  ",
            "PC1 and PC2 explain the most variance.  Populations that cluster tightly ",
            "are genetically similar; spread within a cluster reflects within-population ",
            "diversity."
          ),

          help_h("ADMIXTURE pies (Map tab)"),
          tags$p(style = "font-size:.88em;",
            "Upload a STRUCTURE-format Q-matrix from any external ADMIXTURE / STRUCTURE run ",
            "(first column = sample IDs, remaining columns = K ancestry proportions, ",
            "no header, space-delimited).  The app overlays pie charts on the map showing ",
            "the mean ancestry proportions for each population."
          ),
          help_note("ADMIXTURE computation is not performed in-browser — the Q-matrix must ",
                    "be generated externally and uploaded here.")
        ),

        # ── 5. Phase 3 – Inbreeding & Ne ────────────────────────────────────
        accordion_panel(
          title = tagList(icon("users"), " Phase 3 — F-ROH and Effective Population Size"),
          value = "p3",

          help_h("F-ROH"),
          tags$p(style = "font-size:.88em;",
            "Runs of Homozygosity inbreeding coefficient: fraction of the genome ",
            "covered by consecutive stretches of homozygous SNPs above a minimum ",
            "length threshold.  Captures ", tags$em("recent"), " inbreeding ",
            "(approximately the last 10 generations)."
          ),
          def_list(
            c("F-ROH ≥ 0.25", "Consistent with full-sibling or parent-offspring mating."),
            c("F-ROH ≥ 0.125","Consistent with half-sibling mating."),
            c("F-ROH ≥ 0.0625","Consistent with first-cousin mating.")
          ),
          help_note("Accuracy scales with marker density.  At least ~5,000 SNPs distributed across the genome recommended for reliable ROH detection; sparse or clustered markers inflate ROH length estimates."),

          help_h("Effective Population Size (Ne)"),
          tags$p(style = "font-size:.88em;",
            "Estimated using the single-sample LD-based method of ",
            "Waples & Do (2008).  Ne reflects the size of an idealised random-mating ",
            "population that would drift at the same rate as the sample."
          ),
          def_list(
            c("Ne < 50",  "Immediate extinction risk (IUCN criterion)."),
            c("Ne < 100", "High risk; limited ability to respond to selection."),
            c("Ne < 500", "Long-term viability concern under genetic drift.")
          ),
          help_warn("LD-based Ne estimates assume the sample is from a single ",
                    "randomly mating unit.  Admixed or spatially structured samples ",
                    "produce downward-biased (artificially small) Ne estimates.")
        ),

        # ── 6. Phase 4 – Adaptive Potential ─────────────────────────────────
        accordion_panel(
          title = tagList(icon("leaf"), " Phase 4 — GEA and Genomic Offset"),
          value = "p4",

          help_h("Environmental data"),
          tags$p(style = "font-size:.88em;",
            "Phase 4 analyses require per-population environmental data.  Use one of:"
          ),
          tags$ul(style = "font-size:.88em;",
            tags$li(tags$strong("Auto-fetch (WorldClim):"),
                    " extracts current climate for each population location from ",
                    "WorldClim 2.1 rasters (requires internet).  Select BIO variables in the sidebar."),
            tags$li(tags$strong("Auto-fetch (CMIP6 future):"),
                    " downloads projected climate from the ACCESS-CM2 model under your chosen SSP ",
                    "scenario (2041–2060) for use in the offset calculation."),
            tags$li(tags$strong("Upload CSV:"), " provide your own environmental table — see Input Formats above.")
          ),

          help_h("RDA (Redundancy Analysis)"),
          tags$p(style = "font-size:.88em;",
            "Multivariate regression of population allele frequencies on environmental ",
            "variables.  Geography (latitude + longitude) is partialled out so only ",
            "the environmental signal remains.  The biplot shows populations (dots) and ",
            "environmental vectors; the angle between a vector and a population indicates ",
            "the sign of association."
          ),
          help_note("With a single environmental variable the constrained ordination space is 1D — populations are plotted on one axis only.  Select ≥ 2 variables for a 2D biplot."),

          help_h("LFMM2 (Latent Factor Mixed Model)"),
          tags$p(style = "font-size:.88em;",
            "Per-locus test for genotype-environment association.  Latent factors ",
            "control for population structure; the Genomic Inflation Factor (GIF) ",
            "calibrates p-values.  Multiple testing is corrected with ",
            "Benjamini-Hochberg FDR.  The Manhattan plot shows –log10(adjusted p-value); ",
            "the dashed red line marks FDR = 0.05."
          ),
          def_list(
            c("GIF ≈ 1",   "Well-calibrated; p-values reliable."),
            c("GIF 1.5–2", "Mild inflation — check for batch effects or unmodelled structure."),
            c("GIF > 2",   "Strong undercorrection; increase K or add geographic covariates."),
            c("GIF < 0.8", "Over-correction; consider reducing K.")
          ),

          help_h("Genomic Offset"),
          tags$p(style = "font-size:.88em;",
            "RDA-based climate vulnerability (Fitzpatrick & Keller 2015).  ",
            "For each population, the Euclidean distance is computed between its ",
            "current environmental position and its projected future position in the ",
            "constrained ordination space defined by the RDA.  Larger offset = larger ",
            "adaptive gap under climate change."
          ),
          help_warn(
            "Genomic offset is a relative, not absolute, measure.  ",
            "A high-offset population is not necessarily doomed — it may have high ",
            "standing genetic variation, phenotypic plasticity, or access to migration ",
            "corridors.  Always interpret alongside Ne, connectivity, and local demography."
          ),
          def_list(
            c("Offset (normalised) ≥ 0.75", "Severe climate risk."),
            c("0.50 – 0.75",                 "High climate risk."),
            c("0.25 – 0.50",                 "Moderate climate risk."),
            c("< 0.25",                       "Low climate risk.")
          )
        ),

        # ── 7. Export ────────────────────────────────────────────────────────
        accordion_panel(
          title = tagList(icon("file-export"), " Phase 5 — Export"),
          value = "p5",

          help_h("Conservation Report (HTML / PDF)"),
          tags$p(style = "font-size:.88em;",
            "A one-click report summarising all available analyses.  The HTML version ",
            "can be opened in any browser.  The PDF version requires the ",
            help_ic("pagedown"), " R package and Google Chrome (or Chromium) to be installed ",
            "on the server; if unavailable the HTML version is produced instead with a ",
            "notification."
          ),
          help_note("The report only includes sections for analyses you have run — for example, ",
                    "the GEA and Offset sections appear only if environmental data was uploaded."),

          help_h("Statistics CSV"),
          tags$p(style = "font-size:.88em;",
            "Per-population table with He, Ho, π, F, n_samples, n_loci, FST (mean), ",
            "Ne, and offset (if available).  Compatible with Excel, R, and Python."
          ),

          help_h("Excel workbook (.xlsx)"),
          tags$p(style = "font-size:.88em;",
            "Multi-sheet workbook containing: Diversity Statistics, Pairwise FST, ",
            "LFMM2 Results, and Genomic Offset (only sheets with data are included)."
          ),

          help_h("Map PNG"),
          tags$p(style = "font-size:.88em;",
            "Click the ", tags$strong("camera icon"), " inside the Map tab to save the ",
            "current map viewport as a PNG.  This uses the Leaflet screenshot plugin — ",
            "no data leaves the browser."
          ),

          help_h("Chart PNG"),
          tags$p(style = "font-size:.88em;",
            "Hover over any Plotly chart (PCA, ADMIXTURE bar, F-ROH, Ne, ",
            "GEA Manhattan, Offset bar) and click the ",
            tags$strong("camera icon"), " in the toolbar to download at 2× (print-quality) resolution.  ",
            "Note: the FST heatmap is an HTML table, not a Plotly chart — use your browser's ",
            "Print or screenshot tool to capture it."
          ),

          help_note(
            "All downloads are generated client-side or in the R session — no data is ",
            "sent to any external server."
          )
        ),

        # ── 8. Troubleshooting ───────────────────────────────────────────────
        accordion_panel(
          title = tagList(icon("tools"), " Troubleshooting"),
          value = "trouble",

          help_h("Upload errors"),
          tags$table(
            class = "table table-sm table-bordered",
            style = "font-size:.83em;",
            tags$thead(tags$tr(tags$th("Error"), tags$th("Likely cause"), tags$th("Fix"))),
            tags$tbody(
              tags$tr(
                tags$td("No samples in common"),
                tags$td("sample_id in CSV does not match VCF sample names"),
                tags$td("Check for trailing spaces, case differences, or underscores vs hyphens")
              ),
              tags$tr(
                tags$td("No biallelic SNPs retained"),
                tags$td("VCF contains only multi-allelic or monomorphic sites"),
                tags$td(HTML("Filter with <code>bcftools view -m2 -M2 -v snps</code> before upload"))
              ),
              tags$tr(
                tags$td("Metadata missing required columns"),
                tags$td("Column names misspelled or in wrong case"),
                tags$td(HTML("Columns must be exactly: <code>sample_id, population, latitude, longitude</code>"))
              ),
              tags$tr(
                tags$td("File too large"),
                tags$td("VCF exceeds the upload limit"),
                tags$td(HTML("Thin to ~50,000 SNPs with <code>vcftools --thin 1000</code> or LD pruning in PLINK"))
              )
            )
          ),

          help_h("Analysis errors"),
          tags$table(
            class = "table table-sm table-bordered",
            style = "font-size:.83em;",
            tags$thead(tags$tr(tags$th("Message"), tags$th("Cause"), tags$th("Fix"))),
            tags$tbody(
              tags$tr(
                tags$td("RDA failed: must have ≥ 2 populations"),
                tags$td("Only one population in dataset"),
                tags$td("RDA and FST require ≥ 2 populations")
              ),
              tags$tr(
                tags$td("Ne estimate Inf or negative"),
                tags$td("Insufficient LD signal or very small sample"),
                tags$td("Increase sample size to ≥ 5 per population; increase SNP count")
              ),
              tags$tr(
                tags$td("LFMM2: all p-values = 1"),
                tags$td("Too few loci or K > number of populations − 1"),
                tags$td("Reduce K; ensure ≥ 500 SNPs")
              ),
              tags$tr(
                tags$td("WorldClim fetch failed"),
                tags$td("No internet connection or coordinates out of range"),
                tags$td("Check coordinates; upload environmental CSV manually instead")
              )
            )
          ),

          help_h("Report not generating"),
          tags$ul(style = "font-size:.88em;",
            tags$li("The HTML report requires the ", help_ic("rmarkdown"), " and ", help_ic("knitr"),
                    " packages — run ", help_ic("source('install_packages.R')"), " to install."),
            tags$li("The PDF report additionally requires ", help_ic("pagedown"),
                    " and Google Chrome.  If Chrome is not found, the app falls back to HTML."),
            tags$li("The XLSX export requires the ", help_ic("openxlsx"), " package.")
          )
        ),

        # ── 9. Glossary ──────────────────────────────────────────────────────
        accordion_panel(
          title = tagList(icon("book"), " Glossary"),
          value = "glossary",

          def_list(
            c("Allele",               "One of two or more alternative forms of a gene at a given locus."),
            c("Biallelic SNP",        "A single-nucleotide variant with exactly two alleles in the dataset."),
            c("BIO variable",         paste0(
                "WorldClim bioclimatic variable.  BIO1 = Annual Mean Temperature, ",
                "BIO12 = Annual Precipitation.  Temperature variables are in °C × 10.")),
            c("CMIP6",                paste0(
                "Coupled Model Intercomparison Project Phase 6.  Provides standardised ",
                "future climate projections.  Access-CM2 (Australia) is the default model used here.")),
            c("FST",                  paste0(
                "Fixation index: proportion of total genetic variance due to allele frequency ",
                "differences among populations.  Hudson et al. (1992) ratio-of-averages estimator used here.")),
            c("GIF (Genomic Inflation Factor)",
                                      paste0(
                "Ratio of observed median chi-squared to expected under the null.  ",
                "GIF ≈ 1 indicates well-calibrated tests; GIF > 1 suggests inflation (confounding).")),
            c("Hardy-Weinberg Equilibrium (HWE)",
                                      paste0(
                "Expected genotype frequencies in a large, randomly mating population with no ",
                "selection, mutation, or migration.  Deviations indicate structure or inbreeding.")),
            c("He",                   "Expected heterozygosity under Hardy-Weinberg equilibrium."),
            c("Ho",                   "Observed heterozygosity from the data."),
            c("K (ADMIXTURE)",        "Number of ancestral clusters assumed in an ADMIXTURE analysis."),
            c("K (LFMM2)",            paste0(
                "Number of latent factors used to control for population structure in LFMM2.  ",
                "Typically set to the number of PCs that explain the main structure.")),
            c("LD (Linkage Disequilibrium)",
                                      "Non-random association between alleles at different loci."),
            c("LFMM2",                paste0(
                "Latent Factor Mixed Model (Caye et al. 2019).  A fast frequentist method for ",
                "genotype-environment association testing.")),
            c("MAF",                  "Minor Allele Frequency — frequency of the less common allele."),
            c("Ne",                   paste0(
                "Effective population size.  The size of an idealised randomly mating population ",
                "that drifts at the same rate as the actual population.")),
            c("π (nucleotide diversity)",
                                      "Average proportion of nucleotide sites that differ between two randomly chosen sequences."),
            c("RDA",                  paste0(
                "Redundancy Analysis — a constrained ordination that explains genetic variation ",
                "using environmental predictors while partialling out spatial structure.")),
            c("ROH (Runs of Homozygosity)",
                                      paste0(
                "Long stretches of consecutive homozygous SNPs arising from identity-by-descent ",
                "in recent ancestors.")),
            c("SNP",                  "Single Nucleotide Polymorphism — a single-base difference between individuals."),
            c("SSP",                  paste0(
                "Shared Socioeconomic Pathway.  SSP2-4.5 = intermediate pathway ('middle of the road'); ",
                "SSP5-8.5 = high-emission 'fossil-fuelled development' scenario.")),
            c("VCF",                  paste0(
                "Variant Call Format.  Standard text file describing genomic variant calls, ",
                "produced by GATK, FreeBayes, STACKS, DeepVariant, and other variant callers."))
          )
        )   # end glossary accordion_panel

      )   # end accordion
    )   # end max-width div
  )   # end tagList
}

# ── Server function (no logic needed — pure UI module) ────────────────────────
helpServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # No server-side logic required for a static help tab.
    invisible(NULL)
  })
}
