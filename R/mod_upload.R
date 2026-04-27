# =============================================================================
# mod_upload.R — File upload Shiny module
# Handles VCF + metadata CSV upload, or loads built-in demo data.
# Returns a reactive with parsed genotype matrix, metadata, and validation info.
# =============================================================================

# ── UI ─────────────────────────────────────────────────────────────────────────

#' Upload module UI
#'
#' @param id Module namespace ID
uploadUI <- function(id) {
  ns <- NS(id)

  tagList(
    # ── Demo banner ───────────────────────────────────────────────────────────
    div(
      class = "demo-banner",
      tags$p(
        tags$strong("New here?"),
        " Load the built-in demo dataset to explore the app before uploading your own data."
      ),
      actionButton(
        ns("load_demo"),
        label = tagList(icon("play-circle"), " Load demo dataset"),
        class = "btn btn-success btn-sm w-100"
      )
    ),

    tags$hr(style = "margin: 12px 0;"),

    # ── VCF upload ────────────────────────────────────────────────────────────
    tags$label("VCF file", class = "upload-label"),
    tags$p(
      class = "upload-hint",
      HTML(
        "Biallelic SNPs, any variant caller (GATK, STACKS, FreeBayes).<br>",
        "Plain <code>.vcf</code> or gzipped <code>.vcf.gz</code>. Max 200 MB."
      )
    ),
    fileInput(
      ns("vcf_file"),
      label       = NULL,
      accept      = c(".vcf", ".gz"),
      buttonLabel = "Browse…",
      placeholder = "No file selected"
    ),

    # ── Metadata CSV upload ───────────────────────────────────────────────────
    tags$label("Sample metadata CSV", class = "upload-label"),
    tags$p(
      class = "upload-hint",
      HTML(
        "Required columns: <code>sample_id</code>, <code>population</code>,",
        "<code>latitude</code>, <code>longitude</code>"
      )
    ),
    fileInput(
      ns("meta_file"),
      label       = NULL,
      accept      = ".csv",
      buttonLabel = "Browse…",
      placeholder = "No file selected"
    ),

    # ── Filter options (collapsible) ──────────────────────────────────────────
    tags$details(
      style = "margin-bottom: 8px;",
      tags$summary(
        style = "cursor:pointer; font-size:0.82em; color:#555;",
        "Advanced filter options"
      ),
      div(
        style = "padding: 8px 0 4px;",
        sliderInput(
          ns("maf_threshold"),
          label = HTML("Minor allele frequency (MAF) threshold"),
          min = 0, max = 0.20, value = 0.05, step = 0.01,
          ticks = FALSE
        ),
        tags$p(
          class = "upload-hint",
          style = "margin-top:-8px;",
          HTML("Lower = more loci retained, more rare variants included.<br/>
               Raise if you want only common, high-confidence SNPs.")
        ),
        sliderInput(
          ns("miss_threshold"),
          label = HTML("Max missing data per site"),
          min = 0, max = 0.50, value = 0.20, step = 0.05,
          ticks = FALSE
        ),
        tags$p(
          class = "upload-hint",
          style = "margin-top:-8px;",
          HTML("Higher = more loci retained despite gaps in sequencing coverage.<br/>
               Lower = stricter; use 0.10 for high-quality WGS data.")
        )
      )
    ),

    # ── Analyse + Reset buttons ───────────────────────────────────────────────
    div(
      class = "d-flex gap-2 mt-1",
      actionButton(
        ns("run_analysis"),
        label = tagList(icon("dna"), " Run analysis"),
        class = "btn btn-primary flex-grow-1"
      ),
      actionButton(
        ns("reset_data"),
        label = icon("times"),
        class = "btn btn-outline-secondary",
        title = "Clear loaded data and reset"
      )
    ),

    # ── Status output (progress + validation report) ──────────────────────────
    uiOutput(ns("upload_status")),

    # ── Active dataset indicator ──────────────────────────────────────────────
    uiOutput(ns("active_dataset")),

    tags$hr(style = "margin: 16px 0 8px;"),

    # ── Metadata CSV template download ────────────────────────────────────────
    tags$p(class = "upload-hint", "Need a template?"),
    downloadButton(
      ns("dl_template"),
      label = tagList(icon("download"), " Metadata CSV template"),
      class = "btn btn-outline-secondary btn-sm w-100"
    )
  )
}

# ── Server ──────────────────────────────────────────────────────────────────────

#' Upload module server
#'
#' @param id Module namespace ID
#' @return Reactive list:
#'   $gt           — genotype matrix (samples × loci), or NULL
#'   $metadata     — cleaned data frame, or NULL
#'   $scenario_info — descriptive text string
#'   $is_demo      — logical
#'   $filter_log   — named integer vector (VCF filter step counts), or NULL
#'   $validation   — list($errors, $warnings, $info) from validate_gt_meta(), or NULL
#'   $ready        — logical: TRUE when data is loaded and valid
uploadServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    MAX_VCF_BYTES <- 200 * 1024 * 1024   # 200 MB hard limit

    # Reactive state
    rv <- reactiveValues(
      gt            = NULL,
      metadata      = NULL,
      demo_qmatrix  = NULL,   # demo Q-matrix (K=2) for ADMIXTURE tab, or NULL
      scenario_info = NULL,
      is_demo       = FALSE,
      filter_log    = NULL,
      validation    = NULL,
      steps         = NULL,   # character vector of progress steps completed
      status_type   = NULL,   # "success" | "error" | "warning" | "info"
      error_msg     = NULL
    )

    # ── Helper: push a progress step message ────────────────────────────────
    push_step <- function(msg) {
      rv$steps <- c(rv$steps, msg)
    }

    # ── Load demo data ───────────────────────────────────────────────────────
    observeEvent(input$load_demo, {
      shinyjs::disable("load_demo")
      updateActionButton(session, "load_demo",
                         label = tagList(
                           tags$span(class = "btn-loading", icon("play-circle"), " Loading…")
                         ))
      on.exit({
        shinyjs::enable("load_demo")
        updateActionButton(session, "load_demo",
                           label = tagList(icon("play-circle"), " Load demo dataset"))
      }, add = TRUE)

      rv$gt <- rv$metadata <- rv$filter_log <- rv$validation <- rv$error_msg <- NULL
      rv$steps      <- character(0)
      rv$status_type <- "info"

      push_step("\u23f3 Generating demo dataset…")

      tryCatch({
        demo         <- generate_demo_data()
        rv$gt            <- demo$gt
        rv$metadata      <- demo$metadata
        rv$demo_qmatrix  <- demo$qmatrix
        rv$scenario_info <- demo$scenario_info
        rv$is_demo       <- TRUE
        rv$filter_log    <- NULL
        rv$validation    <- list(
          errors   = character(0),
          warnings = character(0),
          info     = c(
            paste0(nrow(demo$metadata), " samples across ",
                   length(unique(demo$metadata$population)), " populations."),
            paste0(ncol(demo$gt), " polymorphic SNPs (simulated).")
          ),
          matched_ids = demo$metadata$sample_id,
          valid = TRUE
        )
        rv$status_type <- "success"
        rv$steps <- c(
          paste0("\u2713 Demo data loaded: ", nrow(demo$metadata), " samples, ",
                 length(unique(demo$metadata$population)), " populations, ",
                 ncol(demo$gt), " SNPs.")
        )
      }, error = function(e) {
        rv$status_type <- "error"
        rv$error_msg   <- paste("Error generating demo data:", e$message)
      })
    })

    # ── Run analysis on uploaded files ────────────────────────────────────────
    observeEvent(input$run_analysis, {

      # Pre-flight: give an explicit, visible error if either file is missing.
      # req() would silently cancel the observer with no user feedback.
      missing_files <- c(
        if (is.null(input$vcf_file))  "VCF file (.vcf / .vcf.gz)",
        if (is.null(input$meta_file)) "metadata CSV"
      )
      if (length(missing_files) > 0) {
        rv$steps       <- character(0)
        rv$status_type <- "error"
        rv$error_msg   <- paste0(
          "Missing required file", if (length(missing_files) > 1) "s" else "", ": ",
          paste(missing_files, collapse = " and "), ". ",
          "Please upload both files before running the analysis."
        )
        return()
      }

      # Both files present — proceed with analysis
      # Disable button and show loading state while processing
      updateActionButton(session, "run_analysis",
                         label = tagList(
                           tags$span(class = "btn-loading", icon("dna"), " Analysing…")
                         ))
      shinyjs::disable("run_analysis")
      on.exit({
        shinyjs::enable("run_analysis")
        updateActionButton(session, "run_analysis",
                           label = tagList(icon("dna"), " Run analysis"))
      }, add = TRUE)

      # Reset state
      rv$gt <- rv$metadata <- rv$filter_log <- rv$validation <- rv$error_msg <- NULL
      rv$steps       <- character(0)
      rv$is_demo     <- FALSE
      rv$status_type <- "info"

      # ── Step 1: File size guard ─────────────────────────────────────────────
      push_step("\u23f3 Checking file sizes…")
      vcf_size <- file.info(input$vcf_file$datapath)$size
      if (!is.na(vcf_size) && vcf_size > MAX_VCF_BYTES) {
        rv$status_type <- "error"
        rv$error_msg   <- paste0(
          "VCF file is too large (",
          round(vcf_size / 1024 / 1024, 1), " MB). ",
          "Maximum allowed size is 200 MB. ",
          "Please filter your VCF to a representative set of loci ",
          "(e.g. thin to one SNP per 10 kb with vcftools --thin)."
        )
        return()
      }
      push_step(paste0(
        "\u2713 VCF: ", round(vcf_size / 1024 / 1024, 2), " MB | ",
        "CSV: ", round(file.info(input$meta_file$datapath)$size / 1024, 1), " KB"
      ))

      # ── Step 2: Parse metadata CSV ──────────────────────────────────────────
      push_step("\u23f3 Parsing metadata CSV…")
      meta <- tryCatch(
        withCallingHandlers(
          parse_metadata_csv(input$meta_file$datapath),
          warning = function(w) {
            push_step(paste0("\u26a0\ufe0f ", conditionMessage(w)))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          rv$status_type <- "error"
          rv$error_msg   <- paste0("Metadata CSV error: ", e$message)
          NULL
        }
      )
      if (is.null(meta)) return()
      push_step(paste0("\u2713 Metadata: ", nrow(meta), " samples, ",
                       length(unique(meta$population)), " populations."))

      # ── Step 3: Parse VCF ───────────────────────────────────────────────────
      push_step("\u23f3 Parsing VCF (this may take a moment for large files)…")
      vcf_data <- tryCatch(
        parse_vcf(input$vcf_file$datapath,
                  maf         = input$maf_threshold,
                  max_missing = input$miss_threshold),
        error = function(e) {
          rv$status_type <- "error"
          rv$error_msg   <- paste0("VCF parsing error: ", e$message)
          NULL
        }
      )
      if (is.null(vcf_data)) return()

      fl <- vcf_data$filter_log
      push_step(paste0(
        "\u2713 VCF parsed. Filter log: ",
        fl["raw"], " raw sites \u2192 ",
        fl["biallelic_snps"], " biallelic SNPs \u2192 ",
        fl["after_missing"], " after missingness filter \u2192 ",
        fl["after_maf"], " after MAF filter \u2192 ",
        fl["final_polymorphic"], " final polymorphic SNPs retained."
      ))

      # Warn if very few SNPs retained
      if (fl["final_polymorphic"] < 50) {
        push_step(paste0(
          "\u26a0\ufe0f Only ", fl["final_polymorphic"],
          " SNPs retained. Diversity estimates may be unreliable. ",
          "Consider relaxing the MAF or missing-data filter."
        ))
      }

      # ── Step 4: Validate compatibility ─────────────────────────────────────
      push_step("\u23f3 Validating sample ID overlap…")
      val <- validate_gt_meta(vcf_data$gt, meta)

      # Surface all errors, warnings, and info from the structured report
      for (e in val$errors)   push_step(paste0("\u274c ", e))
      for (w in val$warnings) push_step(paste0("\u26a0\ufe0f ", w))
      for (i in val$info)     push_step(paste0("\u2139\ufe0f ", i))

      if (!val$valid) {
        rv$status_type <- "error"
        rv$error_msg   <- paste(val$errors, collapse = "\n")
        return()
      }

      # ── Commit results ─────────────────────────────────────────────────────
      rv$gt            <- vcf_data$gt
      rv$metadata      <- meta
      rv$filter_log    <- vcf_data$filter_log
      rv$validation    <- val

      # Use matched counts (val$matched_ids) so the banner reflects what
      # was actually analysed, not the raw VCF / CSV totals.
      matched_meta   <- meta[meta$sample_id %in% val$matched_ids, ]
      matched_n_pops <- length(unique(matched_meta$population))
      rv$scenario_info <- paste0(
        input$vcf_file$name, " — ",
        length(val$matched_ids), " samples, ",
        vcf_data$n_loci, " polymorphic SNPs, ",
        matched_n_pops, " population", if (matched_n_pops != 1) "s" else ""
      )
      rv$status_type <- if (length(val$warnings) > 0) "warning" else "success"

      push_step(paste0(
        "\u2705 Analysis ready: ", vcf_data$n_loci, " SNPs \u00d7 ",
        length(unique(meta$population)), " populations."
      ))
    })

    # ── Status display ────────────────────────────────────────────────────────
    output$upload_status <- renderUI({
      # Show error panel if there is a hard error
      if (!is.null(rv$error_msg)) {
        return(
          div(
            class = "alert alert-danger",
            style = "margin-top:10px; font-size:0.84em; padding:10px 14px;",
            tags$strong("\u274c Error"), tags$br(),
            rv$error_msg
          )
        )
      }

      # Build step-by-step progress log
      if (is.null(rv$steps) || length(rv$steps) == 0) return(NULL)

      step_html <- lapply(rv$steps, function(s) {
        colour <- if (startsWith(s, "\u2705") || startsWith(s, "\u2713")) "#155724"
                  else if (startsWith(s, "\u26a0")) "#856404"
                  else if (startsWith(s, "\u2139")) "#0c5460"
                  else "#333"
        tags$p(style = paste0("margin:2px 0; color:", colour, ";"), s)
      })

      alert_class <- switch(
        rv$status_type %||% "info",
        "success" = "alert alert-success",
        "warning" = "alert alert-warning",
        "error"   = "alert alert-danger",
        "alert alert-info"
      )

      div(
        class = alert_class,
        style = "margin-top:10px; font-size:0.83em; padding:10px 14px;",
        step_html
      )
    })

    # ── Reset / clear data ────────────────────────────────────────────────────
    observeEvent(input$reset_data, {
      rv$gt <- rv$metadata <- rv$filter_log <- rv$validation <-
        rv$error_msg <- rv$scenario_info <- rv$demo_qmatrix <- NULL
      rv$steps       <- character(0)
      rv$is_demo     <- FALSE
      rv$status_type <- NULL
    })

    # ── Active dataset indicator ──────────────────────────────────────────────
    output$active_dataset <- renderUI({
      req(rv$scenario_info)
      div(
        style = paste0(
          "font-size:0.78em; color:#555; background:#f8f9fa;",
          " border:1px solid #dee2e6; border-radius:4px;",
          " padding:5px 8px; margin-top:8px; line-height:1.4;"
        ),
        icon("circle", style = paste0(
          "color:", if (isTRUE(rv$is_demo)) "#F57C00" else "#388E3C",
          "; font-size:8px; margin-right:4px; vertical-align:middle;"
        )),
        rv$scenario_info
      )
    })

    # ── Metadata CSV template download ────────────────────────────────────────
    output$dl_template <- downloadHandler(
      filename = "popgen_map_metadata_template.csv",
      content  = function(file) {
        template <- data.frame(
          sample_id  = paste0("sample_", 1:5),
          population = c("Pop_A", "Pop_A", "Pop_B", "Pop_B", "Pop_C"),
          latitude   = c(1.500, 1.510, 0.800, 0.815, 1.230),
          longitude  = c(114.500, 114.520, 113.900, 113.920, 116.050),
          stringsAsFactors = FALSE
        )
        write.csv(template, file, row.names = FALSE)
      }
    )

    # ── Return reactive data ──────────────────────────────────────────────────
    return(reactive(list(
      gt            = rv$gt,
      metadata      = rv$metadata,
      demo_qmatrix  = rv$demo_qmatrix,   # non-NULL only in demo mode
      scenario_info = rv$scenario_info,
      is_demo       = rv$is_demo,
      filter_log    = rv$filter_log,
      validation    = rv$validation,
      ready         = !is.null(rv$gt) && !is.null(rv$metadata)
    )))
  })
}
