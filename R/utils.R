# =============================================================================
# utils.R — Shared utility functions
# popgen-map Phase 1
# =============================================================================

#' Null-coalescing operator
#'
#' Returns `a` if it is non-NULL, non-empty, and non-NA; otherwise returns `b`.
#' Defined here (sourced early) so all modules can use it without cross-sourcing.
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

#' Parse and validate the sample metadata CSV
#'
#' @param path File path or dataframe
#' @return Cleaned data frame with columns: sample_id, population, latitude, longitude
parse_metadata_csv <- function(path) {
  if (is.data.frame(path)) {
    df <- path
  } else {
    df <- tryCatch(
      read.csv(path, stringsAsFactors = FALSE, strip.white = TRUE),
      error = function(e) stop(paste("Could not read metadata CSV:", e$message))
    )
  }

  # Normalise column names (lowercase, trim whitespace)
  names(df) <- trimws(tolower(names(df)))

  # Accept common alternative column names
  rename_map <- list(
    sample_id  = c("sample_id", "sampleid", "id", "individual", "ind_id", "ind"),
    population = c("population", "pop", "pop_id", "site", "location"),
    latitude   = c("latitude", "lat", "y"),
    longitude  = c("longitude", "lon", "long", "x")
  )

  for (canonical in names(rename_map)) {
    candidates <- rename_map[[canonical]]
    match      <- intersect(candidates, names(df))[1]
    if (!is.na(match) && match != canonical) {
      names(df)[names(df) == match] <- canonical
    }
  }

  required <- c("sample_id", "population", "latitude", "longitude")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(paste0(
      "Metadata CSV is missing required columns: ", paste(missing, collapse = ", "),
      "\nRequired columns (case-insensitive): sample_id, population, latitude, longitude"
    ))
  }

  # Coerce coordinate columns
  df$latitude  <- suppressWarnings(as.numeric(df$latitude))
  df$longitude <- suppressWarnings(as.numeric(df$longitude))

  bad_coords <- which(is.na(df$latitude) | is.na(df$longitude))
  if (length(bad_coords) > 0) {
    warning(sprintf(
      "%d rows have missing/invalid coordinates and will be excluded: rows %s",
      length(bad_coords), paste(bad_coords, collapse = ", ")
    ))
    df <- df[!is.na(df$latitude) & !is.na(df$longitude), ]
  }

  if (nrow(df) == 0) stop("No valid rows remain after coordinate validation.")

  df[, c("sample_id", "population", "latitude", "longitude")]
}

#' Format a numeric value for display — rounds and adds NA handling
fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "—", formatC(round(x, digits), format = "f", digits = digits))
}

#' Build a plain-language HTML popup for a population
#'
#' @param row A single-row data frame from the pop_stats output.
#'            May optionally contain a column "Ne" (numeric effective pop size)
#'            added by the map module when the Ne metric is active.
#' @return HTML string
build_popup_html <- function(row) {
  risk       <- he_risk_label(row$He)
  risk_col   <- he_colour(row$He)

  # Inbreeding alert — conservation-actionable text at F > 0.1
  inb_note   <- if (!is.na(row$inbreeding_f) && row$inbreeding_f > 0.1) {
    paste0(
      '<p style="color:#B71C1C; font-size:0.85em; margin-top:6px;">',
      "⚠️ Inbreeding F = ", fmt(row$inbreeding_f, 2),
      " — elevated inbreeding detected; ",
      "consider genetic rescue or a managed breeding programme.</p>"
    )
  } else ""

  # Ne row + risk badge: only shown when the Ne column is present and non-NA.
  # The map module augments the df with a "Ne" column before rebuilding popups.
  ne_val      <- if ("Ne" %in% names(row) && !is.na(row$Ne)) {
    format(round(row$Ne), big.mark = ",")
  } else NULL
  ne_row      <- if (!is.null(ne_val)) popup_row("Eff. pop. size (Ne)", ne_val) else ""

  # Second risk badge for Ne (shown only when Ne is available).
  # ne_val is non-NULL only when "Ne" %in% names(row) && !is.na(row$Ne),
  # so a single NULL-check is sufficient.
  ne_badge <- if (!is.null(ne_val)) {
    ne_risk <- ne_risk_label(row$Ne)
    ne_col  <- ne_colour(row$Ne)
    paste0(
      "<span style='display:inline-block; padding:2px 8px; border-radius:12px; background:",
      ne_col, "; color:white; font-size:0.8em; margin-bottom:8px; margin-left:4px;'>",
      "Ne: ", ne_risk, "</span>"
    )
  } else ""

  paste0(
    "<div style='font-family:Arial,sans-serif; min-width:200px;'>",
    "<h4 style='margin:0 0 6px; color:#1B3A5C;'>", htmltools::htmlEscape(row$population), "</h4>",
    "<span style='display:inline-block; padding:2px 8px; border-radius:12px; background:",
    risk_col, "; color:white; font-size:0.8em; margin-bottom:8px;'>", risk, "</span>",
    ne_badge,
    "<table style='width:100%; border-collapse:collapse; font-size:0.88em;'>",
    popup_row("Samples (n)",              row$n),
    popup_row("SNPs analysed",            row$n_loci),
    popup_row("Expected He",              fmt(row$He)),
    popup_row("Observed Ho",              fmt(row$Ho)),
    popup_row("Nucleotide div. (π)", fmt(row$pi)),
    popup_row("Inbreeding F",             fmt(row$inbreeding_f, 3)),
    ne_row,
    "</table>",
    inb_note,
    "<p style='font-size:0.75em; color:#777; margin-top:6px;'>",
    "He: expected heterozygosity | Ho: observed heterozygosity | π: nucleotide diversity",
    "</p></div>"
  )
}

popup_row <- function(label, value) {
  paste0(
    "<tr>",
    "<td style='padding:2px 6px 2px 0; color:#555; font-weight:600;'>", label, "</td>",
    "<td style='padding:2px 0;'>", value, "</td>",
    "</tr>"
  )
}

#' Return a named list of colour-scale breakpoints for the He legend
he_legend_colours <- function() {
  list(
    colours = c("#D32F2F", "#F57C00", "#FBC02D", "#388E3C", "#1565C0"),
    labels  = c("< 0.10 (Critical)", "0.10–0.15", "0.15–0.20", "0.20–0.30", "> 0.30 (Healthy)")
  )
}
