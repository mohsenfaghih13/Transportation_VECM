# =============================================================================
# R/data.R -- loading, validation and windowing
# =============================================================================

# Read the workbook and build a Date column. Fails loudly on the two things
# that have silently corrupted this dataset before: month labels that do not
# parse to 1..12, and interior gaps in a series' history.
load_panel <- function(path = CFG$data_file) {
  raw <- readxl::read_excel(path)

  mth <- suppressWarnings(as.integer(raw$Month))
  yr  <- suppressWarnings(as.integer(raw$Year))
  if (any(is.na(mth)) || any(mth < 1 | mth > 12)) {
    bad <- unique(raw$Month[is.na(mth) | mth < 1 | mth > 12])
    stop("Month column does not parse to 1..12. Offending values: ",
         paste(utils::head(bad, 10), collapse = ", "))
  }
  if (any(is.na(yr))) stop("Year column does not parse to integers.")

  raw$Date <- as.Date(sprintf("%d-%02d-01", yr, mth))
  raw <- raw[order(raw$Date), ]

  if (anyDuplicated(raw$Date)) {
    stop("Duplicate dates: ", paste(format(raw$Date[duplicated(raw$Date)]), collapse = ", "))
  }
  raw
}

month_index <- function(d) as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))

assert_contiguous <- function(dates, what) {
  if (length(dates) < 2) return(invisible(TRUE))
  g <- diff(month_index(dates))
  if (any(g != 1L)) {
    stop("Non-contiguous monthly sample in ", what, ": gap after ",
         format(dates[which(g != 1L)[1]], "%Y-%m"))
  }
  invisible(TRUE)
}

# Coverage of every numeric series: where it starts, ends, and whether it has
# interior holes. Interior holes matter because complete-case filtering would
# silently split the sample.
series_coverage <- function(panel) {
  num <- setdiff(names(panel)[vapply(panel, is.numeric, logical(1))], c("Year", "Month"))
  do.call(rbind, lapply(num, function(v) {
    ok <- which(!is.na(panel[[v]]))
    if (!length(ok)) return(NULL)
    span <- seq(min(ok), max(ok))
    data.frame(series = v, n_obs = length(ok),
               start = format(panel$Date[min(ok)], "%Y-%m"),
               end = format(panel$Date[max(ok)], "%Y-%m"),
               interior_gaps = sum(is.na(panel[[v]][span])),
               stringsAsFactors = FALSE)
  }))
}

# Complete-case sample for one run. The realised start is reported separately
# from the requested one, because the binding constraint is usually the
# inventory series rather than the date. run$end is optional -- only the
# pre-pandemic subsample run uses it; every other run runs through the end
# of the panel.
run_sample <- function(panel, run, modes = CFG$modes) {
  need <- c(unname(modes), run$ic)
  if (isTRUE(run$deflate)) need <- c(need, CFG$deflator_series)
  miss <- setdiff(need, names(panel))
  if (length(miss)) stop("Missing columns for run ", run$tag, ": ", paste(miss, collapse = ", "))

  d <- panel[panel$Date >= as.Date(run$start), , drop = FALSE]
  if (!is.null(run$end)) d <- d[d$Date <= as.Date(run$end), , drop = FALSE]
  keep <- Reduce(`&`, lapply(need, function(v) !is.na(d[[v]])))
  d <- d[keep, , drop = FALSE]
  assert_contiguous(d$Date, run$tag)
  if (nrow(d) < CFG$min_obs) stop("Run ", run$tag, " has only ", nrow(d), " observations.")
  d
}

# Bivariate log system for one mode against the run's inventory series. When
# run$deflate is TRUE, both series are divided by CFG$deflator_series first
# (Waller's real/deflated specification, 2026-08-26): a shared inflation
# trend in a nominal inventory stock and a price index can manufacture
# cointegration on its own, so this guards against that by removing the
# common price-level trend from both sides before logging.
build_system <- function(d, mode_name, run) {
  mcol <- unname(CFG$modes[[mode_name]])
  xm <- as.numeric(d[[mcol]]); xi <- as.numeric(d[[run$ic]])
  suffix <- ""
  if (isTRUE(run$deflate)) {
    defl <- as.numeric(d[[CFG$deflator_series]])
    if (any(defl <= 0)) stop("Non-positive deflator values in ", CFG$deflator_series)
    xm <- xm / defl; xi <- xi / defl
    suffix <- "_real"
  }
  if (any(xm <= 0) || any(xi <= 0)) {
    stop("Non-positive values prevent log transform: ", mcol, " / ", run$ic)
  }
  y <- cbind(log(xm), log(xi))
  colnames(y) <- c(paste0("log_", mode_name, suffix), paste0("log_", run$ic, suffix))
  storage.mode(y) <- "double"
  stopifnot(all(is.finite(y)))
  y
}

# General log system for an arbitrary named set of source columns. cols is a
# named character vector, name -> column in the panel; used by the
# combined (3+ variable) systems in R/grid_combined.R. build_system() above
# stays the dedicated 2-variable path for the pairwise grid.
build_system_n <- function(d, cols) {
  stopifnot(!is.null(names(cols)), all(nzchar(names(cols))))
  x <- lapply(cols, function(cn) as.numeric(d[[cn]]))
  bad <- names(x)[vapply(x, function(v) any(v <= 0), logical(1))]
  if (length(bad)) {
    stop("Non-positive values prevent log transform: ", paste(cols[bad], collapse = ", "))
  }
  y <- do.call(cbind, lapply(x, log))
  colnames(y) <- names(cols)
  storage.mode(y) <- "double"
  stopifnot(all(is.finite(y)))
  y
}
