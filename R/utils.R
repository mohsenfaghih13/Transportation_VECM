# =============================================================================
# R/utils.R -- formatting and small numeric helpers
# =============================================================================

stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.01, "***",
      ifelse(p < 0.05, "**",
        ifelse(p < 0.10, "*", ""))))
}

fmt_num <- function(x, d = 3) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
}

fmt_yn <- function(x) ifelse(is.na(x), "--", ifelse(x, "yes", "no"))

fmt_star <- function(est, p, d = 4) {
  ifelse(is.na(est), "--", paste0(formatC(est, format = "f", digits = d), stars(p)))
}

# Months for half of a disequilibrium to dissipate. Defined only for a
# convergent loading, -1 < alpha < 0; anything else is not a half-life and is
# returned as NA rather than as a large positive number.
half_life <- function(alpha) {
  if (length(alpha) != 1 || is.na(alpha) || alpha >= 0 || alpha <= -1) return(NA_real_)
  log(0.5) / log(1 + alpha)
}

clamp_k <- function(k, lo, hi) {
  if (is.na(k)) return(as.integer(lo))
  max(as.integer(lo), min(as.integer(hi), as.integer(k)))
}

# Ljung-Box test for leftover autocorrelation in one equation's residuals.
# H0: the residuals are white noise (no autocorrelation left). Rejecting H0
# means the model missed real structure and the equation's t-stats/p-values
# should not be trusted, regardless of how significant they look. lag = 12
# matches the monthly seasonal cycle (CFG$season); fitdf = 0 since this is a
# simple diagnostic check, not adjusted for the VECM's own parameter count.
ljung_box <- function(resid_vec, lag = 12) {
  n <- length(resid_vec)
  use_lag <- max(1, min(lag, floor(n / 5)))
  o <- tryCatch(stats::Box.test(resid_vec, lag = use_lag, type = "Ljung-Box", fitdf = 0),
                error = function(e) NULL)
  if (is.null(o)) return(list(stat = NA_real_, p = NA_real_, lag = use_lag))
  list(stat = unname(o$statistic), p = unname(o$p.value), lag = use_lag)
}

write_out <- function(df, ...) {
  path <- file.path(CFG$out_dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

section <- function(title) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
}
