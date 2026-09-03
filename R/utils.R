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

write_out <- function(df, ...) {
  path <- file.path(CFG$out_dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, path, row.names = FALSE)
  invisible(path)
}

section <- function(title) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
}
