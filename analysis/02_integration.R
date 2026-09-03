# =============================================================================
# 02_integration.R -- univariate integration testing
# =============================================================================
# Each series is tested over its own full history, not truncated to the
# common window: a unit root test needs no second series and power rises with
# T. Both windows are reported so any window-dependence is visible.
#
# The VECM stage that follows is entitled to treat a series as I(1) only if
# it appears here as I(1).
# =============================================================================

section("02  INTEGRATION TESTS")

panel <- load_panel()
series <- c(unname(CFG$modes), CFG$extra_series)
series <- series[series %in% names(panel)]

sum_rows <- list(); det_rows <- list()

for (v in series) {
  for (win in c("max", "common")) {
    ok <- !is.na(panel[[v]])
    d <- panel[ok, c("Date", v)]
    if (win == "common") d <- d[d$Date >= CFG$common_start, ]
    if (nrow(d) < CFG$min_obs) next
    assert_contiguous(d$Date, paste(v, win))

    x <- as.numeric(d[[v]])
    if (any(x <= 0)) { warning("Non-positive values in ", v, "; skipped."); next }

    res <- test_series(x, d$Date, v, win)
    sum_rows[[length(sum_rows) + 1]] <- res$summary
    det_rows[[length(det_rows) + 1]] <- cbind(
      series = v, window = win, n_obs = nrow(d), res$detail)
  }
}

tab_sum <- do.call(rbind, sum_rows)
tab_det <- do.call(rbind, det_rows)

write_out(tab_sum, "02_integration", "summary.csv")
write_out(tab_det, "02_integration", "detail.csv")

mx <- tab_sum[tab_sum$window == "max", ]

cat("\nIntegration verdicts, each series on its own maximum sample\n")
print(data.frame(
  series = mx$series, n = mx$n_obs, sample = mx$sample,
  ADF_lvl = fmt_num(mx$adf_level, 2), rej = fmt_yn(mx$adf_level_reject),
  KPSS_lvl = fmt_num(mx$kpss_level, 3), rej. = fmt_yn(mx$kpss_level_reject),
  ADF_dif = fmt_num(mx$adf_diff, 2), rej.. = fmt_yn(mx$adf_diff_reject),
  KPSS_dif = fmt_num(mx$kpss_diff, 3), rej... = fmt_yn(mx$kpss_diff_reject),
  order = mx$order, stringsAsFactors = FALSE), row.names = FALSE)

cat("\nZivot-Andrews, levels, break in intercept and trend\n")
print(data.frame(
  series = mx$series, ZA = fmt_num(mx$za_stat, 2), cv5 = fmt_num(mx$za_cv5, 2),
  break_date = mx$za_break, rejects_unit_root = fmt_yn(mx$za_reject),
  stringsAsFactors = FALSE), row.names = FALSE)

# Does the longer window change any verdict?
cm <- tab_sum[tab_sum$window == "common", c("series", "n_obs", "order")]
names(cm) <- c("series", "n_common", "order_common")
chg <- merge(mx[, c("series", "n_obs", "order")], cm, by = "series")
names(chg)[2:3] <- c("n_max", "order_max")
chg$changed <- chg$order_max != chg$order_common
write_out(chg, "02_integration", "window_comparison.csv")

cat("\nWindow dependence of the verdict\n")
print(chg, row.names = FALSE)

if (any(chg$changed, na.rm = TRUE)) {
  cat("\nVerdict changes with the longer window: ",
      paste(chg$series[which(chg$changed)], collapse = ", "), "\n", sep = "")
}

# Series the VECM stage should not treat as I(1).
bad <- mx$series[mx$order != "I(1)"]
if (length(bad)) {
  cat("\nNOT cleanly I(1), treat downstream results with caution: ",
      paste(bad, collapse = ", "), "\n", sep = "")
}

# AIC hitting the lag ceiling means a truncated search, not a converged one.
ceil <- mx$series[isTRUE(mx$lag_at_ceiling) | mx$lag_at_ceiling]
if (length(ceil)) {
  cat("AIC reached the ", CFG$lag_max, "-lag ceiling for: ",
      paste(ceil, collapse = ", "),
      "\n  -> those ADF statistics are power-deficient; read KPSS alongside.\n", sep = "")
}
