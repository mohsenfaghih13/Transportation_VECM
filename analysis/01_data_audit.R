# =============================================================================
# 01_data_audit.R -- what the data can and cannot support
# =============================================================================
# Run first. Establishes the coverage constraint that governs every later
# stage: the primary pairwise systems are limited to the overlap between the
# modes and Census MTIS inventories, regardless of the requested start date.
# =============================================================================

section("01  DATA AUDIT")

panel <- load_panel()
cat("Workbook:", CFG$data_file, "\n")
cat("Range:", format(min(panel$Date), "%Y-%m"), "..", format(max(panel$Date), "%Y-%m"),
    " rows:", nrow(panel), "\n\n")

cov <- series_coverage(panel)
cat("Series coverage\n")
print(cov, row.names = FALSE)
write_out(cov, "01_data_audit", "series_coverage.csv")

if (any(cov$interior_gaps > 0)) {
  warning("Interior gaps present in: ",
          paste(cov$series[cov$interior_gaps > 0], collapse = ", "))
}

# What each configured run actually resolves to.
res <- do.call(rbind, lapply(CFG$runs, function(run) {
  d <- run_sample(panel, run)
  data.frame(run = run$tag, label = run$label,
             requested_start = run$start, ic_series = run$ic,
             actual_start = format(min(d$Date), "%Y-%m"),
             actual_end = format(max(d$Date), "%Y-%m"),
             n_months = nrow(d),
             binding_constraint = if (as.Date(run$start) < min(d$Date))
               paste0(run$ic, " starts ", format(min(d$Date), "%Y-%m")) else "requested start",
             stringsAsFactors = FALSE)
}))

cat("\nRealised estimation samples\n")
print(res, row.names = FALSE)
write_out(res, "01_data_audit", "run_samples.csv")

cat("\nNote: where binding_constraint names a series rather than the requested\n",
    "start, moving the start date earlier will not add a single observation.\n", sep = "")
