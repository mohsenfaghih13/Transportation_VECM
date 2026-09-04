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

# --- descriptive statistics -------------------------------------------------
# Follows the Maysami & Koh (2000) Table 2 format: mean/std dev/min/max per
# series, in levels and first differences. Computed
# in log units, since that is what actually enters every VECM here (raw-unit
# statistics would not describe the variables the models are estimated on).
# Primary system's common window (run A), the same sample the headline table
# in 04_findings.R uses.
primary_run <- Filter(function(r) isTRUE(r$primary), CFG$runs)[[1]]
dp <- run_sample(panel, primary_run)

desc_vars <- c(unname(CFG$modes), primary_run$ic)
desc_names <- c(names(CFG$modes), "IC (Census MTIS)")

desc_stats <- function(x) {
  c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x))
}

desc <- do.call(rbind, lapply(seq_along(desc_vars), function(i) {
  lv <- log(as.numeric(dp[[desc_vars[i]]]))
  df <- diff(lv)
  rbind(
    data.frame(series = desc_names[i], form = "log level",
               t(desc_stats(lv)), stringsAsFactors = FALSE),
    data.frame(series = desc_names[i], form = "log first-difference",
               t(desc_stats(df)), stringsAsFactors = FALSE))
}))

cat("\nDescriptive statistics, primary system (", primary_run$tag, "), ",
    format(min(dp$Date), "%Y-%m"), "..", format(max(dp$Date), "%Y-%m"), "\n", sep = "")
print(data.frame(series = desc$series, form = desc$form, n = desc$n,
                  mean = fmt_num(desc$mean, 4), sd = fmt_num(desc$sd, 4),
                  min = fmt_num(desc$min, 4), max = fmt_num(desc$max, 4),
                  stringsAsFactors = FALSE), row.names = FALSE)
write_out(desc, "01_data_audit", "descriptive_statistics.csv")
