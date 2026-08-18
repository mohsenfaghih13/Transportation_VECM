# Transportation-inventory VECM analysis (logged series)
# 2009-2024 data: no missing values, ready to analyze

needed <- c("readxl", "ggplot2", "tidyr", "urca", "vars")
for (pkg in needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
library(readxl)
library(ggplot2)
library(tidyr)
library(urca)
library(vars)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
# 2009-2024 data: no missing values, ready to analyze
data_file <- "Transportation_Inventory_Complete_with_WPU30.xlsx"
transport_data <- read_excel(data_file)

print(head(transport_data, 10))
print(colnames(transport_data))

# Summary statistics
# 2009-2024 data: no missing values, ready to analyze
print(summary(transport_data))

# IC  = Total_Inventories (Census MTIS, $ millions)
# TC  = TC_Aggregate_WPU30 (BLS WPU30 transportation services aggregate)
mode_cols <- c(
  "Airfreight_NonScheduled",
  "Airfreight_Scheduled",
  "Trucking_LD_LTL",
  "Trucking_LD_Truckload",
  "Trucking_Local",
  "Warehouse_Construction",
  "Warehousing_Storage"
)

month_num <- as.integer(transport_data$Month)
transport_data$Date <- as.Date(
  sprintf("%d-%02d-01", as.integer(transport_data$Year), month_num)
)
transport_data$TC <- transport_data$TC_Aggregate_WPU30
transport_data$IC <- transport_data$Total_Inventories

# ---------------------------------------------------------------------------
# Log transformation
# ---------------------------------------------------------------------------
# 2009-2024 data: no missing values, ready to analyze
transport_data$log_TC <- log(transport_data$TC)
transport_data$log_IC <- log(transport_data$IC)
for (col in mode_cols) {
  transport_data[[paste0("log_", col)]] <- log(transport_data[[col]])
}
log_mode_cols <- paste0("log_", mode_cols)
log_map <- c(
  log_TC = "log_TC",
  log_IC = "log_IC",
  setNames(log_mode_cols, log_mode_cols)
)

cat("\nLogged series created. Preview:\n")
print(head(transport_data[, c("Year", "Month", "TC", "log_TC", "IC", "log_IC")], 8))

if (!dir.exists("plots")) dir.create("plots")

# ---------------------------------------------------------------------------
# Visualizations (logged levels)
# ---------------------------------------------------------------------------
# 2009-2024 data: no missing values, ready to analyze
log_tc_ic <- pivot_longer(
  transport_data,
  cols = c("log_TC", "log_IC"),
  names_to = "series",
  values_to = "value"
)
log_tc_ic$series <- factor(
  log_tc_ic$series,
  levels = c("log_TC", "log_IC"),
  labels = c("log TC (WPU30 transportation services aggregate)",
             "log IC (Total_Inventories)")
)

p_log_tc_ic <- ggplot(log_tc_ic, aes(x = Date, y = value)) +
  geom_line(color = "#1f4e79", linewidth = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#c45911",
              linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  labs(
    title = "Logged TC and IC, 2009-2024",
    subtitle = "Natural logs. Dashed line = linear trend.",
    x = "Year",
    y = "log(value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

print(p_log_tc_ic)
ggsave("plots/log_tc_ic_over_time.png", p_log_tc_ic, width = 10, height = 5.5, dpi = 150)

log_plot_cols <- c("log_IC", "log_TC", log_mode_cols)
log_long <- pivot_longer(
  transport_data,
  cols = all_of(log_plot_cols),
  names_to = "series",
  values_to = "value"
)
log_long$series <- factor(log_long$series, levels = log_plot_cols)

p_log_modes <- ggplot(log_long, aes(x = Date, y = value)) +
  geom_line(color = "#1f4e79", linewidth = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#c45911",
              linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ series, scales = "free_y", ncol = 2) +
  labs(
    title = "Logged IC, Logged TC Aggregate, and Logged Modes, 2009-2024",
    subtitle = "Natural logs. Red dashed line = linear trend. ADF uses constant + trend.",
    x = "Year",
    y = "log(value)"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

print(p_log_modes)
ggsave("plots/log_ic_tc_modes.png", p_log_modes, width = 11, height = 12, dpi = 150)
cat("\nLog-level plots saved.\n")

# ---------------------------------------------------------------------------
# ADF tests (ur.df on logged data)
# ---------------------------------------------------------------------------
# 2009-2024 data: no missing values, ready to analyze
# Log levels: type = "trend" (constant + linear trend)
# Log differences: type = "drift" (constant only)
# Lags: AIC, maximum 12. Reject H0 at 5% if tau < 5% critical value.

run_urdf <- function(x, name, type = c("trend", "drift")) {
  type <- match.arg(type)
  x <- as.numeric(x)
  fit <- ur.df(x, type = type, lags = 12, selectlags = "AIC")
  tau_name <- if (type == "trend") "tau3" else "tau2"
  spec_label <- if (type == "trend") {
    "trend (constant + linear trend)"
  } else {
    "drift (constant only)"
  }
  teststat <- as.numeric(fit@teststat[1, tau_name])
  cv <- as.numeric(fit@cval[tau_name, ])
  names(cv) <- colnames(fit@cval)
  n_det <- if (type == "trend") 3L else 2L
  lags_sel <- NROW(coef(fit@testreg)) - n_det
  conclusion <- ifelse(teststat < cv["5pct"], "Stationary", "Non-stationary")

  cat("\n========== ADF (ur.df):", name, "==========\n")
  cat("Specification:", spec_label, "\n")
  print(summary(fit))
  cat("tau statistic = ", round(teststat, 4),
      " | 1% CV = ", round(cv["1pct"], 2),
      " | 5% CV = ", round(cv["5pct"], 2),
      " | 10% CV = ", round(cv["10pct"], 2), "\n", sep = "")
  cat("Interpretation: ",
      ifelse(teststat < cv["5pct"],
             "tau < 5% CV → reject unit root (stationary)\n",
             "tau > 5% CV → fail to reject unit root (non-stationary)\n"),
      sep = "")

  data.frame(
    series = name,
    specification = spec_label,
    n = length(x),
    lags_aic = lags_sel,
    tau_stat = teststat,
    crit_1pct = unname(cv["1pct"]),
    crit_5pct = unname(cv["5pct"]),
    crit_10pct = unname(cv["10pct"]),
    conclusion = conclusion
  )
}

adf_log_levels <- do.call(rbind, lapply(names(log_map), function(lbl) {
  run_urdf(transport_data[[log_map[[lbl]]]], lbl, type = "trend")
}))
rownames(adf_log_levels) <- NULL

cat("\n========== ADF summary (LOG levels, type = trend) ==========\n")
print(adf_log_levels, row.names = FALSE, digits = 4)
n_nonstat <- sum(adf_log_levels$conclusion == "Non-stationary")
cat("\n", n_nonstat, " of ", nrow(adf_log_levels),
    " logged series are non-stationary at the 5% level.\n", sep = "")

# First differences of logs: D_log_y = log(y_t) - log(y_{t-1})
diff_log_map <- c(
  D_log_TC = "log_TC",
  D_log_IC = "log_IC",
  D_log_Airfreight_NonScheduled = "log_Airfreight_NonScheduled",
  D_log_Airfreight_Scheduled = "log_Airfreight_Scheduled",
  D_log_Trucking_LD_LTL = "log_Trucking_LD_LTL",
  D_log_Trucking_LD_Truckload = "log_Trucking_LD_Truckload",
  D_log_Trucking_Local = "log_Trucking_Local",
  D_log_Warehouse_Construction = "log_Warehouse_Construction",
  D_log_Warehousing_Storage = "log_Warehousing_Storage"
)

dlog_data <- data.frame(Date = transport_data$Date[-1])
for (dcol in names(diff_log_map)) {
  dlog_data[[dcol]] <- diff(transport_data[[diff_log_map[[dcol]]]])
}

cat("\nFirst differences of logs created. Preview:\n")
print(head(dlog_data[, c("Date", "D_log_IC", "D_log_TC")], 8))

dlog_long <- pivot_longer(
  dlog_data,
  cols = all_of(names(diff_log_map)),
  names_to = "series",
  values_to = "d_value"
)
dlog_long$series <- factor(dlog_long$series, levels = names(diff_log_map))

p_dlog <- ggplot(dlog_long, aes(x = Date, y = d_value)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray50") +
  geom_line(color = "#1f4e79", linewidth = 0.45) +
  facet_wrap(~ series, scales = "free_y", ncol = 2) +
  labs(
    title = "First Differences of Logged Series, 2009-2024",
    subtitle = "D_log_y = log(y_t) - log(y_{t-1}). Gray line = 0. ADF uses constant only (drift).",
    x = "Year",
    y = "First difference of log"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

print(p_dlog)
ggsave("plots/dlog_first_differences.png", p_dlog, width = 11, height = 12, dpi = 150)
cat("\nPlot saved to plots/dlog_first_differences.png\n")

adf_log_diffs <- do.call(rbind, lapply(names(diff_log_map), function(dcol) {
  run_urdf(dlog_data[[dcol]], dcol, type = "drift")
}))
rownames(adf_log_diffs) <- NULL

cat("\n========== ADF summary (first differences of LOGS, type = drift) ==========\n")
print(adf_log_diffs, row.names = FALSE, digits = 4)
n_i0 <- sum(adf_log_diffs$conclusion == "Stationary")
cat("\n", n_i0, " of ", nrow(adf_log_diffs),
    " first-differenced log series are stationary at the 5% level.\n", sep = "")
if (n_i0 == nrow(adf_log_diffs) && n_nonstat == nrow(adf_log_levels)) {
  cat("Log levels are I(1) and log differences are I(0) → series are integrated of order 1 in logs.\n")
}

# ---------------------------------------------------------------------------
# Seasonal dummies
# ---------------------------------------------------------------------------
# 2009-2024 data: no missing values, ready to analyze
# ca.jo(season = 12) uses centered monthly dummies (11 independent columns).
transport_data$month_f <- factor(month_num, levels = 1:12)
season_mm <- model.matrix(~ month_f, data = transport_data)
cat("\n========== Seasonal dummies (month 1-12) ==========\n")
print(head(data.frame(Date = transport_data$Date, month = month_num, season_mm), 5))

# ---------------------------------------------------------------------------
# Johansen tests (ca.jo)
# ---------------------------------------------------------------------------
# 2009-2024 data: no missing values, ready to analyze
# ecdet = "trend": constant + trend in the cointegrating relation
# season = 12; K = AIC lag of VAR in levels (min 2, max 12)
# H0: r = 0 vs r >= 1; H0: r <= 1 vs r = 2

select_var_k <- function(dat) {
  vs <- VARselect(dat, lag.max = 12, type = "both", season = 12)
  cat("VAR lag selection (AIC/HQ/SC/FPE):", paste(vs$selection, collapse = ", "), "\n")
  k <- as.integer(unname(vs$selection["AIC(n)"]))
  if (is.na(k) || k < 2) 2L else k
}

rank_conclusion <- function(stat_r0, cv5_r0, stat_r1, cv5_r1) {
  rej0 <- stat_r0 > cv5_r0
  rej1 <- stat_r1 > cv5_r1
  if (!rej0) {
    "r = 0 (not cointegrated)"
  } else if (!rej1) {
    "r = 1 (cointegrated)"
  } else {
    "r = 2 (full rank)"
  }
}

pick_row <- function(rn, pattern) {
  i <- grep(pattern, rn, perl = TRUE)
  if (length(i) != 1) stop("Could not identify row: ", pattern)
  i
}

run_johansen_pair <- function(pair_label, y_ic, y_other, other_name, dates) {
  dat <- cbind(log_IC = as.numeric(y_ic), other = as.numeric(y_other))
  colnames(dat) <- c("log_IC", other_name)
  storage.mode(dat) <- "double"
  k <- select_var_k(dat)

  jo_tr <- ca.jo(dat, type = "trace", ecdet = "trend", K = k,
                 spec = "longrun", season = 12)
  jo_ei <- ca.jo(dat, type = "eigen", ecdet = "trend", K = k,
                 spec = "longrun", season = 12)

  cat("\n========== Johansen:", pair_label, "==========\n")
  cat("Deterministic: ecdet = 'trend'; season = 12; K(AIC) =", k, "\n")
  cat("\n--- Trace ---\n")
  print(summary(jo_tr))
  cat("\n--- Maximum eigenvalue ---\n")
  print(summary(jo_ei))

  rn_tr <- rownames(jo_tr@cval)
  rn_ei <- rownames(jo_ei@cval)
  i_tr0 <- pick_row(rn_tr, "r = 0")
  i_tr1 <- pick_row(rn_tr, "r <= 1")
  i_ei0 <- pick_row(rn_ei, "r = 0")
  i_ei1 <- pick_row(rn_ei, "r <= 1")

  tr_stat <- as.numeric(jo_tr@teststat)
  ei_stat <- as.numeric(jo_ei@teststat)
  tr_cv <- jo_tr@cval
  ei_cv <- jo_ei@cval

  trace_rank <- rank_conclusion(
    tr_stat[i_tr0], tr_cv[i_tr0, "5pct"],
    tr_stat[i_tr1], tr_cv[i_tr1, "5pct"]
  )
  eigen_rank <- rank_conclusion(
    ei_stat[i_ei0], ei_cv[i_ei0, "5pct"],
    ei_stat[i_ei1], ei_cv[i_ei1, "5pct"]
  )

  cat("Trace rank conclusion (5%):        ", trace_rank, "\n")
  cat("Max-eigenvalue rank conclusion (5%): ", eigen_rank, "\n")

  vec <- jo_tr@V[, 1]
  names(vec) <- rownames(jo_tr@V)
  ic_term <- grep("^log_IC", names(vec), value = TRUE)[1]
  other_term <- grep(paste0("^", other_name), names(vec), value = TRUE)[1]
  trend_term <- grep("trend", names(vec), value = TRUE)[1]
  vec_norm <- vec / vec[[ic_term]]
  cat("Cointegrating vector (raw, col 1):\n")
  print(round(vec, 4))
  cat("Normalized on log_IC = 1:\n")
  print(round(vec_norm, 4))

  ect <- as.numeric(jo_tr@ZK %*% jo_tr@V[, 1, drop = FALSE])
  ect_dates <- tail(dates, length(ect))

  list(
    table = data.frame(
      pair = pair_label,
      K_aic = k,
      trace_r0 = tr_stat[i_tr0],
      trace_r1 = tr_stat[i_tr1],
      trace_cv5_r0 = unname(tr_cv[i_tr0, "5pct"]),
      trace_cv5_r1 = unname(tr_cv[i_tr1, "5pct"]),
      trace_cv1_r0 = unname(tr_cv[i_tr0, "1pct"]),
      trace_cv10_r0 = unname(tr_cv[i_tr0, "10pct"]),
      eigen_r0 = ei_stat[i_ei0],
      eigen_r1 = ei_stat[i_ei1],
      eigen_cv5_r0 = unname(ei_cv[i_ei0, "5pct"]),
      eigen_cv5_r1 = unname(ei_cv[i_ei1, "5pct"]),
      eigen_cv1_r0 = unname(ei_cv[i_ei0, "1pct"]),
      eigen_cv10_r0 = unname(ei_cv[i_ei0, "10pct"]),
      rank_trace = trace_rank,
      rank_maxeigen = eigen_rank,
      coint_at_5pct = grepl("r = 1", trace_rank) || grepl("r = 2", trace_rank),
      beta_log_IC = unname(vec_norm[[ic_term]]),
      beta_other = unname(vec_norm[[other_term]]),
      beta_trend = if (length(trend_term)) unname(vec_norm[[trend_term]]) else NA_real_
    ),
    lambda = data.frame(
      pair = pair_label,
      component = paste0("lambda", seq_along(jo_tr@lambda)),
      eigenvalue = as.numeric(jo_tr@lambda)
    ),
    ect = data.frame(pair = pair_label, Date = ect_dates, ect = ect),
    vector = data.frame(
      pair = pair_label,
      term = names(vec_norm),
      beta_normalized = as.numeric(vec_norm)
    )
  )
}

johansen_pairs <- list(
  list(label = "IC vs TC (WPU30 aggregate)", other = "log_TC", pretty = "log_TC"),
  list(label = "IC vs Airfreight_NonScheduled", other = "log_Airfreight_NonScheduled", pretty = "log_Airfreight_NonScheduled"),
  list(label = "IC vs Airfreight_Scheduled", other = "log_Airfreight_Scheduled", pretty = "log_Airfreight_Scheduled"),
  list(label = "IC vs Trucking_LD_LTL", other = "log_Trucking_LD_LTL", pretty = "log_Trucking_LD_LTL"),
  list(label = "IC vs Trucking_LD_Truckload", other = "log_Trucking_LD_Truckload", pretty = "log_Trucking_LD_Truckload"),
  list(label = "IC vs Trucking_Local", other = "log_Trucking_Local", pretty = "log_Trucking_Local"),
  list(label = "IC vs Warehouse_Construction", other = "log_Warehouse_Construction", pretty = "log_Warehouse_Construction"),
  list(label = "IC vs Warehousing_Storage", other = "log_Warehousing_Storage", pretty = "log_Warehousing_Storage")
)

johansen_out <- lapply(johansen_pairs, function(p) {
  run_johansen_pair(
    pair_label = p$label,
    y_ic = transport_data$log_IC,
    y_other = transport_data[[p$other]],
    other_name = p$pretty,
    dates = transport_data$Date
  )
})

johansen_table <- do.call(rbind, lapply(johansen_out, `[[`, "table"))
johansen_lambda <- do.call(rbind, lapply(johansen_out, `[[`, "lambda"))
johansen_ect <- do.call(rbind, lapply(johansen_out, `[[`, "ect"))
johansen_vec <- do.call(rbind, lapply(johansen_out, `[[`, "vector"))
rownames(johansen_table) <- NULL

cat("\n========== Johansen summary table (5% rank tests) ==========\n")
print(johansen_table[, c("pair", "K_aic",
                         "trace_r0", "trace_cv5_r0", "trace_r1", "trace_cv5_r1",
                         "eigen_r0", "eigen_cv5_r0", "eigen_r1", "eigen_cv5_r1",
                         "rank_trace", "rank_maxeigen")],
      row.names = FALSE, digits = 4)

cat("\nCointegrated pairs (trace test, 5%): r >= 1\n")
print(johansen_table$pair[johansen_table$coint_at_5pct])

cat("\nNormalized cointegrating vectors (log_IC = 1):\n")
print(johansen_vec, row.names = FALSE, digits = 4)

write.csv(johansen_table, "plots/johansen_summary.csv", row.names = FALSE)
write.csv(johansen_vec, "plots/johansen_vectors.csv", row.names = FALSE)

p_eigs <- ggplot(johansen_lambda, aes(x = component, y = eigenvalue)) +
  geom_col(fill = "#1f4e79", width = 0.7) +
  facet_wrap(~ pair, ncol = 2) +
  labs(
    title = "Johansen Eigenvalues by Pair (logged series)",
    subtitle = "ecdet = trend, season = 12. Larger lambda_1 supports a cointegrating relation.",
    x = "Eigenvalue",
    y = "lambda"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 8))

print(p_eigs)
ggsave("plots/johansen_eigenvalues.png", p_eigs, width = 11, height = 12, dpi = 150)

p_ect <- ggplot(johansen_ect, aes(x = Date, y = ect)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray50") +
  geom_line(color = "#1f4e79", linewidth = 0.45) +
  facet_wrap(~ pair, scales = "free_y", ncol = 2) +
  labs(
    title = "Cointegrating Relations (first eigenvector), 2009-2024",
    subtitle = "beta' Y_t from ca.jo (ecdet = trend). Mean-reverting path is consistent with r >= 1.",
    x = "Year",
    y = "Cointegrating combination"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 8))

print(p_ect)
ggsave("plots/johansen_coint_relations.png", p_ect, width = 11, height = 12, dpi = 150)
cat("\nJohansen plots saved to plots/johansen_eigenvalues.png and plots/johansen_coint_relations.png\n")
