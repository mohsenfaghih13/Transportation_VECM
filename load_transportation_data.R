# Load and inspect the transportation inventory dataset

# 1. Install and load the readxl package
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

# 2. Load the Excel file
# Note: the file on disk has a space before the extension
data_file <- "Final_Transportation_Inventory_Dataset .xlsx"
transport_data <- read_excel(data_file)

# 3. Display the first 10 rows
print(head(transport_data, 10))

# 4. Show column names
print(colnames(transport_data))

# 5. Show basic summary statistics
print(summary(transport_data))

# ---------------------------------------------------------------------------
# Missing values and time-series plots (TC / IC, 2003-2024)
# ---------------------------------------------------------------------------
# IC  = Total_Inventories (Census MTIS inventory levels, $ millions)
# TC  = 7 transportation service cost indexes, including Warehousing_Storage

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(ggplot2)
library(tidyr)

# 1. Check for NA values in each column
na_counts <- colSums(is.na(transport_data))
na_pct <- round(100 * na_counts / nrow(transport_data), 2)
na_table <- data.frame(
  column = names(na_counts),
  n_missing = as.integer(na_counts),
  pct_missing = as.numeric(na_pct)
)
print(na_table)

# Which years are missing for columns that have NAs
mode_cols <- c(
  "Airfreight_NonScheduled",
  "Airfreight_Scheduled",
  "Trucking_LD_LTL",
  "Trucking_LD_Truckload",
  "Trucking_Local",
  "Warehouse_Construction",
  "Warehousing_Storage"
)
for (col in c("Total_Inventories", mode_cols)) {
  missing_years <- sort(unique(transport_data$Year[is.na(transport_data[[col]])]))
  if (length(missing_years) > 0) {
    cat("\n", col, " missing in years: ", paste(missing_years, collapse = ", "), "\n", sep = "")
  }
}

# Build a monthly Date from Year and Month (Month is coded M01-M12)
month_num <- as.integer(sub("M", "", transport_data$Month))
transport_data$Date <- as.Date(paste(transport_data$Year, month_num, "01", sep = "-"))

# Average TC index across the 7 transportation modes (ignores NAs in a given month)
tc_cols <- mode_cols
transport_data$TC <- rowMeans(transport_data[, tc_cols], na.rm = TRUE)
transport_data$IC <- transport_data$Total_Inventories

# 2. Plots of TC and IC over time (2003-2024)
# 3. Linear trend overlay to visually inspect upward (non-stationary) drift

if (!dir.exists("plots")) {
  dir.create("plots")
}

# Combined TC vs IC plot (separate panels: TC is an index, IC is $ millions)
tc_ic_long <- pivot_longer(
  transport_data,
  cols = c("TC", "IC"),
  names_to = "series",
  values_to = "value"
)
tc_ic_long$series <- factor(
  tc_ic_long$series,
  levels = c("TC", "IC"),
  labels = c("TC (avg of 7 transportation modes, index)",
             "IC (Total_Inventories, $ millions)")
)

p_tc_ic <- ggplot(tc_ic_long, aes(x = Date, y = value)) +
  geom_line(color = "#1f4e79", linewidth = 0.7, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "#c45911",
              linetype = "dashed", linewidth = 0.5, na.rm = TRUE) +
  facet_wrap(~ series, scales = "free_y", ncol = 1) +
  labs(
    title = "Transportation Cost (TC) and Inventory Cost (IC), 2003-2024",
    subtitle = "Dashed lines are linear trends; an upward slope suggests non-stationarity in levels",
    x = "Year",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

print(p_tc_ic)
ggsave("plots/tc_ic_over_time.png", p_tc_ic, width = 10, height = 5.5, dpi = 150)

# IC + each of the 7 transportation modes in separate panels
plot_cols <- c("Total_Inventories", mode_cols)
plot_long <- pivot_longer(
  transport_data,
  cols = all_of(plot_cols),
  names_to = "series",
  values_to = "value"
)
plot_long$series <- factor(
  plot_long$series,
  levels = plot_cols,
  labels = c("IC (Total_Inventories)", mode_cols)
)

p_modes <- ggplot(plot_long, aes(x = Date, y = value)) +
  geom_line(color = "#1f4e79", linewidth = 0.6, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, color = "#c45911",
              linetype = "dashed", linewidth = 0.5, na.rm = TRUE) +
  facet_wrap(~ series, scales = "free_y", ncol = 2) +
  labs(
    title = "IC (Total_Inventories) and Each Transportation Mode, 2003-2024",
    subtitle = "Red dashed line = linear trend. Gaps are missing values. Upward trends suggest non-stationary levels.",
    x = "Year",
    y = "Value (IC in $ millions; modes are indexes)"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

print(p_modes)
ggsave("plots/ic_and_modes_over_time.png", p_modes, width = 11, height = 11, dpi = 150)

cat("\nPlots saved to plots/tc_ic_over_time.png and plots/ic_and_modes_over_time.png\n")

# ---------------------------------------------------------------------------
# Augmented Dickey-Fuller (ADF) tests for stationarity in levels
# ---------------------------------------------------------------------------
# H0: series has a unit root (non-stationary)
# H1: series is stationary
# tseries::adf.test uses a constant and a linear trend (appropriate given
# the upward plots). Decision rule: p < 0.05 reject H0 (stationary);
# p > 0.05 fail to reject H0 (non-stationary).

if (!requireNamespace("tseries", quietly = TRUE)) {
  install.packages("tseries")
}
library(tseries)

# Critical values from the Banerjee et al. (1993) tables used by tseries
# (ADF with constant + trend). Interpolated to the sample size of each series.
adf_critical_values <- function(n) {
  tableT <- c(25, 50, 100, 250, 500, 1e5)
  # rows: 1%, 5%, 10%
  cv_table <- rbind(
    c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
    c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41),
    c(-3.24, -3.18, -3.15, -3.13, -3.13, -3.12)
  )
  c(
    crit_1pct = unname(approx(tableT, cv_table[1, ], xout = n, rule = 2)$y),
    crit_5pct = unname(approx(tableT, cv_table[2, ], xout = n, rule = 2)$y),
    crit_10pct = unname(approx(tableT, cv_table[3, ], xout = n, rule = 2)$y)
  )
}

run_adf <- function(x, name) {
  x <- as.numeric(na.omit(x))
  test <- adf.test(x, alternative = "stationary")
  cv <- adf_critical_values(length(x))
  cat("\n========== ADF test:", name, "==========\n")
  print(test)
  cat("Critical values (constant + trend): 1% = ",
      round(cv["crit_1pct"], 2), ", 5% = ", round(cv["crit_5pct"], 2),
      ", 10% = ", round(cv["crit_10pct"], 2), "\n", sep = "")
  cat("Interpretation: p-value = ", format.pval(test$p.value, digits = 4),
      ifelse(test$p.value < 0.05,
             " < 0.05 → stationary (reject unit root)\n",
             " > 0.05 → non-stationary (fail to reject unit root)\n"),
      sep = "")
  data.frame(
    series = name,
    n = length(x),
    lags = unname(test$parameter),
    adf_statistic = unname(test$statistic),
    p_value = test$p.value,
    crit_1pct = unname(cv["crit_1pct"]),
    crit_5pct = unname(cv["crit_5pct"]),
    crit_10pct = unname(cv["crit_10pct"]),
    conclusion = ifelse(test$p.value < 0.05, "Stationary", "Non-stationary"),
    stringsAsFactors = FALSE
  )
}

# IC first, then each of the 7 transportation modes (incl. Warehousing_Storage)
adf_series <- c(
  "IC (Total_Inventories)" = "Total_Inventories",
  Airfreight_NonScheduled = "Airfreight_NonScheduled",
  Airfreight_Scheduled = "Airfreight_Scheduled",
  Trucking_LD_LTL = "Trucking_LD_LTL",
  Trucking_LD_Truckload = "Trucking_LD_Truckload",
  Trucking_Local = "Trucking_Local",
  Warehouse_Construction = "Warehouse_Construction",
  Warehousing_Storage = "Warehousing_Storage"
)

adf_results <- do.call(rbind, lapply(names(adf_series), function(lbl) {
  col <- adf_series[[lbl]]
  run_adf(transport_data[[col]], lbl)
}))
rownames(adf_results) <- NULL

cat("\n========== ADF summary (levels, 2003-2024) ==========\n")
print(adf_results, row.names = FALSE, digits = 4)

n_nonstat <- sum(adf_results$conclusion == "Non-stationary")
cat("\n", n_nonstat, " of ", nrow(adf_results),
    " series are non-stationary at the 5% level (p > 0.05).\n", sep = "")
cat("A unit root in levels is the usual justification for differencing and VECM.\n")

# ---------------------------------------------------------------------------
# First differences and ADF tests on D_IC, D_modes
# ---------------------------------------------------------------------------
# D_y_t = y_t - y_{t-1}
# If levels are I(1) and first differences are I(0), the series are
# integrated of order 1 — the VECM setting.

diff_map <- c(
  D_IC = "Total_Inventories",
  D_Airfreight_NonScheduled = "Airfreight_NonScheduled",
  D_Airfreight_Scheduled = "Airfreight_Scheduled",
  D_Trucking_LD_LTL = "Trucking_LD_LTL",
  D_Trucking_LD_Truckload = "Trucking_LD_Truckload",
  D_Trucking_Local = "Trucking_Local",
  D_Warehouse_Construction = "Warehouse_Construction",
  D_Warehousing_Storage = "Warehousing_Storage"
)

for (dcol in names(diff_map)) {
  level_col <- diff_map[[dcol]]
  transport_data[[dcol]] <- c(NA, diff(transport_data[[level_col]]))
}

cat("\nFirst differences created (D_y = y_t - y_{t-1}). Preview:\n")
print(head(transport_data[, c("Year", "Month", "Total_Inventories", "D_IC",
                             "Airfreight_Scheduled", "D_Airfreight_Scheduled")], 8))

# Plots of first differences: look for mean-reversion (stationary), not white noise.
# I(0) series fluctuate around a constant mean with no lasting trend.
# They can still be serially correlated and heteroskedastic; white noise is stronger.
diff_long <- pivot_longer(
  transport_data,
  cols = all_of(names(diff_map)),
  names_to = "series",
  values_to = "d_value"
)
diff_long$series <- factor(diff_long$series, levels = names(diff_map))

p_diffs <- ggplot(diff_long, aes(x = Date, y = d_value)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray50") +
  geom_line(color = "#1f4e79", linewidth = 0.45, na.rm = TRUE) +
  facet_wrap(~ series, scales = "free_y", ncol = 2) +
  labs(
    title = "First Differences of IC and Each Transportation Mode, 2003-2024",
    subtitle = "Gray line = 0. Stationary (I(0)) series wander around a stable mean with no upward trend. They need not look like white noise.",
    x = "Year",
    y = "First difference (D_y = y_t - y_{t-1})"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

print(p_diffs)
ggsave("plots/first_differences.png", p_diffs, width = 11, height = 11, dpi = 150)
cat("\nPlot saved to plots/first_differences.png\n")

run_adf_diff <- function(x, name) {
  x <- as.numeric(na.omit(x))
  test <- adf.test(x, alternative = "stationary")
  cat("\n========== ADF test on first difference:", name, "==========\n")
  print(test)
  cat("Interpretation: p-value = ", format.pval(test$p.value, digits = 4),
      ifelse(test$p.value < 0.05,
             " < 0.05 → first difference is stationary I(0)\n",
             " > 0.05 → first difference is still non-stationary\n"),
      sep = "")
  data.frame(
    series = name,
    n = length(x),
    lags = unname(test$parameter),
    adf_statistic = unname(test$statistic),
    p_value = test$p.value,
    conclusion = ifelse(test$p.value < 0.05,
                        "I(0) stationary",
                        "Not I(0)"),
    stringsAsFactors = FALSE
  )
}

adf_diff_results <- do.call(rbind, lapply(names(diff_map), function(dcol) {
  run_adf_diff(transport_data[[dcol]], dcol)
}))
rownames(adf_diff_results) <- NULL

cat("\n========== ADF summary (first differences) ==========\n")
print(adf_diff_results, row.names = FALSE, digits = 4)

n_i0 <- sum(adf_diff_results$conclusion == "I(0) stationary")
cat("\n", n_i0, " of ", nrow(adf_diff_results),
    " first-differenced series are stationary at the 5% level.\n", sep = "")
if (n_i0 == nrow(adf_diff_results) && n_nonstat == nrow(adf_results)) {
  cat("Levels are I(1) and differences are I(0) → all series are integrated of order 1.\n")
}
