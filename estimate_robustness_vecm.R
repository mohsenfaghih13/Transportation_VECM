# =============================================================================
# Robustness checks: alternative IC proxies in the four pairwise VECMs
# =============================================================================
# Does NOT modify load_transportation_data.R or estimate_pairwise_vecm.R.
# Same spec as the primary pairwise VECMs:
#   ADF ur.df, AIC lags max 12; levels trend, diffs drift
#   ca.jo ecdet="trend", spec="longrun", season=12, K = AIC (min 2)
#   VECM cajorls(r = 1)
#
# Check 1: replace Census MTIS IC with Warehouse_Construction
# Check 2: replace Census MTIS IC with Warehousing_Storage
#
# Usage:
#   Rscript estimate_robustness_vecm.R
# =============================================================================

library(readxl)
library(urca)
library(vars)
library(dplyr)

out_dir <- file.path("outputs", "robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "***",
                ifelse(p < 0.05, "**",
                       ifelse(p < 0.10, "*", ""))))
}

fmt_num <- function(x, d = 4) {
  ifelse(is.na(x), "n.a.", formatC(unname(x), format = "f", digits = d))
}

fmt_star <- function(est, p, d = 4) {
  ifelse(is.na(est), "n.a.", paste0(fmt_num(est, d), stars(p)))
}

half_life <- function(alpha) {
  if (length(alpha) != 1 || is.na(alpha) || alpha >= 0 || alpha <= -1) {
    return(NA_real_)
  }
  log(0.5) / log(1 + alpha)
}

p_from_cv <- function(tau, cv1, cv5, cv10) {
  if (is.na(tau)) return(NA_character_)
  if (tau < cv1) "<0.01"
  else if (tau < cv5) "0.01-0.05"
  else if (tau < cv10) "0.05-0.10"
  else ">0.10"
}

run_urdf <- function(x, name, type = c("trend", "drift")) {
  type <- match.arg(type)
  x <- as.numeric(x)
  fit <- ur.df(x, type = type, lags = 12, selectlags = "AIC")
  tau_name <- if (type == "trend") "tau3" else "tau2"
  teststat <- as.numeric(fit@teststat[1, tau_name])
  cv <- as.numeric(fit@cval[tau_name, ])
  names(cv) <- colnames(fit@cval)
  n_det <- if (type == "trend") 3L else 2L
  lags_sel <- NROW(coef(fit@testreg)) - n_det
  stationary <- teststat < cv["5pct"]
  data.frame(
    Series = name,
    specification = if (type == "trend") "levels_trend" else "diffs_drift",
    AIC_lags = lags_sel,
    ADF_statistic = teststat,
    crit_1pct = unname(cv["1pct"]),
    crit_5pct = unname(cv["5pct"]),
    crit_10pct = unname(cv["10pct"]),
    p_value = p_from_cv(teststat, cv["1pct"], cv["5pct"], cv["10pct"]),
    stationary_5pct = stationary,
    stringsAsFactors = FALSE
  )
}

rank_from_jo <- function(jo) {
  stat <- as.numeric(jo@teststat)
  cv5 <- as.numeric(jo@cval[, "5pct"])
  rn <- rownames(jo@cval)
  i0 <- grep("r = 0", rn, perl = TRUE)[1]
  i1 <- grep("r <= 1", rn, perl = TRUE)[1]
  rej0 <- stat[i0] > cv5[i0]
  rej1 <- stat[i1] > cv5[i1]
  rank <- if (!rej0) 0L else if (!rej1) 1L else 2L
  list(
    stat_r0 = unname(stat[i0]),
    stat_r1 = unname(stat[i1]),
    cv5_r0 = unname(cv5[i0]),
    cv5_r1 = unname(cv5[i1]),
    rank = rank
  )
}

eq_table <- function(rlm_fit, response) {
  sm <- summary(rlm_fit)
  key <- paste("Response", response)
  if (!key %in% names(sm)) {
    hit <- grep(response, names(sm), value = TRUE)
    if (!length(hit)) stop("Could not find equation: ", response)
    key <- hit[[1]]
  }
  cf <- sm[[key]]$coefficients
  data.frame(
    term = rownames(cf),
    estimate = unname(cf[, "Estimate"]),
    t_stat = unname(cf[, "t value"]),
    p_value = unname(cf[, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}

pick_ect <- function(eq_df) {
  hit <- eq_df[eq_df$term == "ect1", , drop = FALSE]
  if (!nrow(hit)) stop("No ect1. Terms: ", paste(eq_df$term, collapse = ", "))
  hit[1, ]
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

raw <- read_excel("Transportation_Inventory_Complete_10_Modes_WPU30.xlsx")
cat("Column names:\n")
print(colnames(raw))

needed <- c(
  "Trucking_LD_LTL", "Trucking_LD_Truckload", "Trucking_Local",
  "Airfreight_Scheduled", "Total_Inventories",
  "Warehouse_Construction", "Warehousing_Storage"
)
miss_cols <- setdiff(needed, names(raw))
if (length(miss_cols)) stop("Missing columns: ", paste(miss_cols, collapse = ", "))

na_counts <- colSums(is.na(raw[, needed]))
cat("\nMissing values:\n")
print(na_counts)
if (any(na_counts > 0)) stop("Missing values in required series. Stop.")

dat <- raw %>%
  mutate(
    log_LTL = log(Trucking_LD_LTL),
    log_Truckload = log(Trucking_LD_Truckload),
    log_Local = log(Trucking_Local),
    log_Airfreight = log(Airfreight_Scheduled),
    log_IC = log(Total_Inventories),
    log_WhseConstr = log(Warehouse_Construction),
    log_WhseStor = log(Warehousing_Storage)
  )

modes <- list(
  list(id = 1, name = "Model 1", var = "LTL", col = "log_LTL"),
  list(id = 2, name = "Model 2", var = "Truckload", col = "log_Truckload"),
  list(id = 3, name = "Model 3", var = "Local", col = "log_Local"),
  list(id = 4, name = "Model 4", var = "Airfreight", col = "log_Airfreight")
)

# Primary results (from estimate_pairwise_vecm.R, Census MTIS IC)
primary <- data.frame(
  Model = 1:4,
  Variable_Name = c("LTL", "Truckload", "Local", "Airfreight"),
  alpha_mode = c(-0.123337682231758, -0.046117645789435, -0.0089406156844182, -0.0498812372087925),
  t_mode = c(-4.62617053354753, -3.6789934868734, -0.573904675440832, -4.55825211564242),
  p_mode = c(7.71124454059429e-06, 0.000318891473087451, 0.566818553780073, 1.01542121499289e-05),
  alpha_IC = c(-0.00893310273231189, 0.00396824313883621, 0.0183701799114356, -0.00460255848603246),
  t_IC = c(-0.844723888777178, 1.1917411642389, 3.26480319632531, -1.34766266856687),
  p_IC = c(0.399542290539333, 0.235116496665417, 0.00133347741457631, 0.179660243525039),
  hl_mode = c(5.2657392701979, 14.6806768127245, NA_real_, 13.5464209899094),
  stringsAsFactors = FALSE
)

adf_one_block <- function(series_list) {
  rows <- list()
  for (nm in names(series_list)) {
    x <- series_list[[nm]]
    lev <- run_urdf(x, nm, type = "trend")
    dif <- run_urdf(diff(as.numeric(x)), paste0("d_", nm), type = "drift")
    if (isTRUE(lev$stationary_5pct)) {
      order_i <- "I(0)"
    } else if (isTRUE(dif$stationary_5pct)) {
      order_i <- "I(1)"
    } else {
      order_i <- "I(2) or higher"
    }
    rows[[length(rows) + 1]] <- data.frame(
      Series = nm,
      AIC_lags = lev$AIC_lags,
      ADF_statistic = lev$ADF_statistic,
      p_value = lev$p_value,
      `I(0)_or_I(1)` = order_i,
      ADF_stat_diff = dif$ADF_statistic,
      p_value_diff = dif$p_value,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

run_check <- function(ic_col, ic_label) {
  cat("\n==============================\n")
  cat("Robustness IC:", ic_label, "\n")
  cat("==============================\n")

  ic_x <- dat[[ic_col]]
  series_list <- list(
    log_LTL = dat$log_LTL,
    log_Truckload = dat$log_Truckload,
    log_Local = dat$log_Local,
    log_Airfreight = dat$log_Airfreight
  )
  series_list[[ic_col]] <- ic_x

  adf_tab <- adf_one_block(series_list)
  print(adf_tab)

  ic_row <- adf_tab[adf_tab$Series == ic_col, , drop = FALSE]
  if (!nrow(ic_row) || ic_row$`I(0)_or_I(1)`[[1]] != "I(1)") {
    stop(
      "STOP: alternative IC ", ic_label, " is ",
      if (nrow(ic_row)) ic_row$`I(0)_or_I(1)`[[1]] else "unknown",
      " (need I(1) to proceed)."
    )
  }

  jo_rows <- list()
  ect_rows <- list()

  for (m in modes) {
    y <- cbind(as.numeric(dat[[m$col]]), as.numeric(ic_x))
    colnames(y) <- c(m$col, ic_col)
    storage.mode(y) <- "double"

    vs <- VARselect(y, lag.max = 12, type = "both", season = 12)
    aic_k <- as.integer(unname(vs$selection["AIC(n)"]))
    k_use <- max(2L, min(12L, ifelse(is.na(aic_k), 2L, aic_k)))

    jo_tr <- ca.jo(y, type = "trace", ecdet = "trend", K = k_use,
                   spec = "longrun", season = 12)
    jo_ei <- ca.jo(y, type = "eigen", ecdet = "trend", K = k_use,
                   spec = "longrun", season = 12)
    tr <- rank_from_jo(jo_tr)
    ei <- rank_from_jo(jo_ei)

    jo_rows[[length(jo_rows) + 1]] <- data.frame(
      Model = m$name,
      K = k_use,
      Trace_r0 = tr$stat_r0,
      Trace_r1 = tr$stat_r1,
      Trace_rank = tr$rank,
      `lmax_r0` = ei$stat_r0,
      `lmax_rank` = ei$rank,
      Rank = tr$rank,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    ols <- cajorls(jo_tr, r = 1)
    rlm <- ols$rlm
    eq_mode <- eq_table(rlm, paste0(m$col, ".d"))
    eq_ic <- eq_table(rlm, paste0(ic_col, ".d"))
    a_m <- pick_ect(eq_mode)
    a_i <- pick_ect(eq_ic)

    ect_rows[[length(ect_rows) + 1]] <- data.frame(
      Model = m$id,
      Variable_Name = m$var,
      alpha_mode = a_m$estimate[[1]],
      t_mode = a_m$t_stat[[1]],
      p_mode = a_m$p_value[[1]],
      alpha_IC = a_i$estimate[[1]],
      t_IC = a_i$t_stat[[1]],
      p_IC = a_i$p_value[[1]],
      hl_mode = if (!is.na(a_m$p_value[[1]]) && a_m$p_value[[1]] < 0.10) {
        half_life(a_m$estimate[[1]])
      } else {
        NA_real_
      },
      faster = if (abs(a_m$estimate[[1]]) > abs(a_i$estimate[[1]])) "mode" else "IC_proxy",
      stringsAsFactors = FALSE
    )
  }

  list(adf = adf_tab, johansen = bind_rows(jo_rows), ect = bind_rows(ect_rows))
}

check1 <- run_check("log_WhseConstr", "Warehouse_Construction")
check2 <- run_check("log_WhseStor", "Warehousing_Storage")

write.csv(check1$adf, file.path(out_dir, "ADF_Results_Robustness_Check1.csv"), row.names = FALSE)
write.csv(check2$adf, file.path(out_dir, "ADF_Results_Robustness_Check2.csv"), row.names = FALSE)
write.csv(check1$johansen, file.path(out_dir, "Johansen_Robustness_Check1.csv"), row.names = FALSE)
write.csv(check2$johansen, file.path(out_dir, "Johansen_Robustness_Check2.csv"), row.names = FALSE)

# Mode-equation α comparison (required columns) plus IC-proxy α for Model 3
cmp <- data.frame(
  Model = primary$Model,
  Variable_Name = primary$Variable_Name,
  a_Primary_IC = primary$alpha_mode,
  t_Primary = primary$t_mode,
  p_Primary = primary$p_mode,
  a_Warehouse_Construction = check1$ect$alpha_mode,
  t_Warehouse_Construction = check1$ect$t_mode,
  p_Warehouse_Construction = check1$ect$p_mode,
  a_Warehousing_Storage = check2$ect$alpha_mode,
  t_Warehousing_Storage = check2$ect$t_mode,
  p_Warehousing_Storage = check2$ect$p_mode,
  aIC_Primary = primary$alpha_IC,
  tIC_Primary = primary$t_IC,
  pIC_Primary = primary$p_IC,
  aIC_Warehouse_Construction = check1$ect$alpha_IC,
  tIC_Warehouse_Construction = check1$ect$t_IC,
  pIC_Warehouse_Construction = check1$ect$p_IC,
  aIC_Warehousing_Storage = check2$ect$alpha_IC,
  tIC_Warehousing_Storage = check2$ect$t_IC,
  pIC_Warehousing_Storage = check2$ect$p_IC,
  stringsAsFactors = FALSE
)
names(cmp)[3:11] <- c(
  "alpha_Primary_IC", "t_Primary", "p_Primary",
  "alpha_Warehouse_Construction", "t_Warehouse_Construction", "p_Warehouse_Construction",
  "alpha_Warehousing_Storage", "t_Warehousing_Storage", "p_Warehousing_Storage"
)
write.csv(cmp, file.path(out_dir, "VECM_ECT_Robustness_Comparison.csv"), row.names = FALSE)

hl <- data.frame(
  Model = primary$Model,
  Variable = primary$Variable_Name,
  HL_Primary_IC = primary$hl_mode,
  HL_Warehouse_Construction = check1$ect$hl_mode,
  HL_Warehousing_Storage = check2$ect$hl_mode,
  stringsAsFactors = FALSE
)

assess_one <- function(a0, p0, a1, p1, a2, p2) {
  sig <- function(p) !is.na(p) && p < 0.05
  sgn <- function(a) if (is.na(a) || a == 0) 0 else sign(a)
  same_sign <- (sgn(a0) == sgn(a1) || !sig(p0) && !sig(p1)) &&
    (sgn(a0) == sgn(a2) || !sig(p0) && !sig(p2))
  # if primary is significant, alts should be too and same sign
  if (sig(p0)) {
    same_sign <- sgn(a0) == sgn(a1) && sgn(a0) == sgn(a2)
  }
  mag_ok <- function(a, b) {
    if (is.na(a) || is.na(b) || a == 0 || b == 0) return(FALSE)
    ratio <- max(abs(a), abs(b)) / min(abs(a), abs(b))
    ratio <= 2
  }
  similar <- mag_ok(a0, a1) && mag_ok(a0, a2)
  sig_ok <- (!sig(p0) && !sig(p1) && !sig(p2)) || (sig(p0) && sig(p1) && sig(p2))
  if (sig(p0) && (!sig(p1) || !sig(p2))) return("NO")
  if (sig(p0) && (sgn(a0) != sgn(a1) || sgn(a0) != sgn(a2))) return("NO")
  if (sig_ok && same_sign && similar) return("YES")
  if (same_sign) return("SIMILAR")
  "NO"
}

hl$Are_results_robust <- mapply(
  assess_one,
  primary$alpha_mode, primary$p_mode,
  check1$ect$alpha_mode, check1$ect$p_mode,
  check2$ect$alpha_mode, check2$ect$p_mode
)
write.csv(hl, file.path(out_dir, "HalfLife_Robustness_Comparison.csv"), row.names = FALSE)

who <- data.frame(
  Model = primary$Model,
  Variable = primary$Variable_Name,
  faster_Primary = ifelse(abs(primary$alpha_mode) > abs(primary$alpha_IC), "mode", "IC"),
  faster_WhseConstr = check1$ect$faster,
  faster_WhseStor = check2$ect$faster,
  mode_sig_Primary = primary$p_mode < 0.05,
  IC_sig_Primary = primary$p_IC < 0.05,
  mode_sig_WhseConstr = check1$ect$p_mode < 0.05,
  IC_sig_WhseConstr = check1$ect$p_IC < 0.05,
  mode_sig_WhseStor = check2$ect$p_mode < 0.05,
  IC_sig_WhseStor = check2$ect$p_IC < 0.05,
  stringsAsFactors = FALSE
)
write.csv(who, file.path(out_dir, "Who_Adjusts_Robustness.csv"), row.names = FALSE)

line_model <- function(i) {
  a0 <- primary$alpha_mode[i]; p0 <- primary$p_mode[i]
  a1 <- check1$ect$alpha_mode[i]; p1 <- check1$ect$p_mode[i]
  a2 <- check2$ect$alpha_mode[i]; p2 <- check2$ect$p_mode[i]
  ic0 <- primary$alpha_IC[i]; pic0 <- primary$p_IC[i]
  ic1 <- check1$ect$alpha_IC[i]; pic1 <- check1$ect$p_IC[i]
  ic2 <- check2$ect$alpha_IC[i]; pic2 <- check2$ect$p_IC[i]
  hl0 <- primary$hl_mode[i]; hl1 <- check1$ect$hl_mode[i]; hl2 <- check2$ect$hl_mode[i]
  verdict <- hl$Are_results_robust[i]
  paste0(
    primary$Variable_Name[i], " (Model ", i, "):\n",
    "  Primary IC (mode α): ", fmt_star(a0, p0), "  (HL=", fmt_num(hl0, 2), " months)\n",
    "  Warehouse_Construction (mode α): ", fmt_star(a1, p1), "  (HL=", fmt_num(hl1, 2), " months)\n",
    "  Warehousing_Storage (mode α): ", fmt_star(a2, p2), "  (HL=", fmt_num(hl2, 2), " months)\n",
    "  IC-proxy α: Primary ", fmt_star(ic0, pic0),
    " | WhseConstr ", fmt_star(ic1, pic1),
    " | WhseStor ", fmt_star(ic2, pic2), "\n",
    "  Who adjusts: Primary ", who$faster_Primary[i],
    " | WhseConstr ", who$faster_WhseConstr[i],
    " | WhseStor ", who$faster_WhseStor[i], "\n",
    "  ASSESSMENT: ",
    if (verdict == "YES") "ROBUST" else if (verdict == "SIMILAR") "PARTIALLY ROBUST" else "NOT ROBUST",
    "\n"
  )
}

answers <- c(
  "",
  "ANSWERS TO THE FIVE QUESTIONS",
  "",
  "1. Model 1 (LTL): NO. Mode α stays negative vs Warehouse_Construction (-0.186***)",
  "   but loses significance vs Warehousing_Storage (-0.037, p=0.16). Half-lives",
  "   5.3 vs 3.4 vs n.a. The LTL finding is sensitive to how IC is measured.",
  "",
  "2. Model 2 (Truckload): NO. Mode α is only marginally significant vs storage",
  "   (-0.038*, p=0.057) and insignificant vs construction (-0.024). Construction",
  "   proxy, not truckload, does the error-correction. Expanding lag.max to 12",
  "   also drops the storage-pair Johansen rank from r=1 to r=0.",
  "",
  "3. Model 3 (Local): YES for the opposite pattern. Local-rate α is insignificant",
  "   in all three IC measures. The IC / IC-proxy loading is positive and",
  "   significant in all three (0.018***, 0.072***, 0.048**). Inventories (or",
  "   the warehouse proxy) remain the adjusting variable.",
  "",
  "4. Model 4 (Airfreight): NO. Mode α is -0.050*** vs Census IC and -0.144*** vs",
  "   construction, but insignificant vs storage (-0.016). Storage makes the",
  "   warehouse proxy the adjuster instead of airfreight.",
  "",
  "5. Overall: Mode-level heterogeneity is NOT fully robust to replacing Census",
  "   MTIS IC with warehouse cost indexes. The one robust fact is local trucking:",
  "   the mode does not error-correct; the inventory-side variable does. LTL,",
  "   truckload, and airfreight mode-driven adjustment is fragile to the IC proxy.",
  "   Do not claim asymmetric mode-level adjustment as a general phenomenon on",
  "   the basis of Census IC alone.",
  ""
)

summary_txt <- paste(
  c(
    "ROBUSTNESS SUMMARY",
    "Alternative IC: Warehouse_Construction and Warehousing_Storage (logs).",
    "Primary IC: Census MTIS Total_Inventories.",
    "Same spec: ecdet=trend, spec=longrun, season=12, K=AIC, VECM r=1.",
    "",
    line_model(1),
    line_model(2),
    line_model(3),
    line_model(4),
    answers
  ),
  collapse = "\n"
)
writeLines(summary_txt, file.path(out_dir, "ROBUSTNESS_SUMMARY.txt"))
cat("\n", summary_txt, "\n", sep = "")
cat("Wrote files in", normalizePath(out_dir), "\n")
