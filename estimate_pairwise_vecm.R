# =============================================================================
# Pairwise VECM estimation: four transportation modes vs inventories
# =============================================================================
# Specification (each bivariate system in logs):
#   VARselect K = 1..12 (constant + trend + monthly seasonals)  # FIXED: Issue #2
#   ca.jo(type = trace/eigen, ecdet = "trend", spec = "longrun",
#         season = 12, K = AIC with K >= 2)
#   VECM: cajorls(r = 1)
#   Diagnostics: equation Ljung-Box(12) and vec2var Portmanteau
#
# Models
#   1. log(Trucking_LD_LTL)       + log(IC)
#   2. log(Trucking_LD_Truckload) + log(IC)
#   3. log(Trucking_Local)        + log(IC)
#   4. log(Airfreight_Scheduled)  + log(IC)
# =============================================================================

library(readxl)
library(urca)
library(vars)
library(dplyr)

stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.01, "***",
                ifelse(p < 0.05, "**",
                       ifelse(p < 0.10, "*", ""))))
}

fmt_num <- function(x, d = 4) {
  ifelse(is.na(x), "", formatC(unname(x), format = "f", digits = d))
}

fmt_star <- function(est, p, d = 4) {
  ifelse(is.na(est), "", paste0(fmt_num(est, d), stars(p)))
}

# gap_t = (1 + alpha) * gap_{t-1}; defined only for -1 < alpha < 0
half_life <- function(alpha) {
  if (length(alpha) != 1 || is.na(alpha) || alpha >= 0 || alpha <= -1) {
    return(NA_real_)
  }
  log(0.5) / log(1 + alpha)
}

# urca stores tests as r <= 1 then r = 0 in the bivariate case
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

rank_label <- function(r) {
  if (r == 0) "r = 0 (not cointegrated)"
  else if (r == 1) "r = 1 (cointegrated)"
  else "r = 2 (full rank)"
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
    std_error = unname(cf[, "Std. Error"]),
    t_stat = unname(cf[, "t value"]),
    p_value = unname(cf[, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}

pick_ect <- function(eq_df) {
  hit <- eq_df[eq_df$term == "ect1", , drop = FALSE]
  if (!nrow(hit)) stop("No ect1 term. Terms: ", paste(eq_df$term, collapse = ", "))
  hit[1, ]
}

info_criteria <- function(resid_mat, n_params) {
  U <- as.matrix(resid_mat)
  Tobs <- nrow(U)
  sigma <- crossprod(U) / Tobs
  ldet <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
  list(
    T = Tobs,
    n_params = n_params,
    AIC = Tobs * ldet + 2 * n_params,
    SBC = Tobs * ldet + log(Tobs) * n_params
  )
}

ljung_box <- function(resid_mat, lag = 12L) {
  U <- as.matrix(resid_mat)
  bind_rows(lapply(seq_len(ncol(U)), function(j) {
    bt <- Box.test(U[, j], lag = lag, type = "Ljung-Box")
    data.frame(
      equation = colnames(U)[j],
      Q = unname(as.numeric(bt$statistic)),
      df = unname(as.numeric(bt$parameter)),
      p_value = unname(bt$p.value),
      no_ac_5pct = unname(bt$p.value) > 0.05,
      stringsAsFactors = FALSE
    )
  }))
}

hl_note <- function(p, hl) {
  if (is.na(p) || p >= 0.10) return("alpha not significant at 10%")
  if (is.na(hl)) return("alpha not in (-1, 0); half-life formula not defined")
  sprintf("%.2f months to close half the gap", hl)
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

raw <- read_excel("Transportation_Inventory_Complete_10_Modes_WPU30.xlsx")
dat <- raw %>%
  mutate(
    Date = as.Date(sprintf("%d-%02d-01", as.integer(Year), as.integer(Month))),
    log_IC = log(Total_Inventories),
    log_Trucking_LD_LTL = log(Trucking_LD_LTL),
    log_Trucking_LD_Truckload = log(Trucking_LD_Truckload),
    log_Trucking_Local = log(Trucking_Local),
    log_Airfreight_Scheduled = log(Airfreight_Scheduled)
  )

models <- list(
  list(name = "Model 1", mode = "Trucking_LD_LTL",
       mode_log = "log_Trucking_LD_LTL",
       title = "log LTL trucking + log IC"),
  list(name = "Model 2", mode = "Trucking_LD_Truckload",
       mode_log = "log_Trucking_LD_Truckload",
       title = "log truckload + log IC"),
  list(name = "Model 3", mode = "Trucking_Local",
       mode_log = "log_Trucking_Local",
       title = "log local trucking + log IC"),
  list(name = "Model 4", mode = "Airfreight_Scheduled",
       mode_log = "log_Airfreight_Scheduled",
       title = "log scheduled airfreight + log IC")
)

dir.create("vecm_output", showWarnings = FALSE)

lag_sel_rows <- list()
johansen_rows <- list()
ect_rows <- list()
beta_rows <- list()
hl_rows <- list()
diag_rows <- list()
lag_coef_rows <- list()
summary_rows <- list()
model_notes <- character()

for (m in models) {
  cat("\n##############################\n")
  cat(m$name, ":", m$title, "\n")
  cat("##############################\n")

  y <- cbind(as.numeric(dat[[m$mode_log]]), as.numeric(dat$log_IC))
  colnames(y) <- c(m$mode_log, "log_IC")
  storage.mode(y) <- "double"

  vs <- VARselect(y, lag.max = 12, type = "both", season = 12)  # FIXED: Issue #2 (was lag.max = 10)
  cat("VARselect AIC/HQ/SC/FPE:", paste(vs$selection, collapse = ", "), "\n")
  print(vs$criteria)

  aic_k <- as.integer(unname(vs$selection["AIC(n)"]))
  sbc_k <- as.integer(unname(vs$selection["SC(n)"]))
  k_use <- max(2L, min(12L, ifelse(is.na(aic_k), 2L, aic_k)))  # FIXED: Issue #2 (was min(10L,))

  lag_sel_rows[[length(lag_sel_rows) + 1]] <- data.frame(
    model = m$name, pair = m$title,
    AIC_K = aic_k, SBC_K = sbc_k, K_used = k_use,
    note = if (!is.na(aic_k) && !is.na(sbc_k) && aic_k != sbc_k) {
      "AIC and SBC disagree; AIC used (ca.jo requires K >= 2)"
    } else {
      "AIC = SBC (ca.jo requires K >= 2)"
    },
    stringsAsFactors = FALSE
  )

  jo_tr <- ca.jo(y, type = "trace", ecdet = "trend", K = k_use,
                 spec = "longrun", season = 12)
  jo_ei <- ca.jo(y, type = "eigen", ecdet = "trend", K = k_use,
                 spec = "longrun", season = 12)
  cat("\n--- Trace ---\n"); print(summary(jo_tr))
  cat("\n--- Max-eigen ---\n"); print(summary(jo_ei))

  tr <- rank_from_jo(jo_tr)
  ei <- rank_from_jo(jo_ei)

  johansen_rows[[length(johansen_rows) + 1]] <- data.frame(
    model = m$name, pair = m$title, K = k_use,
    trace_r0 = tr$stat_r0, trace_cv5_r0 = tr$cv5_r0,
    trace_r1 = tr$stat_r1, trace_cv5_r1 = tr$cv5_r1,
    rank_trace = tr$rank, rank_trace_label = rank_label(tr$rank),
    lmax_r0 = ei$stat_r0, lmax_cv5_r0 = ei$cv5_r0,
    lmax_r1 = ei$stat_r1, lmax_cv5_r1 = ei$cv5_r1,
    rank_lmax = ei$rank, rank_lmax_label = rank_label(ei$rank),
    conclusion_5pct = paste0("Trace: ", rank_label(tr$rank),
                             "; lmax: ", rank_label(ei$rank)),
    stringsAsFactors = FALSE
  )

  ols <- cajorls(jo_tr, r = 1)
  beta <- as.matrix(ols$beta)
  beta_norm <- as.numeric(beta[, 1] / beta[1, 1])
  names(beta_norm) <- rownames(beta)
  beta_mode <- unname(beta_norm[1])
  beta_ic <- unname(beta_norm[2])
  trend_nm <- grep("trend", names(beta_norm), value = TRUE)[1]
  beta_tr <- unname(beta_norm[[trend_nm]])
  lr_elast <- -beta_ic

  beta_rows[[length(beta_rows) + 1]] <- data.frame(
    model = m$name, pair = m$title,
    beta_log_Mode = beta_mode,
    beta_log_IC = beta_ic,
    beta_trend = beta_tr,
    lr_elasticity_mode_wrt_IC = lr_elast,
    interpretation = sprintf(
      "A 1%% rise in inventories is associated with a %.3f%% %s in %s in the long run.",
      abs(lr_elast),
      if (lr_elast >= 0) "increase" else "decrease",
      m$mode
    ),
    stringsAsFactors = FALSE
  )

  rlm <- ols$rlm
  mode_d <- paste0(m$mode_log, ".d")
  ic_d <- "log_IC.d"
  cat("VECM terms:\n"); print(rownames(coef(rlm)))

  eq_mode <- eq_table(rlm, mode_d)
  eq_ic <- eq_table(rlm, ic_d)
  eq_mode$equation <- m$mode_log
  eq_ic$equation <- "log_IC"
  eq_all <- bind_rows(eq_mode, eq_ic)
  eq_all$model <- m$name
  eq_all$pair <- m$title
  lag_coef_rows[[length(lag_coef_rows) + 1]] <- eq_all

  ect_m <- pick_ect(eq_mode)
  ect_i <- pick_ect(eq_ic)
  a_mode <- ect_m$estimate[[1]]
  a_ic <- ect_i$estimate[[1]]
  t_mode <- ect_m$t_stat[[1]]
  t_ic <- ect_i$t_stat[[1]]
  p_mode <- ect_m$p_value[[1]]
  p_ic <- ect_i$p_value[[1]]

  ect_rows[[length(ect_rows) + 1]] <- data.frame(
    model = m$name, pair = m$title,
    alpha_mode = a_mode, t_mode = t_mode, p_mode = p_mode, sig_mode = stars(p_mode),
    alpha_IC = a_ic, t_IC = t_ic, p_IC = p_ic, sig_IC = stars(p_ic),
    stringsAsFactors = FALSE
  )

  hl_mode <- if (!is.na(p_mode) && p_mode < 0.10) half_life(a_mode) else NA_real_
  hl_ic <- if (!is.na(p_ic) && p_ic < 0.10) half_life(a_ic) else NA_real_

  hl_rows[[length(hl_rows) + 1]] <- data.frame(
    model = m$name, pair = m$title,
    variable = c(m$mode_log, "log_IC"),
    alpha = c(a_mode, a_ic),
    t_stat = c(t_mode, t_ic),
    p_value = c(p_mode, p_ic),
    significant_10pct = c(p_mode < 0.10, p_ic < 0.10),
    half_life_months = c(hl_mode, hl_ic),
    formula_valid = c(!is.na(hl_mode), !is.na(hl_ic)),
    note = c(hl_note(p_mode, hl_mode), hl_note(p_ic, hl_ic)),
    stringsAsFactors = FALSE
  )

  abs_m <- abs(a_mode)
  abs_i <- abs(a_ic)
  faster <- if (abs_m > abs_i) m$mode_log else "log_IC"
  one_sided <- xor(p_mode < 0.05, p_ic < 0.05)
  ratio <- max(abs_m, abs_i) / max(min(abs_m, abs_i), 1e-12)
  asymmetric <- if (p_mode >= 0.05 && p_ic >= 0.05) {
    NA
  } else {
    one_sided || (p_mode < 0.05 && p_ic < 0.05 && ratio >= 2)
  }

  U <- residuals(rlm)
  colnames(U) <- c(mode_d, ic_d)
  icrit <- info_criteria(U, length(coef(rlm)))
  lb <- ljung_box(U, lag = 12L)

  v2v <- vec2var(jo_tr, r = 1)
  pt <- tryCatch(serial.test(v2v, lags.pt = 12, type = "PT.adjusted"),
                 error = function(e) NULL)
  pt_p <- if (!is.null(pt)) unname(as.numeric(pt$serial$p.value)) else NA_real_
  pt_stat <- if (!is.null(pt)) unname(as.numeric(pt$serial$statistic)) else NA_real_

  diag_rows[[length(diag_rows) + 1]] <- data.frame(
    model = m$name, pair = m$title,
    AIC = icrit$AIC, SBC = icrit$SBC,
    n_params = icrit$n_params, T = icrit$T,
    LB_Q_mode = lb$Q[[1]], LB_p_mode = lb$p_value[[1]],
    LB_no_ac_mode = lb$no_ac_5pct[[1]],
    LB_Q_IC = lb$Q[[2]], LB_p_IC = lb$p_value[[2]],
    LB_no_ac_IC = lb$no_ac_5pct[[2]],
    PT_stat = pt_stat, PT_p = pt_p,
    PT_no_ac_5pct = ifelse(is.na(pt_p), NA, pt_p > 0.05),
    stringsAsFactors = FALSE
  )

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    model = m$name, pair = m$title, K = k_use,
    rank_trace = tr$rank, rank_lmax = ei$rank,
    beta_IC = beta_ic, lr_elasticity = lr_elast,
    alpha_mode = a_mode, t_mode = t_mode, p_mode = p_mode,
    alpha_IC = a_ic, t_IC = t_ic, p_IC = p_ic,
    abs_alpha_mode = abs_m, abs_alpha_IC = abs_i,
    faster_adjuster = faster,
    asymmetric = ifelse(is.na(asymmetric), "NA", ifelse(asymmetric, "YES", "NO")),
    hl_mode_months = hl_mode, hl_IC_months = hl_ic,
    AIC = icrit$AIC, SBC = icrit$SBC,
    LB_clean = isTRUE(lb$no_ac_5pct[[1]] && lb$no_ac_5pct[[2]]),
    PT_p = pt_p,
    stringsAsFactors = FALSE
  )

  note <- sprintf(
    "%s: K(AIC)=%s (SBC=%s, used=%s). Trace r=%s, lmax r=%s. Faster |alpha|: %s. Asymmetric: %s.",
    m$name, aic_k, sbc_k, k_use, tr$rank, ei$rank, faster,
    ifelse(is.na(asymmetric), "NA", ifelse(asymmetric, "YES", "NO"))
  )
  model_notes <- c(model_notes, note)
  cat(note, "\n")
}

tab_lag <- bind_rows(lag_sel_rows)
tab_joh <- bind_rows(johansen_rows)
tab_ect <- bind_rows(ect_rows)
tab_beta <- bind_rows(beta_rows)
tab_hl <- bind_rows(hl_rows)
tab_diag <- bind_rows(diag_rows)
tab_lags <- bind_rows(lag_coef_rows)
tab_sum <- bind_rows(summary_rows)

n_models <- length(models)

# FIXED: Issue #3 — fail if per-model tables are out of sync before formatting
stopifnot(nrow(tab_joh) == n_models)
tab_joh_print <- tab_joh %>%
  transmute(
    Model = model, Pair = pair, K,
    `Trace r=0` = fmt_num(trace_r0, 3),
    `CV 5% r=0` = fmt_num(trace_cv5_r0, 2),
    `Trace r<=1` = fmt_num(trace_r1, 3),
    `CV 5% r<=1` = fmt_num(trace_cv5_r1, 2),
    `Trace rank` = rank_trace_label,
    `lmax r=0` = fmt_num(lmax_r0, 3),
    `CV 5% lmax r=0` = fmt_num(lmax_cv5_r0, 2),
    `lmax r<=1` = fmt_num(lmax_r1, 3),
    `lmax rank` = rank_lmax_label,
    Conclusion = conclusion_5pct
  )

# FIXED: Issue #3 — join only after confirming one row per model in both tables
stopifnot(nrow(tab_ect) == n_models, nrow(tab_diag) == n_models,
          nrow(tab_ect) == nrow(tab_diag))
# FIXED: Issue #1 — join on model/pair so AIC/SBC/LB are not recycled from an external vector
tab_ect_print <- tab_ect %>%
  left_join(
    tab_diag %>%
      select(model, pair, AIC, SBC, LB_p_mode, LB_p_IC, LB_no_ac_mode, LB_no_ac_IC),
    by = c("model", "pair")
  ) %>%
  transmute(
    Model = model, Pair = pair,
    `alpha_mode` = fmt_star(alpha_mode, p_mode),
    `t_mode` = fmt_num(t_mode, 3),
    `alpha_IC` = fmt_star(alpha_IC, p_IC),
    `t_IC` = fmt_num(t_IC, 3),
    AIC = fmt_num(AIC, 2),
    SBC = fmt_num(SBC, 2),
    `LB p (mode)` = fmt_num(LB_p_mode, 3),
    `LB p (IC)` = fmt_num(LB_p_IC, 3),
    `No AC (both, 5%)` = ifelse(LB_no_ac_mode & LB_no_ac_IC, "YES", "NO")
  )

stopifnot(nrow(tab_beta) == n_models)  # FIXED: Issue #3
tab_beta_print <- tab_beta %>%
  transmute(
    Model = model, Pair = pair,
    `beta_Mode` = fmt_num(beta_log_Mode, 4),
    `beta_IC` = fmt_num(beta_log_IC, 4),
    `beta_trend` = fmt_num(beta_trend, 4),
    `LR elasticity (mode wrt IC)` = fmt_num(lr_elasticity_mode_wrt_IC, 4),
    Interpretation = interpretation
  )

stopifnot(nrow(tab_hl) == 2L * n_models)  # FIXED: Issue #3 (mode + IC per model)
tab_hl_print <- tab_hl %>%
  transmute(
    Model = model, Variable = variable,
    `alpha` = fmt_star(alpha, p_value),
    `t` = fmt_num(t_stat, 3),
    `Half-life (months)` = ifelse(is.na(half_life_months), "n.a.", fmt_num(half_life_months, 2)),
    Note = note
  )

stopifnot(nrow(tab_sum) == n_models)  # FIXED: Issue #3
tab_sum_print <- tab_sum %>%
  transmute(
    Model = model, Pair = pair, K,
    `Trace r` = rank_trace, `lmax r` = rank_lmax,
    `LR elasticity` = fmt_num(lr_elasticity, 3),
    `alpha_mode` = fmt_star(alpha_mode, p_mode),
    `alpha_IC` = fmt_star(alpha_IC, p_IC),
    `|alpha_mode|` = fmt_num(abs_alpha_mode, 4),
    `|alpha_IC|` = fmt_num(abs_alpha_IC, 4),
    `Faster adjuster` = faster_adjuster,
    Asymmetric = asymmetric,
    `HL mode` = ifelse(is.na(hl_mode_months), "n.a.", fmt_num(hl_mode_months, 2)),
    `HL IC` = ifelse(is.na(hl_IC_months), "n.a.", fmt_num(hl_IC_months, 2)),
    `LB clean` = ifelse(LB_clean, "YES", "NO")
  )

tab_lag3 <- tab_lags %>%
  filter(term == "ect1" | grepl("\\.dl[123]$", term)) %>%
  transmute(
    Model = model, Equation = equation, Term = term,
    Estimate = fmt_star(estimate, p_value),
    `t-stat` = fmt_num(t_stat, 3),
    `p-value` = fmt_num(p_value, 4)
  )

cat("\n\n========== TABLE 1: Johansen tests ==========\n")
print(as.data.frame(tab_joh_print), row.names = FALSE)
cat("\n========== TABLE 2: ECT and diagnostics ==========\n")
print(as.data.frame(tab_ect_print), row.names = FALSE)
cat("\n========== TABLE 3: Cointegrating vectors (mode = 1) ==========\n")
print(as.data.frame(tab_beta_print), row.names = FALSE)
cat("\n========== TABLE 4: Half-lives ==========\n")
print(as.data.frame(tab_hl_print), row.names = FALSE)
cat("\n========== TABLE 5: Summary comparison ==========\n")
print(as.data.frame(tab_sum_print), row.names = FALSE)
cat("\n========== VECM coefficients (ECT + lags 1-3) ==========\n")
print(as.data.frame(tab_lag3), row.names = FALSE)

saveRDS(
  list(johansen = tab_joh, ect = tab_ect, beta = tab_beta, half_life = tab_hl,
       summary = tab_sum, diag = tab_diag, lags = tab_lags, lag_sel = tab_lag,
       notes = model_notes,
       johansen_print = tab_joh_print, ect_print = tab_ect_print,
       beta_print = tab_beta_print, hl_print = tab_hl_print,
       sum_print = tab_sum_print, lag3_print = tab_lag3),
  file = "vecm_output/vecm_results.rds"
)

write_tab <- function(x, path) {
  write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8")
}
write_tab(tab_joh, "vecm_output/table1_johansen.csv")
write_tab(tab_ect, "vecm_output/table2_ect.csv")
write_tab(tab_beta, "vecm_output/table3_beta.csv")
write_tab(tab_hl, "vecm_output/table4_halflife.csv")
write_tab(tab_sum, "vecm_output/table5_summary.csv")
write_tab(tab_lags, "vecm_output/vecm_all_coefficients.csv")
write_tab(tab_diag, "vecm_output/table_diagnostics.csv")
write_tab(tab_lag, "vecm_output/table_lag_selection.csv")
write_tab(tab_joh_print, "vecm_output/table1_johansen_print.csv")
write_tab(tab_ect_print, "vecm_output/table2_ect_print.csv")
write_tab(tab_beta_print, "vecm_output/table3_beta_print.csv")
write_tab(tab_hl_print, "vecm_output/table4_halflife_print.csv")
write_tab(tab_sum_print, "vecm_output/table5_summary_print.csv")
write_tab(tab_lag3, "vecm_output/table_lags_1to3_print.csv")
cat("\nCSV tables written to vecm_output/\nDone.\n")
