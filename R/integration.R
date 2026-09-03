# =============================================================================
# R/integration.R -- univariate unit root and stationarity testing
# =============================================================================
# ADF  H0: unit root        (urca::ur.df)
# KPSS H0: stationarity     (urca::ur.kpss)  -> confirmatory pair with ADF
# ZA   H0: unit root, against stationarity around one endogenously dated
#          break in intercept and trend (urca::ur.za)
# =============================================================================

# Lagged differences retained by ur.df after AIC selection.
adf_lags_used <- function(obj) sum(grepl("^z.diff.lag", rownames(obj@testreg$coefficients)))

adf_test <- function(x, type) {
  stat_name <- if (type == "trend") "tau3" else "tau2"
  o <- tryCatch(urca::ur.df(x, type = type, lags = CFG$lag_max, selectlags = "AIC"),
                error = function(e) e)
  if (inherits(o, "error")) {
    return(list(stat = NA_real_, cv5 = NA_real_, lags = NA_integer_, reject = NA))
  }
  st <- unname(o@teststat[1, stat_name])
  cv <- unname(o@cval[stat_name, "5pct"])
  list(stat = st, cv5 = cv, lags = adf_lags_used(o), reject = st < cv)
}

kpss_test <- function(x, type) {
  o <- tryCatch(urca::ur.kpss(x, type = type, lags = CFG$kpss_lags), error = function(e) e)
  if (inherits(o, "error")) return(list(stat = NA_real_, cv5 = NA_real_, reject = NA))
  st <- unname(o@teststat[1]); cv <- unname(o@cval[1, "5pct"])
  # H0 is stationarity, so rejection means NON-stationary
  list(stat = st, cv5 = cv, reject = st > cv)
}

za_test <- function(x, dates, lag) {
  o <- tryCatch(urca::ur.za(x, model = CFG$za_model, lag = lag), error = function(e) e)
  if (inherits(o, "error")) {
    return(list(stat = NA_real_, cv5 = NA_real_, brk = NA_character_, reject = NA))
  }
  bp <- o@bpoint
  list(stat = unname(o@teststat),
       # ur.za returns an UNNAMED vector ordered 1%, 5%, 10%
       cv5 = unname(o@cval[2]),
       brk = if (is.na(bp) || bp < 1 || bp > length(dates)) NA_character_
             else format(dates[bp], "%Y-%m"),
       reject = unname(o@teststat) < unname(o@cval[2]))
}

# ur.df has no seasonal option but the VECMs carry season = 12, so the ADF is
# also run with month-of-year means removed. This keeps the univariate verdict
# under the same seasonal treatment as the system.
seasonal_demean <- function(x, dates) {
  mth <- factor(format(dates, "%m"))
  if (nlevels(mth) < 2) return(x)
  as.numeric(residuals(stats::lm(x ~ mth))) + mean(x)
}

classify_order <- function(adf_lvl, kpss_lvl, adf_dif, kpss_dif) {
  lvl_ns <- (!isTRUE(adf_lvl)) && isTRUE(kpss_lvl)
  dif_st <- isTRUE(adf_dif) && (!isTRUE(kpss_dif))
  if (lvl_ns && dif_st) return("I(1)")
  if (isTRUE(adf_lvl) && !isTRUE(kpss_lvl)) return("I(0)")
  if (dif_st) return("I(1) (level tests disagree)")
  "inconclusive"
}

# Full battery for one series on one window.
test_series <- function(x_raw, dates, series, window) {
  lv <- log(x_raw); df <- diff(lv)
  lv_sa <- seasonal_demean(lv, dates); df_sa <- diff(lv_sa)

  a_lv    <- adf_test(lv, "trend")
  a_lv_dr <- adf_test(lv, "drift")
  a_df    <- adf_test(df, "drift")
  a_lv_sa <- adf_test(lv_sa, "trend")
  a_df_sa <- adf_test(df_sa, "drift")
  k_lv    <- kpss_test(lv, "tau")
  k_lv_mu <- kpss_test(lv, "mu")
  k_df    <- kpss_test(df, "mu")
  z       <- za_test(lv, dates, if (is.na(a_lv$lags)) CFG$za_lag else a_lv$lags)

  lags_all <- c(a_lv$lags, a_lv_sa$lags, a_df$lags, a_df_sa$lags)

  list(
    summary = data.frame(
      series = series, window = window, n_obs = length(lv),
      sample = paste(format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m"), sep = " .. "),
      adf_level = a_lv$stat, adf_level_lags = a_lv$lags, adf_level_reject = a_lv$reject,
      adf_level_sa = a_lv_sa$stat, adf_level_sa_reject = a_lv_sa$reject,
      kpss_level = k_lv$stat, kpss_level_reject = k_lv$reject,
      adf_diff = a_df$stat, adf_diff_lags = a_df$lags, adf_diff_reject = a_df$reject,
      adf_diff_sa = a_df_sa$stat, adf_diff_sa_lags = a_df_sa$lags,
      adf_diff_sa_reject = a_df_sa$reject,
      kpss_diff = k_df$stat, kpss_diff_reject = k_df$reject,
      za_stat = z$stat, za_cv5 = z$cv5, za_break = z$brk, za_reject = z$reject,
      # AIC running to the ceiling means the search was truncated, not
      # converged, and the test loses power by construction.
      lag_at_ceiling = any(lags_all == CFG$lag_max, na.rm = TRUE),
      order = classify_order(a_lv$reject, k_lv$reject, a_df$reject, k_df$reject),
      stringsAsFactors = FALSE),
    detail = rbind(
      data.frame(test = "ADF", spec = "trend", transform = "log level",
                 statistic = a_lv$stat, cv_5pct = a_lv$cv5, lags = a_lv$lags,
                 reject_H0 = a_lv$reject, H0 = "unit root"),
      data.frame(test = "ADF", spec = "drift", transform = "log level",
                 statistic = a_lv_dr$stat, cv_5pct = a_lv_dr$cv5, lags = a_lv_dr$lags,
                 reject_H0 = a_lv_dr$reject, H0 = "unit root"),
      data.frame(test = "ADF", spec = "drift", transform = "log diff",
                 statistic = a_df$stat, cv_5pct = a_df$cv5, lags = a_df$lags,
                 reject_H0 = a_df$reject, H0 = "unit root"),
      data.frame(test = "ADF", spec = "trend", transform = "log level, seas. demeaned",
                 statistic = a_lv_sa$stat, cv_5pct = a_lv_sa$cv5, lags = a_lv_sa$lags,
                 reject_H0 = a_lv_sa$reject, H0 = "unit root"),
      data.frame(test = "ADF", spec = "drift", transform = "log diff, seas. demeaned",
                 statistic = a_df_sa$stat, cv_5pct = a_df_sa$cv5, lags = a_df_sa$lags,
                 reject_H0 = a_df_sa$reject, H0 = "unit root"),
      data.frame(test = "KPSS", spec = "tau", transform = "log level",
                 statistic = k_lv$stat, cv_5pct = k_lv$cv5, lags = NA_integer_,
                 reject_H0 = k_lv$reject, H0 = "stationary"),
      data.frame(test = "KPSS", spec = "mu", transform = "log level",
                 statistic = k_lv_mu$stat, cv_5pct = k_lv_mu$cv5, lags = NA_integer_,
                 reject_H0 = k_lv_mu$reject, H0 = "stationary"),
      data.frame(test = "KPSS", spec = "mu", transform = "log diff",
                 statistic = k_df$stat, cv_5pct = k_df$cv5, lags = NA_integer_,
                 reject_H0 = k_df$reject, H0 = "stationary"),
      data.frame(test = "ZA", spec = CFG$za_model, transform = paste("log level, break", z$brk),
                 statistic = z$stat, cv_5pct = z$cv5, lags = NA_integer_,
                 reject_H0 = z$reject, H0 = "unit root")
    )
  )
}
