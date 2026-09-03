# =============================================================================
# R/grid.R -- the specification grid engine
# =============================================================================
# One function, run_grid(), produces every cointegration result in the
# project. What used to be four separate phase scripts are slices of its
# output:
#   Phase 1  run A, dummy design "none"
#   Phase 2  all runs, designs "none" and "orig"
#   Phase 4  all runs, all five designs
#
# Invariants enforced here:
#   * the cointegrating rank is read from the trace test, never assumed
#   * the VECM is estimated only where that rank is estimable, and the reason
#     is recorded otherwise
#   * ecdet and lag order are grid dimensions, not constants
#   * dummy columns that cannot be identified at a given lag are dropped and
#     reported, not left to produce a singular system
#   * every dummy row carries cv_valid = FALSE, because urca's Johansen
#     critical values are not adjusted for dumvar
# =============================================================================

lag_selection <- function(y, dum) {
  vs <- if (is.null(dum)) {
    vars::VARselect(y, lag.max = CFG$lag_max, type = "both", season = CFG$season)
  } else {
    vars::VARselect(y, lag.max = CFG$lag_max, type = "both", season = CFG$season, exogen = dum)
  }
  pick <- function(nm) suppressWarnings(as.integer(unname(vs$selection[nm])))
  list(AIC = pick("AIC(n)"), SBC = pick("SC(n)"),
       HQ = pick("HQ(n)"), FPE = pick("FPE(n)"))
}

# Implied rank at each conventional level, for one Johansen object.
rank_from_jo <- function(jo) {
  stat <- as.numeric(jo@teststat); cv <- jo@cval; rn <- rownames(cv)
  i0 <- grep("r = 0", rn)[1]; i1 <- grep("r <= 1", rn)[1]
  rank_at <- function(l) {
    if (!(stat[i0] > cv[i0, l])) 0L else if (!(stat[i1] > cv[i1, l])) 1L else 2L
  }
  out <- list(stat_r0 = unname(stat[i0]), stat_r1 = unname(stat[i1]),
              cv_r0 = cv[i0, CFG$cv_levels], cv_r1 = cv[i1, CFG$cv_levels])
  for (l in CFG$cv_levels) out[[paste0("rank_", l)]] <- rank_at(l)
  out
}

eq_table <- function(rlm_fit, response) {
  sm <- summary(rlm_fit); key <- paste("Response", response)
  if (!key %in% names(sm)) {
    hit <- grep(response, names(sm), value = TRUE, fixed = TRUE)
    if (!length(hit)) stop("No equation for ", response)
    key <- hit[[1]]
  }
  cf <- sm[[key]]$coefficients
  out <- data.frame(term = rownames(cf), estimate = unname(cf[, "Estimate"]),
                     se = unname(cf[, "Std. Error"]),
                     t_stat = unname(cf[, "t value"]), p_value = unname(cf[, "Pr(>|t|)"]),
                     stringsAsFactors = FALSE)
  # summary.lm's $df is c(rank, residual df, total params); residual df is
  # what a CI on a coefficient needs.
  attr(out, "resid_df") <- unname(sm[[key]]$df[2])
  out
}

# 95% CI on one coefficient, from its point estimate/SE and the equation's
# residual degrees of freedom (not a normal approximation).
ci95 <- function(eq, nm) {
  est <- grab(eq, nm, "estimate"); se <- grab(eq, nm, "se")
  df <- attr(eq, "resid_df")
  if (is.na(est) || is.na(se) || is.null(df) || is.na(df)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  half <- stats::qt(0.975, df = df) * se
  c(lo = est - half, hi = est + half)
}

grab <- function(eq, nm, field) {
  h <- eq[eq$term == nm, , drop = FALSE]
  if (!nrow(h)) NA_real_ else h[[field]][1]
}

# Estimate at the rank the data gave. For a bivariate system urca::cajorls
# admits r = 1 only, so r = 0 and r = 2 are recorded with their reason rather
# than silently forced to 1.
fit_vecm <- function(jo, r, mode_col, ic_col) {
  if (is.na(r)) return(list(ok = FALSE, r = NA_integer_, status = "rank undetermined"))
  if (r < 1) return(list(ok = FALSE, r = r, status = "r = 0: no cointegration"))
  if (r >= ncol(jo@x)) return(list(ok = FALSE, r = r, status = sprintf("r = %d: full rank", r)))

  fit <- tryCatch(urca::cajorls(jo, r = r), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, r = r, status = paste("cajorls error:", conditionMessage(fit))))
  }
  eq_m <- eq_table(fit$rlm, paste0(mode_col, ".d"))
  eq_i <- eq_table(fit$rlm, paste0(ic_col, ".d"))
  if (!"ect1" %in% eq_m$term) return(list(ok = FALSE, r = r, status = "ect1 missing"))

  shock_terms <- intersect(c("pandemic", "pan_step", "pan_win", "brk_mode"), eq_m$term)

  # Normalized cointegrating vector. cajorls fixes beta[1,1] == 1 (the mode
  # column, always first in y from build_system()), so the informative part
  # is beta_ic -- the long-run elasticity of the mode with respect to IC,
  # from the relation  mode_l + beta_ic*IC_l (+ beta_det*det) ~ stationary,
  # i.e.  mode_l = -beta_ic*IC_l (- beta_det*det) in the long run. A third
  # beta row (constant or trend) exists only when ecdet != "none".
  beta <- fit$beta
  beta_ic  <- unname(beta[2, 1])
  beta_det <- if (nrow(beta) >= 3) unname(beta[3, 1]) else NA_real_
  det_label <- if (nrow(beta) >= 3) rownames(beta)[3] else NA_character_

  # Formal weak-exogeneity tests (Waller's request: alrtest, not t-stats on
  # alpha). A is a (p x (p-1)) restriction matrix with alpha = A %*% psi;
  # putting a 1 in row j and zeros elsewhere forces every OTHER row's alpha
  # to zero while leaving row j free, i.e. it tests "every variable except j
  # is weakly exogenous." In a bivariate system that reduces to a single
  # scalar restriction per test: A = (1,0)' tests IC weakly exogenous (mode
  # free); A = (0,1)' tests mode weakly exogenous (IC free). df = r*(p-1).
  we_test <- function(which_free) {
    A <- matrix(if (which_free == 1) c(1, 0) else c(0, 1), nrow = 2)
    t <- tryCatch(urca::alrtest(jo, A = A, r = r), error = function(e) e)
    if (inherits(t, "error")) return(list(LR = NA_real_, p = NA_real_))
    list(LR = unname(t@teststat), p = unname(t@pval[1]))
  }
  we_ic   <- we_test(1)  # restricts IC's alpha to 0, mode left free
  we_mode <- we_test(2)  # restricts mode's alpha to 0, IC left free

  ci_mode <- ci95(eq_m, "ect1"); ci_IC <- ci95(eq_i, "ect1")

  list(ok = TRUE, r = r, status = sprintf("estimated at r = %d", r),
       alpha_mode = grab(eq_m, "ect1", "estimate"),
       t_mode = grab(eq_m, "ect1", "t_stat"),
       p_mode = grab(eq_m, "ect1", "p_value"),
       alpha_mode_lo = unname(ci_mode["lo"]), alpha_mode_hi = unname(ci_mode["hi"]),
       alpha_IC = grab(eq_i, "ect1", "estimate"),
       t_IC = grab(eq_i, "ect1", "t_stat"),
       p_IC = grab(eq_i, "ect1", "p_value"),
       alpha_IC_lo = unname(ci_IC["lo"]), alpha_IC_hi = unname(ci_IC["hi"]),
       beta_ic = beta_ic, beta_det = beta_det, det_label = det_label,
       we_mode_LR = we_mode$LR, we_mode_p = we_mode$p,
       we_IC_LR = we_ic$LR, we_IC_p = we_ic$p,
       shock_term = if (length(shock_terms)) shock_terms[1] else NA_character_,
       shock_coef = if (length(shock_terms)) grab(eq_m, shock_terms[1], "estimate") else NA_real_,
       shock_p = if (length(shock_terms)) grab(eq_m, shock_terms[1], "p_value") else NA_real_,
       crisis_coef = grab(eq_m, "crisis", "estimate"),
       crisis_p = grab(eq_m, "crisis", "p_value"))
}

# ---------------------------------------------------------------------------
# The grid
# ---------------------------------------------------------------------------

run_grid <- function(panel,
                     runs = CFG$runs,
                     dummy_sets = CFG$dummy_sets,
                     ecdet_set = CFG$ecdet_set,
                     lag_rules = CFG$lag_rules,
                     verbose = TRUE) {

  cells <- list(); lags <- list(); joh <- list(); brks <- list(); cover <- list()

  for (run in runs) {
    # A run may specify its own dummy_sets (run E's pre-pandemic window: see
    # config.R) to override the global default -- otherwise every run shares
    # the same list passed in above.
    run_dummy_sets <- if (!is.null(run$dummy_sets)) run$dummy_sets else dummy_sets

    d <- run_sample(panel, run)
    brk_ic <- za_break_date(log(as.numeric(d[[run$ic]])), d$Date)

    cover[[length(cover) + 1]] <- data.frame(
      run = run$tag, label = run$label, requested_start = run$start,
      ic_series = run$ic, n_months = nrow(d),
      actual_start = format(min(d$Date), "%Y-%m"),
      actual_end = format(max(d$Date), "%Y-%m"),
      crisis_months = sum(d$Date >= CFG$crisis_from & d$Date <= CFG$crisis_to),
      ic_break = if (is.na(brk_ic)) NA_character_ else format(brk_ic, "%Y-%m"),
      stringsAsFactors = FALSE)

    if (verbose) {
      cat(sprintf("  %-20s %s..%s  n=%3d  ic_break=%s\n", run$tag,
                  format(min(d$Date), "%Y-%m"), format(max(d$Date), "%Y-%m"), nrow(d),
                  if (is.na(brk_ic)) "none" else format(brk_ic, "%Y-%m")))
    }

    for (mode_name in names(CFG$modes)) {
      y <- build_system(d, mode_name, run)
      mode_col <- colnames(y)[1]; ic_col <- colnames(y)[2]
      brk_mode <- za_break_date(y[, 1], d$Date)

      brks[[length(brks) + 1]] <- data.frame(
        run = run$tag, Model = mode_name, n_months = nrow(d),
        mode_break = if (is.na(brk_mode)) NA_character_ else format(brk_mode, "%Y-%m"),
        ic_series = run$ic,
        ic_break = if (is.na(brk_ic)) NA_character_ else format(brk_ic, "%Y-%m"),
        stringsAsFactors = FALSE)

      for (dset in run_dummy_sets) {
        dum_raw <- build_dummies(dset, d$Date, brk_mode, brk_ic)
        sel <- lag_selection(y, finalize_dummies(dum_raw, 0L)$mat)
        k_by_rule <- vapply(lag_rules, function(r) clamp_k(sel[[r]], CFG$lag_min, CFG$lag_max),
                            integer(1))

        lags[[length(lags) + 1]] <- data.frame(
          run = run$tag, Model = mode_name, dummies = dset,
          AIC_K = sel$AIC, SBC_K = sel$SBC, HQ_K = sel$HQ, FPE_K = sel$FPE,
          stringsAsFactors = FALSE)

        for (ed in ecdet_set) {
          for (rule in lag_rules) {
            k_use <- unname(k_by_rule[rule])
            fd <- finalize_dummies(dum_raw, k_use)
            dum <- fd$mat

            mk <- function(tp) {
              if (is.null(dum)) {
                urca::ca.jo(y, type = tp, ecdet = ed, K = k_use,
                            spec = CFG$jo_spec, season = CFG$season)
              } else {
                urca::ca.jo(y, type = tp, ecdet = ed, K = k_use,
                            spec = CFG$jo_spec, season = CFG$season, dumvar = dum)
              }
            }

            base <- data.frame(
              run = run$tag, Model = mode_name, ic_series = run$ic, dummies = dset,
              dummies_active = if (is.null(dum)) "" else paste(colnames(dum), collapse = ";"),
              dummies_dropped = paste(fd$dropped, collapse = ";"),
              ecdet = ed, lag_rule = rule, K = k_use, n_months = nrow(d),
              cv_valid = (dset == "none"), stringsAsFactors = FALSE)

            jo_tr <- tryCatch(mk("trace"), error = function(e) e)
            jo_ei <- tryCatch(mk("eigen"), error = function(e) e)

            if (inherits(jo_tr, "error") || inherits(jo_ei, "error")) {
              msg <- conditionMessage(if (inherits(jo_tr, "error")) jo_tr else jo_ei)
              cells[[length(cells) + 1]] <- cbind(base, empty_result(paste("ca.jo error:", msg)))
              next
            }

            tr <- rank_from_jo(jo_tr); ei <- rank_from_jo(jo_ei)

            for (tag in c("trace", "eigen")) {
              z <- if (tag == "trace") tr else ei
              for (hyp in c("r = 0", "r <= 1")) {
                st <- if (hyp == "r = 0") z$stat_r0 else z$stat_r1
                cvv <- if (hyp == "r = 0") z$cv_r0 else z$cv_r1
                joh[[length(joh) + 1]] <- cbind(
                  base[, c("run", "Model", "ic_series", "dummies", "ecdet", "lag_rule", "K", "cv_valid")],
                  data.frame(test = tag, hypothesis = hyp, statistic = st,
                             cv_10pct = unname(cvv["10pct"]), cv_5pct = unname(cvv["5pct"]),
                             cv_1pct = unname(cvv["1pct"]),
                             reject_5pct = st > cvv["5pct"], stringsAsFactors = FALSE))
              }
            }

            r_det <- tr[[paste0("rank_", CFG$rank_level)]]
            fit <- fit_vecm(jo_tr, r_det, mode_col, ic_col)
            cells[[length(cells) + 1]] <- cbind(base, result_row(fit, tr, ei, r_det))
          }
        }
      }
    }
  }

  list(cells = do.call(rbind, cells),
       lags = do.call(rbind, lags),
       johansen = do.call(rbind, joh),
       breaks = do.call(rbind, brks),
       coverage = do.call(rbind, cover))
}

empty_result <- function(status) {
  data.frame(rank_10pct = NA_integer_, rank_5pct = NA_integer_, rank_1pct = NA_integer_,
             eigen_rank_5pct = NA_integer_, tests_agree_5pct = NA,
             rank_used = NA_integer_, vecm_status = status,
             alpha_mode = NA_real_, t_mode = NA_real_, p_mode = NA_real_, sig_mode = "",
             alpha_mode_lo = NA_real_, alpha_mode_hi = NA_real_,
             alpha_IC = NA_real_, t_IC = NA_real_, p_IC = NA_real_, sig_IC = "",
             alpha_IC_lo = NA_real_, alpha_IC_hi = NA_real_,
             beta_ic = NA_real_, beta_det = NA_real_, det_label = NA_character_,
             we_mode_LR = NA_real_, we_mode_p = NA_real_,
             we_IC_LR = NA_real_, we_IC_p = NA_real_,
             hl_mode = NA_real_, shock_term = NA_character_,
             shock_coef = NA_real_, shock_p = NA_real_,
             crisis_coef = NA_real_, crisis_p = NA_real_, stringsAsFactors = FALSE)
}

result_row <- function(fit, tr, ei, r_det) {
  ok <- isTRUE(fit$ok)
  data.frame(
    rank_10pct = tr$rank_10pct, rank_5pct = tr$rank_5pct, rank_1pct = tr$rank_1pct,
    eigen_rank_5pct = ei$rank_5pct,
    tests_agree_5pct = tr$rank_5pct == ei$rank_5pct,
    rank_used = if (ok) fit$r else NA_integer_,
    vecm_status = fit$status,
    alpha_mode = if (ok) fit$alpha_mode else NA_real_,
    t_mode = if (ok) fit$t_mode else NA_real_,
    p_mode = if (ok) fit$p_mode else NA_real_,
    sig_mode = if (ok) stars(fit$p_mode) else "",
    alpha_mode_lo = if (ok) fit$alpha_mode_lo else NA_real_,
    alpha_mode_hi = if (ok) fit$alpha_mode_hi else NA_real_,
    alpha_IC = if (ok) fit$alpha_IC else NA_real_,
    t_IC = if (ok) fit$t_IC else NA_real_,
    p_IC = if (ok) fit$p_IC else NA_real_,
    sig_IC = if (ok) stars(fit$p_IC) else "",
    alpha_IC_lo = if (ok) fit$alpha_IC_lo else NA_real_,
    alpha_IC_hi = if (ok) fit$alpha_IC_hi else NA_real_,
    beta_ic = if (ok) fit$beta_ic else NA_real_,
    beta_det = if (ok) fit$beta_det else NA_real_,
    det_label = if (ok) fit$det_label else NA_character_,
    we_mode_LR = if (ok) fit$we_mode_LR else NA_real_,
    we_mode_p = if (ok) fit$we_mode_p else NA_real_,
    we_IC_LR = if (ok) fit$we_IC_LR else NA_real_,
    we_IC_p = if (ok) fit$we_IC_p else NA_real_,
    hl_mode = if (ok) half_life(fit$alpha_mode) else NA_real_,
    shock_term = if (ok) fit$shock_term else NA_character_,
    shock_coef = if (ok) fit$shock_coef else NA_real_,
    shock_p = if (ok) fit$shock_p else NA_real_,
    crisis_coef = if (ok) fit$crisis_coef else NA_real_,
    crisis_p = if (ok) fit$crisis_p else NA_real_,
    stringsAsFactors = FALSE)
}
