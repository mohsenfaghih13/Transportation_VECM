# =============================================================================
# R/grid_combined.R -- specification grid for the 3-variable combined systems
# =============================================================================
# Generalizes the bivariate engine in grid.R to an arbitrary number of
# variables p. Two structural differences from the pairwise grid:
#   * rank can be 0..p-1, not just 0/1, so it is read from the trace-test
#     rows generically instead of the hardcoded "r = 0" / "r <= 1" lookup
#     rank_from_jo() uses.
#   * a rank r >= 2 produces r cointegrating relations, so an equation can
#     carry more than one ect term. Loadings are reported long -- one row
#     per (equation, ect term) -- rather than forced into fixed
#     alpha_mode/alpha_IC columns the way the bivariate grid does.
# =============================================================================

# Rank at every configured critical-value level, for a Johansen object with
# an arbitrary number of variables. rownames(jo@cval) are "r = 0  |",
# "r <= 1 |", "r <= 2 |", ... in decreasing-r print order; grep finds each
# by pattern rather than assuming that order.
rank_from_jo_n <- function(jo) {
  stat <- as.numeric(jo@teststat); cv <- jo@cval; rn <- rownames(cv)
  p <- ncol(jo@x)
  idx <- vapply(0:(p - 1), function(i) {
    pat <- if (i == 0) "r = 0" else sprintf("r <= %d", i)
    grep(pat, rn, fixed = TRUE)[1]
  }, integer(1))
  rank_at <- function(l) {
    for (k in seq_along(idx)) if (!(stat[idx[k]] > cv[idx[k], l])) return(k - 1L)
    p
  }
  out <- list()
  for (l in CFG$cv_levels) out[[paste0("rank_", l)]] <- rank_at(l)
  out
}

# Fit at the rank the trace test gave, for an arbitrary number of variables.
# r = 0 (no cointegration) and r = p (full rank, system already stationary)
# are recorded with their reason rather than estimated. Loadings come back
# long: one row per (equation variable, ect term) -- ect2 etc. exist only
# when r >= 2.
fit_vecm_n <- function(jo, r, var_names) {
  p <- length(var_names)
  if (is.na(r)) {
    return(list(ok = FALSE, r = NA_integer_, status = "rank undetermined", loadings = NULL))
  }
  if (r < 1) return(list(ok = FALSE, r = r, status = "r = 0: no cointegration", loadings = NULL))
  if (r >= p) return(list(ok = FALSE, r = r, status = sprintf("r = %d: full rank", r), loadings = NULL))

  fit <- tryCatch(urca::cajorls(jo, r = r), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, r = r, status = paste("cajorls error:", conditionMessage(fit)),
                loadings = NULL))
  }

  sm <- summary(fit$rlm)
  rows <- lapply(var_names, function(v) {
    key <- paste("Response", paste0(v, ".d"))
    if (!key %in% names(sm)) return(NULL)
    cf <- sm[[key]]$coefficients
    ect_rows <- grep("^ect[0-9]+$", rownames(cf))
    if (!length(ect_rows)) return(NULL)
    data.frame(equation = v, ect_term = rownames(cf)[ect_rows],
               estimate = unname(cf[ect_rows, "Estimate"]),
               t_stat = unname(cf[ect_rows, "t value"]),
               p_value = unname(cf[ect_rows, "Pr(>|t|)"]),
               stringsAsFactors = FALSE)
  })
  loadings <- do.call(rbind, rows)

  # Normalized cointegrating vector(s), long format: one row per (row of
  # beta -- a variable or a deterministic term, ect_term). cajorls
  # normalizes each ect column so a different one of the first r variables
  # equals 1 in that column (Phillips normalization), so unlike the
  # bivariate grid's single beta_ic there is no one "the interesting
  # coefficient" here -- report the whole matrix and let the reader see
  # which variable each relation is normalized on.
  bmat <- fit$beta
  beta_long <- do.call(rbind, lapply(seq_len(ncol(bmat)), function(j) {
    data.frame(ect_term = colnames(bmat)[j], row = rownames(bmat),
               coefficient = unname(bmat[, j]), stringsAsFactors = FALSE)
  }))

  # Formal weak-exogeneity test per variable (Waller's request, generalized
  # from the bivariate grid's we_mode/we_IC to p variables). To test whether
  # variable i is weakly exogenous, restrict alpha = A %*% psi with A the
  # (p x (p-1)) identity matrix with row i removed: every OTHER variable's
  # alpha stays free, row i is forced to zero. df = r*(p-(p-1)) = r.
  we_long <- do.call(rbind, lapply(seq_along(var_names), function(i) {
    A <- diag(p)[, -i, drop = FALSE]
    t <- tryCatch(urca::alrtest(jo, A = A, r = r), error = function(e) e)
    if (inherits(t, "error")) {
      return(data.frame(variable = var_names[i], LR = NA_real_, df = r, p_value = NA_real_,
                         stringsAsFactors = FALSE))
    }
    data.frame(variable = var_names[i], LR = unname(t@teststat), df = r,
               p_value = unname(t@pval[1]), stringsAsFactors = FALSE)
  }))

  list(ok = TRUE, r = r, status = sprintf("estimated at r = %d", r),
       loadings = loadings, beta = beta_long, weak_exo = we_long)
}

run_combined_grid <- function(panel,
                              systems = CFG$combined_systems,
                              dummy_sets = CFG$combined_dummy_sets,
                              ecdet_set = CFG$ecdet_set,
                              lag_rules = CFG$lag_rules,
                              verbose = TRUE) {

  cells <- list(); loadings <- list(); betas <- list(); weak_exo <- list(); lags <- list()

  for (sysdef in systems) {
    d <- run_sample(panel, sysdef, modes = sysdef$vars)
    var_names <- c(names(sysdef$vars), "IC")
    y <- build_system_n(d, c(sysdef$vars, IC = sysdef$ic))

    if (verbose) {
      cat(sprintf("  %-14s %s..%s  n=%3d  vars=%s\n", sysdef$tag,
                  format(min(d$Date), "%Y-%m"), format(max(d$Date), "%Y-%m"),
                  nrow(d), paste(var_names, collapse = "+")))
    }

    for (dset in dummy_sets) {
      dum_raw <- build_dummies(dset, d$Date)
      sel <- lag_selection(y, finalize_dummies(dum_raw, 0L)$mat)
      k_by_rule <- vapply(lag_rules, function(r) clamp_k(sel[[r]], CFG$lag_min, CFG$lag_max),
                          integer(1))

      lags[[length(lags) + 1]] <- data.frame(
        system = sysdef$tag, dummies = dset,
        AIC_K = sel$AIC, SBC_K = sel$SBC, HQ_K = sel$HQ, FPE_K = sel$FPE,
        stringsAsFactors = FALSE)

      for (ed in ecdet_set) {
        for (rule in lag_rules) {
          k_use <- unname(k_by_rule[rule])
          fd <- finalize_dummies(dum_raw, k_use)
          dum <- fd$mat

          mk <- function(tp) {
            if (is.null(dum)) {
              urca::ca.jo(y, type = tp, ecdet = ed, K = k_use, spec = CFG$jo_spec, season = CFG$season)
            } else {
              urca::ca.jo(y, type = tp, ecdet = ed, K = k_use, spec = CFG$jo_spec,
                          season = CFG$season, dumvar = dum)
            }
          }

          base <- data.frame(
            system = sysdef$tag, label = sysdef$label, dummies = dset,
            dummies_active = if (is.null(dum)) "" else paste(colnames(dum), collapse = ";"),
            dummies_dropped = paste(fd$dropped, collapse = ";"),
            ecdet = ed, lag_rule = rule, K = k_use, n_months = nrow(d),
            cv_valid = (dset == "none"), stringsAsFactors = FALSE)

          jo_tr <- tryCatch(mk("trace"), error = function(e) e)
          jo_ei <- tryCatch(mk("eigen"), error = function(e) e)

          if (inherits(jo_tr, "error") || inherits(jo_ei, "error")) {
            msg <- conditionMessage(if (inherits(jo_tr, "error")) jo_tr else jo_ei)
            cells[[length(cells) + 1]] <- cbind(base, data.frame(
              rank_10pct = NA_integer_, rank_5pct = NA_integer_, rank_1pct = NA_integer_,
              eigen_rank_5pct = NA_integer_, tests_agree_5pct = NA,
              rank_used = NA_integer_, vecm_status = paste("ca.jo error:", msg),
              stringsAsFactors = FALSE))
            next
          }

          tr <- rank_from_jo_n(jo_tr); ei <- rank_from_jo_n(jo_ei)
          r_det <- tr[[paste0("rank_", CFG$rank_level)]]
          fit <- fit_vecm_n(jo_tr, r_det, var_names)

          cells[[length(cells) + 1]] <- cbind(base, data.frame(
            rank_10pct = tr$rank_10pct, rank_5pct = tr$rank_5pct, rank_1pct = tr$rank_1pct,
            eigen_rank_5pct = ei$rank_5pct,
            tests_agree_5pct = tr$rank_5pct == ei$rank_5pct,
            rank_used = if (isTRUE(fit$ok)) fit$r else NA_integer_,
            vecm_status = fit$status, stringsAsFactors = FALSE))

          if (isTRUE(fit$ok) && !is.null(fit$loadings)) {
            meta <- base[, c("system", "dummies", "ecdet", "lag_rule", "K", "cv_valid")]
            loadings[[length(loadings) + 1]] <- cbind(meta, fit$loadings)
            betas[[length(betas) + 1]] <- cbind(meta, fit$beta)
            weak_exo[[length(weak_exo) + 1]] <- cbind(meta, fit$weak_exo)
          }
        }
      }
    }
  }

  list(cells = do.call(rbind, cells),
       loadings = do.call(rbind, loadings),
       beta = do.call(rbind, betas),
       weak_exo = do.call(rbind, weak_exo),
       lags = do.call(rbind, lags))
}
