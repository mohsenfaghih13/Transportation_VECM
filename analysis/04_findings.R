# =============================================================================
# 04_findings.R -- what survives, and what only held under one specification
# =============================================================================
# Pure reporting. Reads outputs/03_grid/cells.csv and lag_selection.csv;
# estimates nothing.
#
# A finding counts as robust when it holds across the specification grid, not
# when it holds in the one cell that was originally hardcoded. The headline
# table is therefore built from a survival count, not from a single model.
# =============================================================================

section("04  FINDINGS")

cells <- read.csv(file.path(CFG$out_dir, "03_grid", "cells.csv"),
                  stringsAsFactors = FALSE)
lags <- read.csv(file.path(CFG$out_dir, "03_grid", "lag_selection.csv"),
                 stringsAsFactors = FALSE)

primary_tag <- Filter(function(r) isTRUE(r$primary), CFG$runs)[[1]]$tag
prim <- cells[cells$run == primary_tag, ]

# --- 1. sensitivity of the headline count to the dummy design --------------
cat("Significant mode-side adjustment in the primary system, by dummy design\n")
by_design <- do.call(rbind, lapply(CFG$dummy_sets, function(ds) {
  s <- prim[prim$dummies == ds, ]
  hit <- s[!is.na(s$p_mode) & s$p_mode < CFG$alpha_sig, ]
  data.frame(design = ds, cells = nrow(s),
             identified = sum(!is.na(s$rank_used)),
             mode_sig = nrow(hit),
             models = paste(sort(unique(hit$Model)), collapse = ", "),
             stringsAsFactors = FALSE)
}))
print(by_design, row.names = FALSE)
write_out(by_design, "04_findings", "sensitivity_to_dummy_design.csv")

# --- 2. survival count -----------------------------------------------------
# For each mode, the share of identified specifications showing significant
# negative adjustment. Reported per run: pooling the primary Census IC system
# with the warehouse-construction robustness runs would let a result carried
# by a different inventory measure inflate the primary count.
survival_for <- function(df) {
  do.call(rbind, lapply(names(CFG$modes), function(m) {
    s <- df[df$Model == m & !is.na(df$rank_used), ]
    sig <- s[!is.na(s$p_mode) & s$p_mode < CFG$alpha_sig & s$alpha_mode < 0, ]
    data.frame(
      Model = m, identified = nrow(s), sig_negative = nrow(sig),
      share = if (nrow(s)) round(nrow(sig) / nrow(s), 3) else NA_real_,
      alpha_min = if (nrow(sig)) round(min(sig$alpha_mode), 4) else NA_real_,
      alpha_max = if (nrow(sig)) round(max(sig$alpha_mode), 4) else NA_real_,
      hl_min = if (nrow(sig)) round(min(sig$hl_mode, na.rm = TRUE), 1) else NA_real_,
      hl_max = if (nrow(sig)) round(max(sig$hl_mode, na.rm = TRUE), 1) else NA_real_,
      stringsAsFactors = FALSE)
  }))
}

cat("\nSurvival by run: share of identified cells with significant negative alpha\n")
surv <- do.call(rbind, lapply(CFG$runs, function(r) {
  s <- survival_for(cells[cells$run == r$tag, ])
  cbind(run = r$tag, primary = isTRUE(r$primary), s)
}))
surv <- surv[order(surv$run, -surv$share), ]
print(surv, row.names = FALSE)
write_out(surv, "04_findings", "survival_by_mode.csv")

# Excluding the design known to be mistimed, within the primary system only.
cat("\nSurvival in the primary system, excluding the mistimed 'orig' design\n")
sp <- survival_for(prim[prim$dummies != "orig", ])
print(sp[order(-sp$share), ], row.names = FALSE)
write_out(sp, "04_findings", "survival_primary_excl_orig.csv")

# --- 3. headline table under the preferred design --------------------------
cat("\nHeadline: primary system under dummy design '", CFG$headline_dummy, "'\n", sep = "")
hd <- prim[prim$dummies == CFG$headline_dummy, ]
hd <- hd[order(match(hd$Model, names(CFG$modes)), hd$ecdet, hd$lag_rule), ]
print(data.frame(
  Model = hd$Model, ecdet = hd$ecdet, lag = hd$lag_rule, K = hd$K,
  rank = hd$rank_5pct,
  alpha_mode = fmt_star(hd$alpha_mode, hd$p_mode),
  alpha_IC = fmt_star(hd$alpha_IC, hd$p_IC),
  beta_ic = fmt_num(hd$beta_ic, 3),
  half_life = fmt_num(hd$hl_mode, 1),
  status = hd$vecm_status, stringsAsFactors = FALSE), row.names = FALSE)
write_out(hd, "04_findings", "headline_primary.csv")

# beta_ic is the long-run elasticity of the mode w.r.t. IC, normalized so
# the mode's own coefficient is 1: mode_l = -beta_ic * IC_l (- beta_det *
# det). Reported here, alongside alpha, so a reader can verify a
# significant alpha actually points back toward a sensibly-signed long-run
# relationship rather than just being statistically significant in
# isolation (Waller's 2026-08-26 request).
cat("\nNormalized cointegrating vector, headline design (mode_l = -beta_ic * IC_l)\n")
print(data.frame(Model = hd$Model, ecdet = hd$ecdet, lag = hd$lag_rule,
                  beta_ic = fmt_num(hd$beta_ic, 4),
                  implied_mode_on_IC = fmt_num(-hd$beta_ic, 4),
                  det_label = hd$det_label, beta_det = fmt_num(hd$beta_det, 4),
                  stringsAsFactors = FALSE), row.names = FALSE)

# 95% CI on alpha, from the equation's own residual df via qt(), not a
# normal approximation. Waller's request: treat half-lives as descriptive,
# not as precise ratios comparable across models, until CIs are in view.
cat("\nAlpha with 95% confidence intervals, headline design\n")
print(data.frame(Model = hd$Model, ecdet = hd$ecdet, lag = hd$lag_rule,
                  alpha_mode = fmt_num(hd$alpha_mode, 4),
                  mode_CI = paste0("[", fmt_num(hd$alpha_mode_lo, 4), ", ",
                                    fmt_num(hd$alpha_mode_hi, 4), "]"),
                  alpha_IC = fmt_num(hd$alpha_IC, 4),
                  IC_CI = paste0("[", fmt_num(hd$alpha_IC_lo, 4), ", ",
                                 fmt_num(hd$alpha_IC_hi, 4), "]"),
                  stringsAsFactors = FALSE), row.names = FALSE)

# --- 4. which mode adjusts, and which side ---------------------------------
cat("\nWho adjusts (informal, t-stat on alpha), counted over identified cells",
    "in the primary system\n")
who <- do.call(rbind, lapply(names(CFG$modes), function(m) {
  s <- prim[prim$Model == m & !is.na(prim$rank_used), ]
  data.frame(Model = m, identified = nrow(s),
             mode_adjusts = sum(!is.na(s$p_mode) & s$p_mode < CFG$alpha_sig),
             IC_adjusts = sum(!is.na(s$p_IC) & s$p_IC < CFG$alpha_sig),
             both = sum(!is.na(s$p_mode) & s$p_mode < CFG$alpha_sig &
                        !is.na(s$p_IC) & s$p_IC < CFG$alpha_sig),
             neither = sum((is.na(s$p_mode) | s$p_mode >= CFG$alpha_sig) &
                           (is.na(s$p_IC) | s$p_IC >= CFG$alpha_sig)),
             stringsAsFactors = FALSE)
}))
print(who, row.names = FALSE)
write_out(who, "04_findings", "who_adjusts.csv")

# Formal version of the same question: reject weak exogeneity via alrtest
# (a likelihood-ratio test comparing the restricted-alpha VECM to the
# unrestricted one) rather than reading significance off an individual
# alpha's t-statistic in the unrestricted fit. Waller's 2026-08-26 request.
# The two tables answer a related but distinct question each: the informal
# one asks "is this alpha far from zero"; alrtest asks "does forcing this
# alpha to zero measurably hurt the fit of the whole system." They can and
# sometimes do disagree on individual cells.
cat("\nWho adjusts (formal, alrtest weak-exogeneity LR test), counted over",
    "identified cells in the primary system\n")
who_formal <- do.call(rbind, lapply(names(CFG$modes), function(m) {
  s <- prim[prim$Model == m & !is.na(prim$rank_used), ]
  data.frame(Model = m, identified = nrow(s),
             mode_rejects_exogeneity = sum(!is.na(s$we_mode_p) & s$we_mode_p < CFG$alpha_sig),
             IC_rejects_exogeneity = sum(!is.na(s$we_IC_p) & s$we_IC_p < CFG$alpha_sig),
             both = sum(!is.na(s$we_mode_p) & s$we_mode_p < CFG$alpha_sig &
                        !is.na(s$we_IC_p) & s$we_IC_p < CFG$alpha_sig),
             neither = sum((is.na(s$we_mode_p) | s$we_mode_p >= CFG$alpha_sig) &
                           (is.na(s$we_IC_p) | s$we_IC_p >= CFG$alpha_sig)),
             stringsAsFactors = FALSE)
}))
print(who_formal, row.names = FALSE)
write_out(who_formal, "04_findings", "who_adjusts_formal.csv")

# Where the two methods land on opposite sides of the significance
# threshold for the same cell -- worth knowing which cells these are before
# trusting either table blindly.
disagree_mode <- prim[!is.na(prim$rank_used) &
                       (prim$p_mode < CFG$alpha_sig) != (prim$we_mode_p < CFG$alpha_sig), ]
disagree_IC <- prim[!is.na(prim$rank_used) &
                     (prim$p_IC < CFG$alpha_sig) != (prim$we_IC_p < CFG$alpha_sig), ]
n_disagree <- nrow(disagree_mode) + nrow(disagree_IC)
cat("\nInformal vs. formal test disagree on significance in", n_disagree,
    "of", 2 * sum(!is.na(prim$rank_used)), "(mode, IC) pairs across",
    "identified primary cells.\n")

# --- 5. shock coefficients -------------------------------------------------
cat("\nShock coefficient in the mode equation, by design (identified cells)\n")
sh <- cells[!is.na(cells$shock_p), ]
if (nrow(sh)) {
  shs <- do.call(rbind, lapply(split(sh, list(sh$run, sh$dummies), drop = TRUE), function(g) {
    data.frame(run = g$run[1], design = g$dummies[1], n = nrow(g),
               median_coef = round(median(g$shock_coef, na.rm = TRUE), 5),
               median_p = round(median(g$shock_p, na.rm = TRUE), 4),
               sig_5pct = sum(g$shock_p < 0.05), stringsAsFactors = FALSE)
  }))
  shs <- shs[order(shs$run, match(shs$design, union(CFG$dummy_sets, unique(cells$dummies)))), ]
  print(shs, row.names = FALSE)
  write_out(shs, "04_findings", "shock_coefficients.csv")
}

crisis_hits <- sum(!is.na(cells$crisis_p) & cells$crisis_p < 0.05)
cat("\nCrisis dummy significant at 5% in", crisis_hits, "of",
    sum(!is.na(cells$crisis_p)), "cells where it is identified.\n")

# --- 6. do HQ and FPE agree with AIC/SBC on lag order? ----------------------
# Sparsh's 2026-08-26 email: AIC and SBC disagree on every model (6 lags vs
# 2, typically); check whether HQ and FPE side with one or the other, as a
# corroborating signal. This does not change which lag is used anywhere --
# AIC and SBC remain CFG$lag_rules -- it is purely informational.
lags$HQ_matches_AIC  <- lags$HQ_K == lags$AIC_K
lags$HQ_matches_SBC  <- lags$HQ_K == lags$SBC_K
lags$FPE_matches_AIC <- lags$FPE_K == lags$AIC_K
lags$FPE_matches_SBC <- lags$FPE_K == lags$SBC_K

cat("\nDo HQ and FPE agree with AIC or SBC on lag order? (all runs, all",
    "dummy designs, one row per mode x design)\n")
agree <- data.frame(
  n = nrow(lags),
  `HQ==AIC` = sum(lags$HQ_matches_AIC), `HQ==SBC` = sum(lags$HQ_matches_SBC),
  `FPE==AIC` = sum(lags$FPE_matches_AIC), `FPE==SBC` = sum(lags$FPE_matches_SBC),
  check.names = FALSE)
print(agree, row.names = FALSE)
write_out(lags, "04_findings", "lag_criteria_agreement.csv")

if (all(!lags$HQ_matches_AIC & !lags$HQ_matches_SBC) &&
    all(!lags$FPE_matches_AIC & !lags$FPE_matches_SBC)) {
  cat("Neither HQ nor FPE agrees with AIC or SBC in any cell -- both sit",
      "at a lag order distinct from either rule used in the grid.\n")
}

# --- 7. caveats that must travel with these numbers ------------------------
cat("\nCaveats\n")
cat("  * Johansen critical values are not adjusted for dumvar. Rows with\n")
cat("    cv_valid = FALSE (", sum(!cells$cv_valid), " of ", nrow(cells),
    ") have indicative ranks only.\n", sep = "")
cat("  * A bivariate system admits r = 1 only, so r = 0 and r = 2 cells are\n")
cat("    reported unestimated rather than forced.\n")
cat("  * Half-lives are defined only for -1 < alpha < 0.\n")
