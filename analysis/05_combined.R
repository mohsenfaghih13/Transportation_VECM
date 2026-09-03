# =============================================================================
# 05_combined.R -- three-variable combined systems (Truckload/LTL + Air + IC)
# =============================================================================
# Sparsh's 2026-08-26 request: a system of Truckload + Airfreight + IC, with
# an LTL-substituted alternate -- never both trucking variables together.
# Reads CFG$combined_systems; nothing here estimates outside run_combined_grid().
#
# Rank can be 0, 1 or 2 in this three-variable system (not just 0/1 as in
# the bivariate grid): a rank of 2 means TWO cointegrating relations, so an
# equation carries ect1 AND ect2 loadings, both reported.
# =============================================================================

section("05  COMBINED SYSTEMS")

panel <- load_panel()

n_expected <- length(CFG$combined_systems) * length(CFG$combined_dummy_sets) *
              length(CFG$ecdet_set) * length(CFG$lag_rules)
cat("Cells to estimate:", n_expected, "\n\n")

COMB <- run_combined_grid(panel)

stopifnot(nrow(COMB$cells) == n_expected)

write_out(COMB$cells,    "05_combined", "cells.csv")
write_out(COMB$loadings, "05_combined", "loadings.csv")
write_out(COMB$beta,     "05_combined", "cointegrating_vectors.csv")
write_out(COMB$weak_exo, "05_combined", "weak_exogeneity.csv")
write_out(COMB$lags,     "05_combined", "lag_selection.csv")

cat("\nRank and status by system and dummy design\n")
print(COMB$cells[, c("system", "dummies", "ecdet", "lag_rule", "K",
                      "rank_5pct", "eigen_rank_5pct", "tests_agree_5pct",
                      "rank_used", "vecm_status")],
      row.names = FALSE)

failed <- COMB$cells[grepl("error", COMB$cells$vecm_status), ]
if (nrow(failed)) {
  cat("\nCells that errored (", nrow(failed), ")\n", sep = "")
  print(failed[, c("system", "dummies", "ecdet", "lag_rule", "K", "vecm_status")],
        row.names = FALSE)
} else {
  cat("\nAll cells resolved without error.\n")
}

drops <- COMB$cells[nzchar(COMB$cells$dummies_dropped), ]
if (nrow(drops)) {
  cat("\nDummy columns dropped as unidentifiable\n")
  print(as.data.frame(table(design = drops$dummies, dropped = drops$dummies_dropped)),
        row.names = FALSE)
}

if (!is.null(COMB$loadings) && nrow(COMB$loadings)) {
  cat("\nError-correction loadings, all identified cells (long format)\n")
  ld <- COMB$loadings
  print(data.frame(system = ld$system, dummies = ld$dummies, ecdet = ld$ecdet,
                    lag = ld$lag_rule, K = ld$K, equation = ld$equation,
                    ect_term = ld$ect_term,
                    estimate = fmt_star(ld$estimate, ld$p_value),
                    p = fmt_num(ld$p_value, 4), stringsAsFactors = FALSE),
        row.names = FALSE)

  cat("\nWho adjusts, counted over identified cells with a single cointegrating",
      "relation (r = 1)\n")
  ld1 <- ld[ld$ect_term == "ect1", ]
  who <- do.call(rbind, lapply(split(ld1, list(ld1$system, ld1$equation), drop = TRUE),
                                function(g) {
    data.frame(system = g$system[1], equation = g$equation[1], n = nrow(g),
               sig = sum(g$p_value < CFG$alpha_sig, na.rm = TRUE),
               sig_negative = sum(g$p_value < CFG$alpha_sig & g$estimate < 0, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  who <- who[order(who$system, who$equation), ]
  print(who, row.names = FALSE)
  write_out(who, "05_combined", "who_adjusts_r1.csv")
} else {
  cat("\nNo cells identified a cointegrating rank in [1, p-1].\n")
}

if (!is.null(COMB$weak_exo) && nrow(COMB$weak_exo)) {
  cat("\nFormal weak-exogeneity tests (alrtest), all identified cells\n")
  we <- COMB$weak_exo
  print(data.frame(system = we$system, dummies = we$dummies, ecdet = we$ecdet,
                    lag = we$lag_rule, K = we$K, variable = we$variable,
                    LR = fmt_num(we$LR, 3), df = we$df,
                    p = fmt_num(we$p_value, 4), stringsAsFactors = FALSE),
        row.names = FALSE)

  cat("\nWho adjusts (formal, alrtest), counted over all identified cells",
      "regardless of rank\n")
  who_formal <- do.call(rbind, lapply(split(we, list(we$system, we$variable), drop = TRUE),
                                       function(g) {
    data.frame(system = g$system[1], variable = g$variable[1], n = nrow(g),
               rejects_exogeneity = sum(g$p_value < CFG$alpha_sig, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  who_formal <- who_formal[order(who_formal$system, who_formal$variable), ]
  print(who_formal, row.names = FALSE)
  write_out(who_formal, "05_combined", "who_adjusts_formal.csv")
}

if (!is.null(COMB$beta) && nrow(COMB$beta)) {
  cat("\nNormalized cointegrating vector(s), all identified cells\n")
  cat("(cajorls normalizes each ect column on a different one of the first\n")
  cat(" r variables -- read the row equal to 1.000 in a column as that\n")
  cat(" column's normalizing variable, not as a result.)\n")
  bt <- COMB$beta
  print(data.frame(system = bt$system, dummies = bt$dummies, ecdet = bt$ecdet,
                    lag = bt$lag_rule, K = bt$K, ect_term = bt$ect_term,
                    row = bt$row, coefficient = fmt_num(bt$coefficient, 4),
                    stringsAsFactors = FALSE), row.names = FALSE)
}

cat("\nCaveats\n")
cat("  * za_series dummy design is not run here: it dates a break per mode,\n")
cat("    and this system has two mode-side variables, so \"the\" mode break\n")
cat("    is not well defined without extending build_dummies().\n")
cat("  * Johansen critical values are not adjusted for dumvar. Rows with\n")
cat("    cv_valid = FALSE have indicative ranks only.\n")
cat("  * A rank of 2 here means TWO cointegrating relations; ect1 and ect2\n")
cat("    are both reported for the relevant equations, neither privileged.\n")
