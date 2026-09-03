# =============================================================================
# 03_grid.R -- the full specification grid
# =============================================================================
# Every cointegration and VECM result in the project comes from this one
# call. The output is a single long table; later stages slice it rather than
# re-estimating anything.
#
# Dimensions: runs x modes x ecdet x lag rules x dummy designs.
# =============================================================================

section("03  SPECIFICATION GRID")

panel <- load_panel()

# Most runs share CFG$dummy_sets, but a run may specify its own (run E's
# pre-pandemic window uses a reduced set -- see config.R) so this sums
# per-run rather than assuming one uniform count across all runs.
per_cell <- length(CFG$modes) * length(CFG$ecdet_set) * length(CFG$lag_rules)
n_expected <- sum(vapply(CFG$runs, function(r) {
  ds <- if (!is.null(r$dummy_sets)) r$dummy_sets else CFG$dummy_sets
  per_cell * length(ds)
}, numeric(1)))
cat("Cells to estimate:", n_expected, "\n\n")

GRID <- run_grid(panel)

stopifnot(nrow(GRID$cells) == n_expected)

write_out(GRID$cells,    "03_grid", "cells.csv")
write_out(GRID$lags,     "03_grid", "lag_selection.csv")
write_out(GRID$johansen, "03_grid", "johansen_full.csv")
write_out(GRID$breaks,   "03_grid", "break_dates.csv")
write_out(GRID$coverage, "03_grid", "coverage.csv")

cat("\nEndogenous break dates, estimated on each run's own window\n")
print(GRID$breaks, row.names = FALSE)

cat("\nOutcome by run and dummy design\n")
agg <- aggregate(
  cbind(identified = !is.na(GRID$cells$rank_used),
        mode_sig = !is.na(GRID$cells$p_mode) & GRID$cells$p_mode < CFG$alpha_sig,
        IC_sig = !is.na(GRID$cells$p_IC) & GRID$cells$p_IC < CFG$alpha_sig,
        shock_sig = !is.na(GRID$cells$shock_p) & GRID$cells$shock_p < 0.05,
        crisis_sig = !is.na(GRID$cells$crisis_p) & GRID$cells$crisis_p < 0.05) ~
    run + dummies,
  data = GRID$cells, FUN = sum)
dummy_order <- union(CFG$dummy_sets, unique(GRID$cells$dummies))
agg <- agg[order(agg$run, match(agg$dummies, dummy_order)), ]
agg$cells <- per_cell  # constant per (run, design): modes x ecdet x lag rules
print(agg, row.names = FALSE)
write_out(agg, "03_grid", "outcome_by_design.csv")

# Anything that failed outright rather than returning a rank.
failed <- GRID$cells[grepl("error", GRID$cells$vecm_status), ]
if (nrow(failed)) {
  cat("\nCells that errored (", nrow(failed), ")\n", sep = "")
  print(failed[, c("run", "Model", "dummies", "ecdet", "lag_rule", "K", "vecm_status")],
        row.names = FALSE)
} else {
  cat("\nAll cells resolved without error.\n")
}

# Dummy columns that could not be identified at a given lag.
drops <- GRID$cells[nzchar(GRID$cells$dummies_dropped), ]
if (nrow(drops)) {
  cat("\nDummy columns dropped as unidentifiable\n")
  print(as.data.frame(table(design = drops$dummies, dropped = drops$dummies_dropped)),
        row.names = FALSE)
}
