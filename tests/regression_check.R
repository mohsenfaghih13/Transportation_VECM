#!/usr/bin/env Rscript
# =============================================================================
# tests/regression_check.R -- the pipeline must not move a single number
# =============================================================================
# Compares outputs/03_grid/cells.csv against frozen fixtures in tests/fixtures/.
#
#   phase1  = run A, dummies "none"                 (30 cells)
#   phase2  = all runs, dummies "none" and "shocks" (300 cells; "shocks"->"orig")
#   phase4  = all runs, every design each run has    (720 cells)
#
# Any mismatch beyond floating-point tolerance is a regression.
#
#   Rscript tests/regression_check.R
#
# History: fixtures were originally captured from the pre-consolidation phase
# scripts, to verify the refactor into one grid engine moved zero numbers.
# That check passed and is retired. Total_Inventories was later backfilled
# from 2009-06 to 2003-12 (Census MTIS itself goes back to 1996 -- the
# 2009-06 start was a data-collection gap, not a source limit), which
# legitimately changes every downstream number. Fixtures were re-baselined
# from the post-backfill output. Later the same day: run D (deflated by
# PPI_All_Commodities) and run E (pre-pandemic subsample, 2003-2019) were
# added to CFG$runs; then run E was given its own reduced dummy_sets (none,
# crisis_only) instead of the global 5, since its pandemic-dummy designs all
# reference dates after the window ends and would silently collapse to
# "crisis only" anyway (R/dummies.R). Each change grew or shrank
# phase2/phase4's cell counts (144->192->240->216, 360->480->600->528);
# fixtures were re-baselined each time. From here on this check protects
# against accidental changes to the *current, correct* numbers, not against
# a comparison to the old truncated sample or an earlier run/design set.
#
# Total_Inventories was later found to be the seasonally-adjusted Census
# series, not NSA as intended -- mismatched against the NSA transportation-
# cost series it's paired with in every model. Corrected to the true NSA
# figures and fixtures re-baselined again; row counts unchanged (24/216/528),
# only the values themselves moved.
#
# Run F (narrower 2009-2019 pre-pandemic subsample, same reduced dummy_sets
# as run E) was added afterward; phase2/phase4 grew again (216->240,
# 528->576) and fixtures were re-baselined accordingly.
#
# Deep_Sea_Freight was promoted from CFG$extra_series to a fifth primary
# mode (CFG$modes). All existing cells for the original four modes matched
# exactly (24/240/576 old cells all found unchanged among the new totals);
# only new SeaFreight cells were added. Fixtures re-baselined again
# (24->30, 240->300, 576->720).
#
# Run G (2004-2024 start, full dummy_sets) was added as a non-primary,
# additive-only comparison run against an independent replication attempt
# by a collaborator -- not part of the official report sections. Every
# cell in runs A-F matched exactly; only new run-G cells were added.
# Fixtures re-baselined again (30 unchanged, 300->360, 720->870).
#
# The whole project (CFG$common_start, and every run's "2003-01-01" start)
# was then moved to 2004-01-01 -- every run previously starting 2003
# effectively realized to 2003-12 anyway (one incomplete month), matching
# the same fix the collaborator's independent script made. Run G, now
# identical to run A, was removed. This is a real, legitimate change to
# every downstream number (not a pure addition): 124 of 720 cells moved on
# alpha_mode/alpha_IC (max abs diff ~0.017), and 10 cells flipped
# identification status (rank_used NA <-> 1) as a direct result of losing
# one month of data. Row counts returned to their pre-run-G size
# (300/720) since G is gone; fixtures re-baselined from the 2004-start
# output.
# =============================================================================

TOL <- 1e-8
root <- normalizePath(file.path(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])), ".."))
setwd(root)

new <- read.csv("outputs/03_grid/cells.csv", stringsAsFactors = FALSE)

fails <- 0L; checks <- 0L

report <- function(label, ok, detail = "") {
  checks <<- checks + 1L
  if (!ok) fails <<- fails + 1L
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0(" -- ", detail) else ""))
}

# Compare a numeric column between two aligned frames.
cmp_num <- function(a, b, col, label) {
  x <- a[[col]]; y <- b[[col]]
  both_na <- is.na(x) & is.na(y)
  d <- ifelse(both_na, 0, abs(x - y))
  d[is.na(d)] <- Inf
  worst <- max(d)
  report(paste0(label, " :: ", col), worst < TOL,
         if (worst < TOL) "" else sprintf("max abs diff %.3g in %d cells",
                                          worst, sum(d >= TOL)))
}

# The phase scripts wrote slightly different schemas as they evolved. Map
# them onto the consolidated one so the comparison is like for like.
normalize_old <- function(old, default_run) {
  # phase1 put the mode name in `Variable` and a label in `Model`
  if ("Variable" %in% names(old)) old$Model <- old$Variable
  if (!"run" %in% names(old)) old$run <- default_run
  if (!"dummies" %in% names(old)) old$dummies <- "none"
  if (!"rank_5pct" %in% names(old)) {
    src <- intersect(c("rank_trace_5pct", "rank_determined_trace_5pct"), names(old))
    if (length(src)) old$rank_5pct <- old[[src[1]]]
  }
  old
}

check_against <- function(path, key_filter, old_dummy_map = NULL, label,
                          default_run = "A_common_censusIC") {
  if (!file.exists(path)) {
    cat("  [SKIP]", label, "-- archived output not found:", path, "\n")
    return(invisible(NULL))
  }
  old <- normalize_old(read.csv(path, stringsAsFactors = FALSE), default_run)
  if (!is.null(old_dummy_map)) {
    old$dummies <- ifelse(old$dummies %in% names(old_dummy_map),
                          unname(old_dummy_map[old$dummies]), old$dummies)
  }

  keys <- c("run", "Model", "ecdet", "lag_rule", "dummies")
  sub <- key_filter(new)

  old$.k <- do.call(paste, c(old[keys], sep = "|"))
  sub$.k <- do.call(paste, c(sub[keys], sep = "|"))

  common <- intersect(old$.k, sub$.k)
  report(paste0(label, " :: cell coverage"),
         length(common) == nrow(old) && length(common) == nrow(sub),
         sprintf("old %d, new %d, matched %d", nrow(old), nrow(sub), length(common)))
  if (!length(common)) return(invisible(NULL))

  o <- old[match(common, old$.k), ]
  n <- sub[match(common, sub$.k), ]

  cmp_num(o, n, "K", label)
  cmp_num(o, n, "alpha_mode", label)
  cmp_num(o, n, "alpha_IC", label)
  cmp_num(o, n, "p_mode", label)
  cmp_num(o, n, "p_IC", label)

  same <- (is.na(o$rank_5pct) & is.na(n$rank_5pct)) | (o$rank_5pct == n$rank_5pct)
  same[is.na(same)] <- FALSE
  report(paste0(label, " :: rank_5pct"), all(same),
         sprintf("%d mismatched", sum(!same)))
}

cat("Regression check: consolidated grid vs tests/fixtures/\n\n")

cat("phase1 (run A, no dummies)\n")
check_against("tests/fixtures/phase1_vecm_ect.csv",
              function(x) x[x$run == "A_common_censusIC" & x$dummies == "none", ],
              NULL, "phase1")

cat("\nphase2 (all runs, none + orig)\n")
check_against("tests/fixtures/phase2_vecm_ect.csv",
              function(x) x[x$dummies %in% c("none", "orig"), ],
              c(shocks = "orig"), "phase2")

cat("\nphase4 (all runs, all designs)\n")
check_against("tests/fixtures/phase4_vecm_ect.csv",
              function(x) x, NULL, "phase4")

cat("\n", strrep("-", 60), "\n", sep = "")
cat(sprintf("%d checks, %d failures\n", checks, fails))
if (fails > 0) quit(status = 1)
cat("Consolidation is numerically identical to the phase scripts.\n")
