# =============================================================================
# R/dummies.R -- endogenous break dating and shock dummy construction
# =============================================================================

# Zivot-Andrews break date for one series, estimated on the sample passed in
# rather than on the full history. The distinction matters: local trucking
# breaks at 2020-06 on the 187-month window but 2016-01 on the 241-month one.
za_break_date <- function(x, dates, lag = CFG$za_lag) {
  o <- tryCatch(urca::ur.za(x, model = CFG$za_model, lag = lag), error = function(e) e)
  if (inherits(o, "error")) return(as.Date(NA))
  bp <- o@bpoint
  if (is.na(bp) || bp < 1 || bp > length(dates)) return(as.Date(NA))
  dates[bp]
}

# Raw dummy matrix for a design. Degenerate columns are NOT removed here;
# that depends on the lag order and is handled by finalize_dummies().
#
# crisis_only exists for windows that end before the pandemic (run E,
# 2003-2019): orig/za_step/za_window/za_series all use fixed 2020+ dates for
# their pandemic component, which fall entirely outside such a window and so
# get dropped by finalize_dummies() anyway -- they'd all silently collapse to
# "crisis only" under the hood. crisis_only makes that explicit up front
# instead of running four designs that turn out to be identical.
build_dummies <- function(design, dates, brk_mode = as.Date(NA), brk_ic = as.Date(NA)) {
  if (design == "none") return(NULL)
  crisis <- as.numeric(dates >= CFG$crisis_from & dates <= CFG$crisis_to)

  switch(design,
    crisis_only = cbind(crisis = crisis),
    orig = cbind(crisis = crisis,
                 pandemic = as.numeric(dates >= CFG$orig_pan_from & dates <= CFG$orig_pan_to)),
    za_step = cbind(crisis = crisis,
                    pan_step = as.numeric(dates >= CFG$za_pan_from)),
    za_window = cbind(crisis = crisis,
                      pan_win = as.numeric(dates >= CFG$za_pan_from & dates <= CFG$za_win_to)),
    za_series = {
      m <- cbind(crisis = crisis)
      if (!is.na(brk_mode)) m <- cbind(m, brk_mode = as.numeric(dates >= brk_mode))
      if (!is.na(brk_ic))   m <- cbind(m, brk_IC   = as.numeric(dates >= brk_ic))
      m
    },
    stop("Unknown dummy design: ", design)
  )
}

# Drop columns that cannot be estimated at this lag order.
#
# Two failure modes, both of which produce an exactly singular system in
# ca.jo rather than a warning:
#   1. A column that is constant on the EFFECTIVE sample. ca.jo discards the
#      first K observations, so a dummy whose variation lies entirely inside
#      that burn-in becomes all-zero once estimation starts. The crisis dummy
#      on the 187-month window has only 7 non-zero months and hits this at
#      K >= 7.
#   2. Two identical columns, which happens under za_series when a mode and
#      the inventory series share a break date.
finalize_dummies <- function(m, k_drop) {
  if (is.null(m) || !ncol(m)) return(list(mat = NULL, dropped = character(0)))

  eff <- if (k_drop > 0 && k_drop < nrow(m)) m[-seq_len(k_drop), , drop = FALSE] else m
  keep <- apply(eff, 2, function(x) length(unique(x)) > 1)

  if (sum(keep) > 1) {
    sub <- m[, keep, drop = FALSE]
    dup <- duplicated(t(sub))
    keep[which(keep)[dup]] <- FALSE
  }

  list(mat = if (any(keep)) m[, keep, drop = FALSE] else NULL,
       dropped = colnames(m)[!keep])
}
