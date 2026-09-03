# =============================================================================
# config.R -- every analytical choice in the project, in one place
# =============================================================================
# Nothing below this file should hardcode a sample date, a lag bound, a
# deterministic term, a cointegrating rank or a shock window. If you want to
# change what the analysis does, change it here and re-run run_all.R.
# =============================================================================

CFG <- list()

# --- input -----------------------------------------------------------------
CFG$data_file <- "Transportation_Inventory_Complete_10_Modes_WPU30.xlsx"
CFG$out_dir   <- "outputs"

# --- variables -------------------------------------------------------------
# The four modes that carry the primary pairwise systems.
CFG$modes <- c(LTL        = "Trucking_LD_LTL",
               Truckload  = "Trucking_LD_Truckload",
               Local      = "Trucking_Local",
               Airfreight = "Airfreight_Scheduled")

# Additional series carried through the univariate integration tests only.
CFG$extra_series <- c("Airfreight_NonScheduled", "Rail_Transportation",
                      "Inland_Water_Freight", "Deep_Sea_Freight",
                      "TC_Aggregate_WPU30", "Total_Inventories",
                      "Warehouse_Construction", "Warehousing_Storage",
                      "PPI_All_Commodities")

# Deflator for run D. Dividing an I(1) nominal series by an I(1) deflator
# does not automatically produce another I(1) series -- it does only if the
# nominal series and the deflator are themselves 1:1 cointegrated. Worth
# checking PPI_All_Commodities's own I(1) status (it's in extra_series
# above precisely so 02_integration tests it) before trusting run D's
# order-of-integration verdicts at face value.
CFG$deflator_series <- "PPI_All_Commodities"

# --- estimation runs -------------------------------------------------------
# Each run pairs an inventory-side series with a requested start date. The
# realised sample is whatever complete-case filtering allows, which is not
# the same thing: Warehouse_Construction begins 2004-12, so run C resolves
# to fewer months than requested regardless of how early the start is set.
#
#   A  the primary system: mode + Census MTIS inventories
#   B  same window, different inventory measure  -> isolates the measure
#   C  same measure, longer window               -> isolates the window
#   D  same as A, deflated by an overall PPI      -> guards against a shared
#      inflation trend manufacturing spurious cointegration (Waller,
#      2026-08-26; deflator confirmed by Sparsh + Waller, 2026-09-02)
#   E  same as A, ends 2019-12                    -> pre-pandemic subsample;
#      tests whether the mode-level pattern survives excluding 2020 on
#      (Waller, 2026-08-26; scope confirmed lighter-weight by Sparsh,
#      2026-09-02). Uses its own reduced dummy_sets (below): orig/za_step/
#      za_window/za_series all reference fixed 2020+ dates that fall outside
#      this window and would silently collapse to "crisis only" anyway (see
#      R/dummies.R), so run E is given crisis_only explicitly instead of
#      running four designs that turn out to be identical.
CFG$runs <- list(
  list(tag = "A_common_censusIC",   start = "2003-01-01", ic = "Total_Inventories",
       label = "Census MTIS IC, common window", primary = TRUE),
  list(tag = "B_common_whseconstr", start = "2009-06-01", ic = "Warehouse_Construction",
       label = "Warehouse construction, common window", primary = FALSE),
  list(tag = "C_full_whseconstr",   start = "2003-01-01", ic = "Warehouse_Construction",
       label = "Warehouse construction, full window", primary = FALSE),
  list(tag = "D_deflated_censusIC", start = "2003-01-01", ic = "Total_Inventories",
       label = "Census MTIS IC, deflated by All-Commodities PPI",
       primary = FALSE, deflate = TRUE),
  list(tag = "E_prepandemic_censusIC", start = "2003-01-01", end = "2019-12-31",
       ic = "Total_Inventories", label = "Pre-pandemic subsample (2003-2019), Census MTIS IC",
       primary = FALSE, dummy_sets = c("none", "crisis_only"))
)

# --- specification grid ----------------------------------------------------
# The rank is never fixed. It is read from the trace test at CFG$rank_level
# and the VECM is estimated only where that rank is estimable.
CFG$ecdet_set  <- c("none", "const", "trend")
CFG$lag_rules  <- c("AIC", "SBC")
CFG$lag_max    <- 12L
CFG$lag_min    <- 2L          # ca.jo requires K >= 2
CFG$season     <- 12L
CFG$jo_spec    <- "longrun"
CFG$rank_level <- "5pct"      # level at which the reported rank is read
CFG$cv_levels  <- c("10pct", "5pct", "1pct")

# --- combined systems --------------------------------------------------------
# Three-variable systems, per Sparsh's 2026-08-26 request: Truckload +
# Airfreight + Census IC as the primary combined model, with an LTL-
# substituted alternate. Never both trucking variables in the same system --
# TL dominates truck spend and is primary, LTL is supplemental.
CFG$combined_systems <- list(
  list(tag = "TL_Air_IC",  vars = c(Truckload  = "Trucking_LD_Truckload",
                                     Airfreight = "Airfreight_Scheduled"),
       ic = "Total_Inventories", start = "2003-01-01",
       label = "Truckload + Airfreight + Census IC", primary = TRUE),
  list(tag = "LTL_Air_IC", vars = c(LTL        = "Trucking_LD_LTL",
                                     Airfreight = "Airfreight_Scheduled"),
       ic = "Total_Inventories", start = "2003-01-01",
       label = "LTL + Airfreight + Census IC (alternate)", primary = FALSE)
)

# za_series is excluded here: it dates a break per mode, and a combined
# system has two mode-side variables, so "the" mode break is not well
# defined without extending build_dummies(). The other four designs use
# fixed, series-independent dates and are unaffected.
CFG$combined_dummy_sets <- c("none", "orig", "za_step", "za_window")

# --- shock dummies ---------------------------------------------------------
# Phase 3 dated the breaks endogenously and found the shift sits in
# mid-to-late 2021, not 2020, and that nothing breaks in 2008-09. All five
# designs are retained so the sensitivity to that choice stays visible.
CFG$crisis_from <- as.Date("2008-01-01")
CFG$crisis_to   <- as.Date("2009-12-01")
CFG$orig_pan_from <- as.Date("2020-01-01")   # the original, mistimed window
CFG$orig_pan_to   <- as.Date("2022-12-01")
CFG$za_pan_from   <- as.Date("2021-07-01")   # earliest of the 2021 break cluster
CFG$za_win_to     <- as.Date("2022-12-01")

CFG$dummy_sets <- c("none", "orig", "za_step", "za_window", "za_series")

# The design used for headline numbers. za_step is the faithful translation
# of the Zivot-Andrews alternative, which is a shift in intercept and trend.
CFG$headline_dummy <- "za_step"

# --- univariate tests ------------------------------------------------------
CFG$common_start <- as.Date("2003-12-01")  # window the primary VECMs use
CFG$kpss_lags    <- "long"
CFG$za_model     <- "both"
CFG$za_lag       <- 4L
CFG$min_obs      <- 40L

# --- reporting -------------------------------------------------------------
CFG$alpha_sig  <- 0.10   # threshold for "significantly error-correcting"
