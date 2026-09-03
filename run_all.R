#!/usr/bin/env Rscript
# =============================================================================
# run_all.R -- single entry point for the whole analysis
# =============================================================================
#   Rscript run_all.R              run every stage
#   Rscript run_all.R 03 04        run selected stages only
#
# Stages are ordered and each writes to outputs/<stage>/. Stage 04 reads
# stage 03's output, so running 04 alone requires 03 to have run before.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(urca)
  library(vars)
})

ROOT <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(dirname(f)) else normalizePath(getwd())
}, error = function(e) normalizePath(getwd()))
setwd(ROOT)

source("config.R")
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

STAGES <- c("01_data_audit", "02_integration", "03_grid", "04_findings", "05_combined")

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  sel <- STAGES[vapply(STAGES, function(s) any(startsWith(s, args)), logical(1))]
  if (!length(sel)) stop("No stage matches: ", paste(args, collapse = ", "),
                         "\nAvailable: ", paste(STAGES, collapse = ", "))
} else {
  sel <- STAGES
}

t0 <- Sys.time()
for (s in sel) {
  source(file.path("analysis", paste0(s, ".R")), local = new.env())
}

section("DONE")
cat("Stages run: ", paste(sel, collapse = ", "), "\n", sep = "")
cat("Elapsed: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n", sep = "")
cat("Results: ", normalizePath(CFG$out_dir), "\n", sep = "")
