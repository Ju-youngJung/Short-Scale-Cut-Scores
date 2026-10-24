## =============================================================================
## 00_setup.R -- exact CSEM study: setup
##   setwd("<project root>"); source("00_setup.R")
## =============================================================================
for (d in c("R", "scripts", "results", "docs")) dir.create(d, showWarnings = FALSE)
if (!requireNamespace("mirt", quietly = TRUE)) message("Optional: install.packages('mirt') for external checks of WLE/EAP.")
for (f in c("R/01_model.R", "R/02_estimators.R", "R/03_exact_ci.R", "R/04_summary.R", "R/05_checks.R", "R/06_report.R", "R/07_figures_png.R", "R/08_decisions.R", "R/09_screening.R", "R/10_final_summary.R", "R/11_final_report.R")) {
  if (!file.exists(f)) stop("Missing file: ", f); source(f)
}
CSEM_STUDY_VERSION <- 1L
cat("Setup complete. Run: source('scripts/run_everything.R')\n")
