## =============================================================================
## scripts/run_everything.R -- ONE-SHOT: exact study, empirical ICAR study, figures,
## and the supplementary workbook. Restart R, then:
##   setwd("<project folder>");  source("scripts/run_everything.R")
## Outputs:
##   results/final_v2/        ALL_summary.csv, ALL_marginal_summary.csv, ALL_screen_size.csv,
##                            ALL_cell_info.csv, REPORT_final.md, run_log.txt, figures/*.png
##   results/empirical_icar/  icar_bank_2pl.csv, ALL_summary.csv, ALL_marginal_summary.csv,
##                            respondent_scores.csv, bootstrap_calibration.csv, run_log.txt, figures/*.png
##   results/Supplementary_Materials.xlsx
## Runtime: about 15 min (exact study) + 2 min (ICAR) + bootstrap (B = 50, ~30-60 min; set B_BOOT below)
## Requirements: packages mirt, psychTools, openxlsx (installed automatically if missing)
## =============================================================================
t_all <- Sys.time()
for (p in c("mirt", "psychTools", "openxlsx")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
B_BOOT <- 50L                        # calibration bootstrap draws for ICAR (0 to skip)

## 1. exact (enumeration) study -------------------------------------------------
cat("\n==== 1/4 Exact study ====\n")
source("scripts/run_final_study.R")              # resumable; skips completed cells

## 2. empirical ICAR study ------------------------------------------------------
cat("\n==== 2/4 Empirical ICAR study ====\n")
rm(list = setdiff(ls(), c("t_all", "B_BOOT")))   # clean workspace between studies
icar_B <- B_BOOT
source("scripts/run_empirical_icar.R")

## 3. publication figures from the final summary --------------------------------
cat("\n==== 3/4 Figures ====\n")
rm(list = setdiff(ls(), c("t_all", "B_BOOT")))
in_file <- "results/final_v2/ALL_summary.csv"; out_dir <- "results/final_v2/figures"
source("scripts/make_final_figures.R")

## 4. supplementary workbook ------------------------------------------------------
cat("\n==== 4/4 Supplementary workbook ====\n")
rm(list = setdiff(ls(), c("t_all", "B_BOOT")))
root <- getwd(); out_file <- file.path(root, "results", "Supplementary_Materials.xlsx"); refit_icar <- TRUE
source("scripts/make_supplement_xlsx.R")

cat(sprintf("\nALL DONE in %.1f min. Figures for the manuscript:\n", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
print(list.files("results/final_v2/figures", pattern = "^G[123]_", full.names = TRUE))
print(list.files("results/empirical_icar/figures", pattern = "^G[123]_", full.names = TRUE))
