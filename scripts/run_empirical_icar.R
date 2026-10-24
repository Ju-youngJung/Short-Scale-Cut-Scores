## =============================================================================
## scripts/run_empirical_icar.R -- empirical illustration: ICAR items (psych::ability)
##   1. calibrate the 2PL, report fit;  2. run the exact pipeline on the estimated bank;
##   3. score actual respondents with every rule;  4. bootstrap the calibration.
## Usage: restart R, setwd(<project>/csem), source("scripts/run_empirical_icar.R")
## =============================================================================
source("00_setup.R")
for (f in c("R/08_decisions.R", "R/09_screening.R", "R/10_final_summary.R", "R/11_final_report.R", "R/12_empirical.R")) source(f)

cfg <- list(version = "empirical_icar_v1", estimators = c("WLE", "EAP", "HU"),
            se_types = list(WLE = "classical", EAP = "psd", HU = "sandwich"),
            theta_grid = seq(-3, 3, by = 0.5), gh_nodes = 15L, eps_levels = c(0.10, 0.20), conds = c("up", "down", "guess"),
            alpha = 0.05, screen_alphas = c(0.05, 0.20),
            cuts = c(0.67, 1.28),                 # ~75th and 90th percentiles of N(0,1): illustrative selection cuts
            L = 4, grid_step = 0.02, prior_mean = 0, prior_sd = 1, eap_nodes = 61L, eap_lim = 6, hu_k = 1, bw_c = 4, width_lim = 6,
            B_boot = if (exists("icar_B")) icar_B else 50, cores = max(1L, parallel::detectCores() - 1L))
out_root <- "results/empirical_icar"; dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_root, "run_log.txt"), open = "wt"); sink(log_con, split = TRUE)
writeLines(c(paste("date:", format(Sys.time())), paste("R:", R.version.string), "config:", capture.output(dput(cfg))), file.path(out_root, "protocol_snapshot.txt"))

## ---- 1. data and calibration --------------------------------------------------
dat <- icar_load(); cat("ICAR items:", ncol(dat), " respondents:", nrow(dat), " complete cases:", sum(complete.cases(dat)), "\n")
cal <- calibrate_2pl(dat); a <- cal$a; b <- cal$b
bank <- data.frame(item = cal$items, a = round(a, 3), b = round(b, 3)); print(bank, row.names = FALSE)
write.csv(bank, file.path(out_root, "icar_bank_2pl.csv"), row.names = FALSE)
fs <- fit_summary(cal$mod)
if (!is.null(fs$M2)) { cat("\nM2 fit:\n"); print(fs$M2) }
if (!is.null(fs$itemfit)) { cat("\nItem fit (S-X2):\n"); print(fs$itemfit) }
if (!is.null(fs$Q3)) cat(sprintf("\nQ3: max |Q3| = %.3f; pairs with |Q3| > 0.2: %d of %d\n", max(abs(fs$Q3), na.rm = TRUE), sum(abs(fs$Q3) > 0.2, na.rm = TRUE), sum(is.finite(fs$Q3))))

## ---- 2. exact pipeline on the estimated bank (treated as known) ---------------
spec <- list(layout = "ICAR", J = length(a), a = a, b = b)
res <- run_cell_final(spec, cfg, out_root)
agg <- list(summary = res$summary, info = res$info, screen_size = res$screen_size)
ms <- with_pass_ratios(marginal_final(res$summary, res$gh), cfg)
write.csv(res$summary, file.path(out_root, "ALL_summary.csv"), row.names = FALSE)
write.csv(ms, file.path(out_root, "ALL_marginal_summary.csv"), row.names = FALSE)
write_report_final(agg, ms, cfg, file.path(out_root, "REPORT_icar_exact.md"))
plot_final(res$summary, file.path(out_root, "figures"), cfg, eps_show = 0.20, cond_show = "up")

## ---- 3. actual respondents (complete cases) -------------------------------------
sc_enum <- readRDS(file.path(out_root, cell_key(spec), "scores.rds"))
Y <- all_patterns(length(a)); ci <- exact_ci_setup(Y, a, b, cfg$alpha)
R <- dat[complete.cases(dat), ]
rs <- score_respondents(R, a, b, cfg, Y, ci, sc_enum$lz)
rt <- respondent_table(rs, cfg)
write.csv(rs, file.path(out_root, "respondent_scores.csv"), row.names = FALSE)
cat("\n==== Respondents (complete cases, N =", nrow(R), "): decisions by rule ====\n"); print(rt$decisions, row.names = FALSE, digits = 3)
cat("\n==== Disagreement between rules (proportion of respondents) ====\n"); print(rt$disagreement, row.names = FALSE, digits = 3)
cat("\n==== Screen pass rates among respondents ====\n"); print(round(rt$pass_rates, 3))
cat(sprintf("\nWLE: mean %.2f, SD %.2f; exact CI unbounded on one side for %.1f%% of respondents\n",
            mean(rs$wle), sd(rs$wle), 100 * mean(!is.finite(rs$exact_L) | !is.finite(rs$exact_U))))

## ---- 4. calibration bootstrap -----------------------------------------------------
if (cfg$B_boot > 0) {
  cat("\n==== Calibration bootstrap (B =", cfg$B_boot, ") ====\n")
  bt <- bootstrap_calibration(dat, R, cfg, B = cfg$B_boot)
  write.csv(bt, file.path(out_root, "bootstrap_calibration.csv"), row.names = FALSE)
  q <- apply(bt[, -1], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  cat("\nBootstrap 2.5% / 50% / 97.5% of marginal rates and respondent proportions:\n"); print(round(t(q), 3))
}
cat("\nDone. Outputs in", out_root, "\n"); sink(); close(log_con)
