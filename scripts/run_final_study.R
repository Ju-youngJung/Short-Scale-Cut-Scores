## =============================================================================
## scripts/run_final_study.R -- FINAL DESIGN: cut-score decisions and aberrant
## responding in short tests (exact; one-shot; resumable)
##   RQ1 decisions: raw-score rule (TCC-mapped), WLE point, WLE Wald, EAP, exact interval;
##       exact accuracy/consistency; Rudner's normal approximation for comparison
##   RQ2 contamination: coverage and decisions; Huber as the robust alternative
##   RQ3 treatment: l_z* screening calibrated by its exact null distribution
## Usage: restart R, setwd(<project root>/csem), source("scripts/run_final_study.R")
## =============================================================================

source("00_setup.R")
for (f in c("R/08_decisions.R", "R/09_screening.R", "R/10_final_summary.R", "R/11_final_report.R")) source(f)

cfg <- list(
  version      = "final_design_v2",
  Js           = c(8, 12, 16, 20),
  layouts      = c("wide", "peaked"),
  estimators   = c("WLE", "EAP", "HU"),
  se_types     = list(WLE = "classical", EAP = "psd", HU = "sandwich"),
  theta_grid   = seq(-3, 3, by = 0.5),
  gh_nodes     = 15L,
  eps_levels   = c(0.10, 0.20),
  conds        = c("up", "down", "guess"),
  alpha        = 0.05,
  screen_alphas = c(0.05, 0.20),        # exact-null sizes of the l_z* screen (detection-cost trade-off)
  cuts         = c(0.25, 1.25, 2.25),
  L            = 4, grid_step = 0.02,
  prior_mean   = 0, prior_sd = 1, eap_nodes = 61L, eap_lim = 6,
  hu_k         = 1,                     # VERIFY against Schuster & Yuan (2011)
  bw_c         = 4,                     # unused (biweight excluded)
  width_lim    = 6,
  cores        = max(1L, parallel::detectCores() - 1L)
)
cfg$cells <- unlist(lapply(cfg$layouts, function(l) lapply(cfg$Js, function(J) list(layout = l, J = J))), recursive = FALSE)
out_root <- "results/final_v2"; dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_root, "run_log.txt"), open = "wt"); sink(log_con, split = TRUE)
f <- list.files("R", full.names = TRUE)
writeLines(c(paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), paste("R:", R.version.string), "config:", capture.output(dput(cfg)),
             "code md5:", paste(f, tools::md5sum(f))), file.path(out_root, "protocol_snapshot.txt"))
cat("Cores:", cfg$cores, " Cells:", length(cfg$cells), "\n")

## ---- 0. checks on a small cell -----------------------------------------------
bk <- make_bank(8, "wide"); Y8 <- all_patterns(8)
sc8 <- score_all(Y8, bk$a, bk$b, cfg); ci8 <- exact_ci_setup(Y8, bk$a, bk$b, cfg$alpha)
L8 <- lz_star_wle(Y8, sc8$WLE$theta, bk$a, bk$b); scr8 <- screen_pass(L8, sc8$WLE$theta, Y8, bk$a, bk$b, cfg$screen_alphas[1])
## exact size of the screen at each calibration point t: P_t(L < q_t), q_t = cutoff used for patterns with theta_hat = t
size_at_cal <- max(vapply(unique(round(sc8$WLE$theta, 10)), function(t) {
  pr <- exp(pattern_logprob(Y8, p2pl(t, bk$a, bk$b))); pr <- pr / sum(pr)
  q <- scr8$cutoff[round(sc8$WLE$theta, 10) == t][1]; sum(pr[is.finite(L8) & L8 < q]) }, numeric(1)))
lw8 <- max(abs(as.numeric(tapply(exp(pattern_logprob(Y8, p2pl(0.3, bk$a, bk$b))), factor(rowSums(Y8), levels = 0:8), sum)) - lord_wingersky(0.3, bk$a, bk$b)))
chk <- list(expfam = check_expfam(Y8, bk$a, bk$b, ci8), mlr = check_mlr(ci8), ci_coverage = check_ci_coverage(Y8, bk$a, bk$b, ci8, cfg$alpha),
            robust_moments = check_robust_moments(Y8, bk$a, bk$b, cfg), roots = check_roots(Y8, bk$a, bk$b, cfg, sc8),
            lord_wingersky = list(pass = lw8 < 1e-12, value = lw8, note = "raw-score distribution: enumeration = Lord-Wingersky"),
            lz_null_finite = list(pass = all(is.finite(L8)), value = sum(!is.finite(L8)), note = "l_z* finite for all patterns (J = 8)"),
            mirt = check_mirt(Y8, bk$a, bk$b, cfg, sc8))
run_checks(chk)
cat(sprintf("[REPORT] exact screen size at calibration points (should be <= %.2f): max = %.4f\n", cfg$screen_alphas[1], size_at_cal))
fails <- names(chk)[vapply(chk, function(r) isFALSE(r$pass), logical(1))]
if (length(fails)) { sink(); close(log_con); stop("Checks failed: ", paste(fails, collapse = ", ")) }

## ---- 1. cells ---------------------------------------------------------------
T0 <- Sys.time(); cells <- list()
for (spec in cfg$cells) cells[[cell_key(spec)]] <- run_cell_final(spec, cfg, out_root)

## ---- 2. aggregate, post-checks, report, figures --------------------------------
agg <- list(summary = do.call(rbind, lapply(cells, `[[`, "summary")), info = do.call(rbind, lapply(cells, `[[`, "info")),
            screen_size = do.call(rbind, lapply(cells, `[[`, "screen_size")))
write.csv(agg$summary, file.path(out_root, "ALL_summary.csv"), row.names = FALSE)
write.csv(agg$info, file.path(out_root, "ALL_cell_info.csv"), row.names = FALSE)
write.csv(agg$screen_size, file.path(out_root, "ALL_screen_size.csv"), row.names = FALSE)
m0 <- agg$summary[agg$summary$eps == 0 & agg$summary$estimator == "EXACT_CI", ]
post <- list(exact_ci_all_cells = list(pass = min(m0$coverage) >= 1 - cfg$alpha - 1e-10, value = min(m0$coverage), note = "min exact-CI coverage under the model, all cells"),
             lord_wingersky_all = list(pass = max(agg$info$lw_max_diff) < 1e-12, value = max(agg$info$lw_max_diff), note = "Lord-Wingersky check, all cells"),
             integrity = { ip <- agg$summary[agg$summary$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI"), c("coverage", grep("^wrong3_|^undet_", names(agg$summary), value = TRUE))]
                           list(pass = !anyNA(ip), value = sum(is.na(ip)), note = "no NA in coverage / three-way columns for interval procedures") })
run_checks(post)
ms <- with_pass_ratios(marginal_final(agg$summary, cells[[1]]$gh), cfg)
write.csv(ms, file.path(out_root, "ALL_marginal_summary.csv"), row.names = FALSE)
write_report_final(agg, ms, cfg, file.path(out_root, "REPORT_final.md"))
plot_final(agg$summary, file.path(out_root, "figures"), cfg, eps_show = 0.20, cond_show = "up")

cat("\n==== Cell information ====\n"); print(agg$info, row.names = FALSE)
cat("\n==== Marginal, model, cut 1.25: wrong (2-outcome) / consistency / wrong3 / undetermined ====\n")
k <- ms[ms$eps == 0 & ms$J %in% c(8, 20), c("layout", "J", "estimator", "wrong_point_1.25", "cons_point_1.25", "wrong3_1.25", "undet_1.25")]
print(k[order(k$layout, k$J, k$estimator), ], row.names = FALSE, digits = 3)
cat("\n==== Screening, eps = 0.2 up, J = 20: pass and exact-interval coverage all vs given pass, by screen size ====\n")
for (al in cfg$screen_alphas) { sx <- paste0("_s", al)
  k <- ms[ms$eps == 0.2 & ms$cond == "up" & ms$J == 20 & ms$estimator == "EXACT_CI", c("layout", paste0("pass", sx), "coverage", paste0("cov_given_pass", sx))]
  cat("screen size", al, "\n"); print(k, row.names = FALSE, digits = 3) }
cat(sprintf("\nTotal %.1f min. Outputs in %s: REPORT_final.md, ALL_*.csv, figures/\n", as.numeric(difftime(Sys.time(), T0, units = "mins")), out_root))
sink(); close(log_con)
