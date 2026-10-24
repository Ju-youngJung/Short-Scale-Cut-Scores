## =============================================================================
## 11_final_report.R (v2) -- tables (REPORT_final.md) and PNG figures for the final design
## =============================================================================

with_pass_ratios <- function(ms, cfg) {
  for (al in cfg$screen_alphas) { sx <- paste0("_s", al)
    ms[[paste0("cov_given_pass", sx)]] <- ms[[paste0("cov_pass", sx)]] / ms[[paste0("pass", sx)]]
    for (ct in cfg$cuts) {
      ms[[paste0("wrong3_given_pass", sx, "_", ct)]] <- ms[[paste0("wrong3_pass", sx, "_", ct)]] / ms[[paste0("pass", sx)]]
      ms[[paste0("undet_given_pass", sx, "_", ct)]]  <- ms[[paste0("undet_pass", sx, "_", ct)]]  / ms[[paste0("pass", sx)]]
      ms[[paste0("wrong_point_given_pass", sx, "_", ct)]] <- ms[[paste0("wrong_point_pass", sx, "_", ct)]] / ms[[paste0("pass", sx)]]
    }
  }
  ms
}

write_report_final <- function(agg, ms, cfg, path) {
  g <- agg$summary[agg$summary$on_grid, ]; cuts <- cfg$cuts
  L <- c("# Cut-score decisions and aberrant responding in short tests: exact results", "", paste("Generated:", format(Sys.time())), "",
         "## Table 0. Cell information (raw cuts mapped through the TCC; Lord-Wingersky check)", "", md_table(agg$info), "")
  m0 <- ms[ms$eps == 0, ]
  for (ct in cuts) {
    t1 <- m0[, c("layout", "J", "estimator", "se_type", paste0("wrong_point_", ct), paste0("cons_point_", ct), paste0("wrong3_", ct), paste0("undet_", ct))]
    names(t1) <- c("layout", "J", "estimator", "se_type", "wrong (2-outcome)", "consistency", "wrong (3-way)", "undetermined")
    t1$accuracy <- 1 - t1[["wrong (2-outcome)"]]
    L <- c(L, paste0("## Table 1 (cut ", ct, "). Marginal decision rates under the model (RAW_lo/RAW_hi: raw-score rule; RUDNER: normal-approximation index for the WLE point rule)"), "",
           md_table(t1[order(t1$layout, t1$J, t1$estimator), ]), "")
  }
  Js <- range(g$J)
  for (ct in cuts) {
    k <- g[g$eps == 0 & g$J %in% Js & abs(g$theta - ct) <= 1.3 & g$estimator %in% c("RAW_lo", "RAW_hi", "RUDNER", "WLE", "EAP", "EXACT_CI"),
           c("layout", "J", "theta", "estimator", "se_type", paste0("wrong_point_", ct), paste0("wrong3_", ct), paste0("undet_", ct))]
    L <- c(L, paste0("## Table 2 (cut ", ct, "). Conditional decision rates near the cut"), "", md_table(k[order(k$layout, k$J, k$theta, k$estimator), ]), "")
  }
  t3 <- ms[ms$eps > 0 & ms$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI", "RAW_hi"),
           c("layout", "J", "eps", "cond", "estimator", "coverage", "rmse", paste0("wrong_point_", cuts[1]), paste0("wrong3_", cuts[1]), paste0("undet_", cuts[1]))]
  L <- c(L, "## Table 3. Marginal coverage, RMSE and decision rates (first cut) under contamination", "", md_table(t3[order(t3$layout, t3$J, t3$eps, t3$cond, t3$estimator), ]), "")
  for (al in cfg$screen_alphas) { sx <- paste0("_s", al)
    t4 <- ms[ms$estimator %in% c("WLE", "EXACT_CI"), c("layout", "J", "eps", "cond", "estimator", paste0("pass", sx), "coverage", paste0("cov_given_pass", sx),
                                                        paste0("wrong3_", cuts[1]), paste0("wrong3_given_pass", sx, "_", cuts[1]), paste0("undet_", cuts[1]), paste0("undet_given_pass", sx, "_", cuts[1]))]
    L <- c(L, paste0("## Table 4 (screen size ", al, "). l_z* screening with exact-null plug-in cutoff: pass rate, coverage and first-cut decisions, all vs given pass"), "",
           md_table(t4[order(t4$layout, t4$J, t4$eps, t4$cond, t4$estimator), ]), "")
  }
  ssz <- agg$screen_size[, c("layout", "J", "screen_alpha", "theta", "size")]
  L <- c(L, "## Table 5. Exact conditional size of the screen under the model, P(fail | theta)", "",
         md_table(reshape(ssz, idvar = c("layout", "J", "screen_alpha"), timevar = "theta", direction = "wide")), "")
  writeLines(L, path); invisible(path)
}

plot_final <- function(summ, out_dir, cfg, eps_show = 0.20, cond_show = "up") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE); s <- summ[summ$on_grid, ]
  banks <- unique(s$layout); Js <- sort(unique(s$J)); jc <- setNames(hcl.colors(length(Js), "Zissou 1"), Js); th <- fig_theme
  curve <- function(d, y, lty = 1) for (J in Js) { k <- d[d$J == J, ]; k <- k[order(k$theta), ]; if (nrow(k)) lines(k$theta, k[[y]], col = jc[as.character(J)], lwd = th$lwd, lty = lty) }
  open_png <- function(f, ncol, nrow = 1) { png(f, width = th$width_in, height = nrow * th$height_per_row + 0.6, units = "in", res = th$res)
    layout(rbind(matrix(1:(nrow * ncol), nrow, ncol, byrow = TRUE), rep(nrow * ncol + 1, ncol)), heights = c(rep(1, nrow), 0.22)); par(mar = th$mar, mgp = th$mgp) }
  frame <- function(main, ylab, ylim, hline = NULL) { plot(NA, xlim = c(-3, 3), ylim = ylim, xlab = expression(theta), ylab = ylab, main = main, cex.main = th$cex_main, cex.axis = th$cex_axis, cex.lab = th$cex_lab, las = 1)
    if (!is.null(hline)) abline(h = hline, lty = 3, col = "grey45") }
  leg <- function(extra = NULL, lty_extra = NULL) { par(mar = c(0, 0, 0, 0)); plot.new()
    legend("center", legend = c(paste("J =", Js), extra), col = c(jc, rep("black", length(extra))), lty = c(rep(1, length(Js)), lty_extra), lwd = th$lwd, horiz = TRUE, bty = "n", cex = 0.85, seg.len = 2.5) }
  pick <- function(est, se = NULL, eps = 0, cond = "none", lay) { d <- s[s$estimator == est & s$eps == eps & s$cond == cond & s$layout == lay, ]; if (!is.null(se)) d <- d[d$se_type == se, ]; d }
  for (lay in banks) {
    tag <- gsub("[^A-Za-z0-9]+", "_", lay); m0 <- function(e, se = NULL) pick(e, se, lay = lay)
    f <- file.path(out_dir, paste0("G1_wrong_two_outcome_", tag, ".png")); open_png(f, length(cfg$cuts))
    for (ct in cfg$cuts) { v <- paste0("wrong_point_", ct); frame(paste0(lay, ": wrong decision, cut ", ct), "probability", c(0, 1))
      curve(m0("RAW_hi"), v, 1); curve(m0("WLE", "classical"), v, 2); curve(m0("RUDNER"), v, 3) }
    leg(c("raw score (k+)", "WLE point", "Rudner approx."), c(1, 2, 3)); dev.off()
    f <- file.path(out_dir, paste0("G2_three_way_", tag, ".png")); open_png(f, length(cfg$cuts), 2)
    for (ct in cfg$cuts) { v <- paste0("undet_", ct); frame(paste0(lay, ": undetermined, cut ", ct), "probability", c(0, 1))
      curve(m0("EXACT_CI"), v, 1); curve(m0("WLE", "classical"), v, 2); curve(m0("EAP"), v, 3) }
    for (ct in cfg$cuts) { v <- paste0("wrong3_", ct); frame(paste0(lay, ": wrong 3-way decision, cut ", ct), "probability", c(0, 0.3), cfg$alpha / 2)
      curve(m0("EXACT_CI"), v, 1); curve(m0("WLE", "classical"), v, 2); curve(m0("EAP"), v, 3) }
    leg(c("exact", "WLE Wald", "EAP"), c(1, 2, 3)); dev.off()
    f <- file.path(out_dir, paste0("G3_coverage_", cond_show, eps_show, "_", tag, ".png")); open_png(f, 4)
    for (e in list(c("WLE", "classical"), c("EAP", "psd"), c("HU", "sandwich"), c("EXACT_CI", "exact"))) {
      frame(paste0(lay, ": ", e[1], " (", cond_show, ", eps = ", eps_show, ")"), "coverage", c(0, 1), 1 - cfg$alpha); curve(pick(e[1], e[2], eps_show, cond_show, lay), "coverage") }
    leg(); dev.off()
    for (al in cfg$screen_alphas) { sx <- paste0("_s", al)
      f <- file.path(out_dir, paste0("G4_screening_s", al, "_", cond_show, eps_show, "_", tag, ".png")); open_png(f, 3)
      d <- pick("EXACT_CI", NULL, eps_show, cond_show, lay); d$cgp <- d[[paste0("cov_pass", sx)]] / d[[paste0("pass", sx)]]
      frame(paste0(lay, ": exact interval, ", cond_show, " eps = ", eps_show, " (screen ", al, "): coverage all (dashed) vs given pass (solid)"), "coverage", c(0, 1), 1 - cfg$alpha)
      curve(d, "cgp", 1); curve(d, "coverage", 2)
      frame(paste0(lay, ": pass rate under contamination"), "P(pass)", c(0, 1)); curve(d, paste0("pass", sx))
      d0 <- m0("EXACT_CI"); frame(paste0(lay, ": pass rate under the model"), "P(pass)", c(0.5, 1), 1 - al); curve(d0, paste0("pass", sx))
      leg(); dev.off()
    }
  }
  invisible(out_dir)
}
