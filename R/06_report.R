## =============================================================================
## 06_report.R -- report tables (markdown) and figures
## =============================================================================

md_table <- function(df, digits = 3) {
  f <- function(x) if (is.numeric(x)) ifelse(is.na(x), "", formatC(x, digits = digits, format = "f")) else as.character(x)
  d <- as.data.frame(lapply(df, f), stringsAsFactors = FALSE)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

write_report <- function(agg, ms, rs, cfg, path) {
  s <- agg$summary; g <- s[s$on_grid, ]; m0 <- g[g$eps == 0, ]
  L <- c("# Exact evaluation of reported CSEMs and confidence intervals", "", paste("Generated:", format(Sys.time())), "")
  ## Table 1: SE ratio at selected theta
  t1 <- m0[m0$estimator != "EXACT_CI" & m0$theta %in% c(-3, -2, 0, 2, 3), c("layout", "J", "estimator", "se_type", "theta", "se_ratio")]
  t1w <- reshape(t1, idvar = c("layout", "J", "estimator", "se_type"), timevar = "theta", direction = "wide")
  names(t1w) <- sub("se_ratio.", "theta=", names(t1w), fixed = TRUE)
  L <- c(L, "## Table 1. Mean reported SE / exact SD under the model", "", md_table(t1w[order(t1w$layout, t1w$estimator, t1w$se_type, t1w$J), ]), "")
  ## Table 2: marginal coverage and width under the model
  t2 <- ms[ms$eps == 0, c("layout", "J", "estimator", "se_type", "coverage", "width", "rmse", "mass_bound", "mass_unbounded")]
  L <- c(L, "## Table 2. N(0,1)-marginal coverage, expected width, RMSE (model)", "", md_table(t2[order(t2$layout, t2$J, t2$estimator), ]), "")
  ## Table 3: minimum conditional coverage over the theta grid
  t3 <- aggregate(coverage ~ layout + J + estimator + se_type + eps + cond, data = g, FUN = min)
  L <- c(L, "## Table 3. Minimum conditional coverage over the theta grid", "", md_table(t3[order(t3$layout, t3$J, t3$eps, t3$cond, t3$estimator), ]), "")
  ## Table 4: decisions at cut scores (marginal, model): point rule vs interval-based three-way rules
  dc <- grep("^(wrong3_|undet_|wrong_point_)", names(ms), value = TRUE)
  t4 <- ms[ms$eps == 0, c("layout", "J", "estimator", "se_type", dc)]
  L <- c(L, "## Table 4. Marginal decision rates at cut scores (model)",
         "wrong_point: estimate on the wrong side of the cut; wrong3 / undet: wrong / undetermined under the three-way interval rule", "",
         md_table(t4[order(t4$layout, t4$J, t4$estimator), ]), "")
  ## Table 4b: conditional decision rates near the cuts for WLE-Wald vs exact CI
  g4 <- g[g$eps == 0 & ((g$estimator == "WLE" & g$se_type == "classical") | g$estimator == "EXACT_CI"), ]
  for (ct in cfg$cuts) {
    k <- g4[abs(g4$theta - ct) <= 1.3, c("layout", "J", "estimator", "theta", paste0("wrong3_", ct), paste0("undet_", ct), paste0("wrong_point_", ct))]
    L <- c(L, paste0("## Table 4b. Conditional decision rates near cut ", ct, " (WLE Wald vs exact CI)"), "", md_table(k[order(k$layout, k$J, k$theta, k$estimator), ]), "")
  }
  ## Table 5: contamination, marginal coverage
  t5 <- ms[ms$eps > 0, c("layout", "J", "estimator", "se_type", "eps", "cond", "coverage", "rmse")]
  L <- c(L, "## Table 5. Marginal coverage and RMSE under contamination", "", md_table(t5[order(t5$layout, t5$J, t5$eps, t5$cond, t5$estimator), ]), "")
  ## Table 6: robust estimators: multiple roots and truncation
  L <- c(L, "## Table 6. Patterns with multiple roots / truncated / missing estimates", "", md_table(agg$roots), "")
  ## Table 7: descriptive convergence slopes
  L <- c(L, "## Table 7. Slope of log|SE ratio - 1| on log J (descriptive; SE = lack of fit)", "",
         md_table(rs[!is.na(rs$slope), ]), "")
  writeLines(L, path); invisible(path)
}

plot_report <- function(summ, path, alpha = 0.05, eps_show = 0.10, cut_show = NULL) {
  summ <- summ[summ$on_grid, ]; pdf(path, width = 12, height = 8); on.exit(dev.off())
  pan <- function(d, y, main, ylab, Js, jc, hline = NULL, ylim = NULL, lty = 1, add = FALSE) {
    if (!add) { if (is.null(ylim)) ylim <- range(d[[y]], na.rm = TRUE); if (!all(is.finite(ylim))) ylim <- c(0, 1)
      plot(NA, xlim = range(d$theta), ylim = ylim, xlab = expression(theta), ylab = ylab, main = main, cex.main = 0.85)
      if (!is.null(hline)) abline(h = hline, lty = 3, col = "grey40") }
    for (J in Js) { k <- d[d$J == J, ]; k <- k[order(k$theta), ]; if (nrow(k)) lines(k$theta, k[[y]], col = jc[as.character(J)], lwd = 2, lty = lty) }
  }
  for (lay in unique(summ$layout)) {
    s <- summ[summ$layout == lay, ]; Js <- sort(unique(s$J)); jc <- setNames(hcl.colors(length(Js), "Zissou 1"), Js)
    m <- s[s$eps == 0, ]
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
    for (e in list(c("ML", "classical"), c("WLE", "classical"), c("WLE", "sandwich"), c("MAP", "classical"), c("EAP", "psd"), c("HU", "sandwich")))
      pan(m[m$estimator == e[1] & m$se_type == e[2], ], "se_ratio", paste0(lay, ": ", e[1], " ", e[2], " SE / exact SD"), "ratio", Js, jc, 1)
    legend("topright", legend = paste("J =", Js), col = jc, lwd = 2, bty = "n", cex = 0.8)
    par(mfrow = c(2, 3))
    for (e in list(c("WLE", "classical"), c("EAP", "psd"), c("EXACT_CI", "exact")))
      pan(m[m$estimator == e[1] & m$se_type == e[2], ], "coverage", paste0(lay, ": coverage, ", e[1], " ", e[2], " (model)"), "coverage", Js, jc, 1 - alpha, c(0.8, 1))
    u <- s[s$eps == eps_show & s$cond == "up", ]
    for (e in list(c("WLE", "classical"), c("HU", "sandwich"), c("EXACT_CI", "exact")))
      pan(u[u$estimator == e[1] & u$se_type == e[2], ], "coverage", paste0(lay, ": coverage, ", e[1], " ", e[2], " (up, eps=", eps_show, ")"), "coverage", Js, jc, 1 - alpha, c(0, 1))
    legend("bottomleft", legend = paste("J =", Js), col = jc, lwd = 2, bty = "n", cex = 0.8)
    par(mfrow = c(1, 2))
    pan(m[m$estimator == "WLE" & m$se_type == "classical", ], "width", paste0(lay, ": expected width, Wald (solid) vs exact (dashed)"), "width", Js, jc, NULL, c(0, 8))
    pan(m[m$estimator == "EXACT_CI", ], "width", "", "", Js, jc, lty = 2, add = TRUE)
    if (!is.null(cut_show)) {
      uc <- paste0("undet_", cut_show); wc <- paste0("wrong3_", cut_show); pc <- paste0("wrong_point_", cut_show)
      kw <- m[m$estimator == "WLE" & m$se_type == "classical", ]; ke <- m[m$estimator == "EXACT_CI", ]
      pan(ke, uc, paste0(lay, ": cut ", cut_show, ": undetermined, exact (solid) vs WLE Wald (dashed)"), "probability", Js, jc, NULL, c(0, 1))
      pan(kw, uc, "", "", Js, jc, lty = 2, add = TRUE)
    }
    par(mfrow = c(1, 2))
    if (!is.null(cut_show)) {
      pan(ke, wc, paste0(lay, ": cut ", cut_show, ": wrong decisions, exact (solid) vs WLE Wald (dashed)"), "probability", Js, jc, NULL, c(0, 0.3))
      pan(kw, wc, "", "", Js, jc, lty = 2, add = TRUE)
      pan(kw, pc, paste0(lay, ": cut ", cut_show, ": wrong point classification, WLE"), "probability", Js, jc, NULL, c(0, 1))
    }
  }
  invisible(path)
}
