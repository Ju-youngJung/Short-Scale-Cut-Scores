## =============================================================================
## 07_figures_png.R -- publication figures as separate PNG files
##
## Design rules
##  * one PNG per figure; every panel in a figure shares the same y-axis limits,
##    computed from the data shown in that figure (optionally across banks too);
##  * colours encode J (one palette for the whole paper), line type encodes method;
##  * a single legend row below the panels;
##  * data are reshaped once; each figure is a small spec evaluated by save_figure().
## Usage (after source("00_setup.R")):
##   summ <- read.csv("results/csem/ALL_summary.csv")
##   make_all_pngs(summ, "results/csem/figures", common_across_banks = TRUE)
## =============================================================================

fig_theme <- list(width_in = 10, height_per_row = 3.1, res = 300, lwd = 1.8, cex_main = 0.9,
                  cex_axis = 0.85, cex_lab = 0.95, mar = c(3.6, 3.8, 2.4, 0.8), mgp = c(2.2, 0.7, 0))

j_palette <- function(Js) setNames(hcl.colors(length(Js), "Zissou 1"), Js)

## one series = subset of rows (estimator, se_type, eps, cond) -> curve per J
series_df <- function(summ, estimator, se_type, eps = 0, cond = "none", layout = NULL) {
  d <- summ[summ$on_grid & summ$estimator == estimator & summ$se_type == se_type & summ$eps == eps & summ$cond == cond, ]
  if (!is.null(layout)) d <- d[d$layout == layout, ]
  d
}

draw_panel <- function(d_list, yvar, main, ylab, ylim, Js, jc, ltys, hline = NULL, th = fig_theme) {
  plot(NA, xlim = c(-3, 3), ylim = ylim, xlab = expression(theta), ylab = ylab, main = main,
       cex.main = th$cex_main, cex.axis = th$cex_axis, cex.lab = th$cex_lab, las = 1)
  if (!is.null(hline)) abline(h = hline, lty = 3, col = "grey45")
  for (i in seq_along(d_list)) for (J in Js) {
    k <- d_list[[i]][d_list[[i]]$J == J, ]; k <- k[order(k$theta), ]
    if (nrow(k)) lines(k$theta, k[[yvar]], col = jc[as.character(J)], lwd = th$lwd, lty = ltys[i])
  }
}

## panels: list of list(series = list(df, df, ...), main, hline); all share yvar and ylim
save_figure <- function(file, panels, yvar, ylab, Js, ltys = 1, lty_labels = NULL, ylim = NULL,
                        ncol = 3, pad = 0.02, th = fig_theme) {
  nrow <- ceiling(length(panels) / ncol); jc <- j_palette(Js)
  if (is.null(ylim)) {
    v <- unlist(lapply(panels, function(p) lapply(p$series, function(d) d[[yvar]])))
    v <- v[is.finite(v)]; r <- range(v); ylim <- r + c(-1, 1) * pad * diff(r)
  }
  png(file, width = th$width_in, height = nrow * th$height_per_row + 0.6, units = "in", res = th$res)
  on.exit(dev.off())
  layout(rbind(matrix(seq_len(nrow * ncol), nrow, ncol, byrow = TRUE), rep(nrow * ncol + 1, ncol)),
         heights = c(rep(1, nrow), 0.22))
  par(mar = th$mar, mgp = th$mgp)
  for (p in panels) draw_panel(p$series, yvar, p$main, ylab, ylim, Js, jc, ltys, p$hline, th)
  if (length(panels) < nrow * ncol) for (i in seq_len(nrow * ncol - length(panels))) plot.new()
  par(mar = c(0, 0, 0, 0)); plot.new()
  leg <- paste("J =", Js); lcol <- jc; llty <- rep(1, length(Js))
  if (!is.null(lty_labels)) { leg <- c(leg, lty_labels); lcol <- c(lcol, rep("black", length(lty_labels))); llty <- c(llty, ltys) }
  legend("center", legend = leg, col = lcol, lty = llty, lwd = th$lwd, horiz = TRUE, bty = "n", cex = 0.85, seg.len = 2.5)
  invisible(file)
}

## ---------------------------------------------------------------------------
## Figure specifications
## ---------------------------------------------------------------------------
se_panels <- list(c("ML", "classical", "ML, classical SE"), c("WLE", "classical", "WLE, classical SE"),
                  c("WLE", "sandwich", "WLE, sandwich SE"), c("MAP", "classical", "MAP, classical SE"),
                  c("EAP", "psd", "EAP, posterior SD"), c("HU", "sandwich", "Huber, sandwich SE"))

make_all_pngs <- function(summ, out_dir, common_across_banks = TRUE, eps_show = 0.10, cond_show = "up",
                          cuts = NULL, alpha = 0.05) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  banks <- setdiff(unique(summ$layout), grep("a = ", unique(summ$layout), value = TRUE))
  Js <- sort(unique(summ$J[summ$layout %in% banks]))
  if (is.null(cuts)) cuts <- as.numeric(sub("undet_", "", grep("^undet_", names(summ), value = TRUE)))
  index <- list()
  ## shared y-limits across banks, if requested, computed per figure type
  lim <- function(fn) if (common_across_banks) { v <- unlist(lapply(banks, fn)); v <- v[is.finite(v)]; r <- range(v); r + c(-1, 1) * 0.02 * diff(r) } else NULL
  ylim_se  <- lim(function(b) unlist(lapply(se_panels, function(p) series_df(summ, p[1], p[2], layout = b)$se_ratio)))
  ylim_cov <- c(0.80, 1.00); ylim_cov_c <- c(0, 1)
  ylim_w   <- lim(function(b) c(series_df(summ, "WLE", "classical", layout = b)$width, series_df(summ, "EXACT_CI", "exact", layout = b)$width))
  for (b in banks) {
    tag <- gsub("[^A-Za-z0-9]+", "_", b)
    ## F1: SE accuracy
    f <- file.path(out_dir, paste0("F1_se_ratio_", tag, ".png"))
    save_figure(f, lapply(se_panels, function(p) list(series = list(series_df(summ, p[1], p[2], layout = b)), main = paste0(b, ": ", p[3]), hline = 1)),
                "se_ratio", "reported SE / exact SD", Js, ylim = ylim_se); index[[f]] <- "SE accuracy under the model"
    ## F2: coverage under the model
    f <- file.path(out_dir, paste0("F2_coverage_model_", tag, ".png"))
    save_figure(f, list(list(series = list(series_df(summ, "WLE", "classical", layout = b)), main = paste0(b, ": WLE Wald"), hline = 1 - alpha),
                        list(series = list(series_df(summ, "EAP", "psd", layout = b)), main = paste0(b, ": EAP posterior-SD interval"), hline = 1 - alpha),
                        list(series = list(series_df(summ, "EXACT_CI", "exact", layout = b)), main = paste0(b, ": exact interval"), hline = 1 - alpha)),
                "coverage", "coverage probability", Js, ylim = ylim_cov); index[[f]] <- "Coverage under the model"
    ## F3: coverage under contamination
    f <- file.path(out_dir, paste0("F3_coverage_", cond_show, eps_show, "_", tag, ".png"))
    save_figure(f, list(list(series = list(series_df(summ, "WLE", "classical", eps_show, cond_show, b)), main = paste0(b, ": WLE Wald (", cond_show, ", eps = ", eps_show, ")"), hline = 1 - alpha),
                        list(series = list(series_df(summ, "HU", "sandwich", eps_show, cond_show, b)), main = paste0(b, ": Huber sandwich"), hline = 1 - alpha),
                        list(series = list(series_df(summ, "EXACT_CI", "exact", eps_show, cond_show, b)), main = paste0(b, ": exact interval"), hline = 1 - alpha)),
                "coverage", "coverage probability", Js, ylim = ylim_cov_c); index[[f]] <- "Coverage under contamination"
    ## F4: width, Wald vs exact (same panel, two line types)
    f <- file.path(out_dir, paste0("F4_width_", tag, ".png"))
    save_figure(f, list(list(series = list(series_df(summ, "WLE", "classical", layout = b), series_df(summ, "EXACT_CI", "exact", layout = b)),
                             main = paste0(b, ": expected width"), hline = NULL)),
                "width", "expected width", Js, ltys = c(1, 2), lty_labels = c("WLE Wald", "exact"), ylim = ylim_w, ncol = 1); index[[f]] <- "Expected width"
    ## F5: decisions at each cut: undetermined and wrong, exact vs WLE Wald, plus point rule
    for (ct in cuts) {
      f <- file.path(out_dir, paste0("F5_decisions_cut", ct, "_", tag, ".png"))
      w <- series_df(summ, "WLE", "classical", layout = b); e <- series_df(summ, "EXACT_CI", "exact", layout = b)
      ## the three panels use different variables: drawn explicitly; probability scale [0,1] shared for undetermined/point, [0,0.3] for wrong decisions
      png(f, width = fig_theme$width_in, height = fig_theme$height_per_row + 0.6, units = "in", res = fig_theme$res)
      layout(rbind(1:3, c(4, 4, 4)), heights = c(1, 0.22)); par(mar = fig_theme$mar, mgp = fig_theme$mgp); jc <- j_palette(Js)
      draw_panel(list(e, w), paste0("undet_", ct), paste0(b, ": undetermined, cut ", ct), "probability", c(0, 1), Js, jc, c(1, 2))
      draw_panel(list(e, w), paste0("wrong3_", ct), paste0(b, ": wrong 3-way decision, cut ", ct), "probability", c(0, 0.3), Js, jc, c(1, 2), hline = alpha / 2)
      draw_panel(list(w), paste0("wrong_point_", ct), paste0(b, ": wrong point classification (WLE), cut ", ct), "probability", c(0, 1), Js, jc, 1)
      par(mar = c(0, 0, 0, 0)); plot.new()
      legend("center", legend = c(paste("J =", Js), "exact", "WLE Wald"), col = c(jc, "black", "black"), lty = c(rep(1, length(Js)), 1, 2),
             lwd = fig_theme$lwd, horiz = TRUE, bty = "n", cex = 0.85, seg.len = 2.5)
      dev.off(); index[[f]] <- paste("Decisions at cut", ct)
    }
  }
  ## F6: cross-bank comparison of the WLE SE ratio and exact/Wald coverage (shared scale)
  f <- file.path(out_dir, "F6_banks_WLE_se_ratio.png")
  save_figure(f, lapply(banks, function(b) list(series = list(series_df(summ, "WLE", "classical", layout = b)), main = paste0(b, ": WLE classical SE / exact SD"), hline = 1)),
              "se_ratio", "reported SE / exact SD", Js, ylim = ylim_se, ncol = length(banks)); index[[f]] <- "WLE SE accuracy by bank"
  idx <- data.frame(file = basename(names(index)), content = unlist(index), row.names = NULL)
  write.csv(idx, file.path(out_dir, "figure_index.csv"), row.names = FALSE); idx
}
