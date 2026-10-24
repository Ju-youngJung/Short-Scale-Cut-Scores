## =============================================================================
## make_final_figures.R -- exact-study figures from ALL_summary.csv (standalone)
## Everything (banks, J, cuts, screen sizes) is read from the CSV, so the figures
## always match the file. Panels within a figure share their y-axis range.
## Usage: setwd(<project>); source("scripts/make_final_figures.R")
## =============================================================================

if (!exists("in_file")) in_file <- "results/final_v2/ALL_summary.csv"   # final-design summary
if (!exists("out_dir")) out_dir <- "results/final_v2/figures"
eps_show <- 0.20; cond_show <- "up"; alpha <- 0.05

summ <- read.csv(in_file, check.names = FALSE, stringsAsFactors = FALSE)
summ$on_grid <- as.logical(summ$on_grid)
s <- summ[summ$on_grid, ]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- what the file contains --------------------------------------------------
banks  <- unique(s$layout)
Js     <- sort(unique(s$J))
cuts   <- sort(as.numeric(sub("^undet_", "", grep("^undet_[0-9.]+$", names(s), value = TRUE))))
sizes  <- sub("^pass_s", "", grep("^pass_s[0-9.]+$", names(s), value = TRUE))
cat("Banks:", banks, "| J:", Js, "| cuts:", cuts, "| screen sizes:", if (length(sizes)) sizes else "none", "\n")
if (!length(cuts)) stop("No decision columns (undet_*) found: is this the final-design summary?")

## ---- plotting helpers ---------------------------------------------------------
th  <- list(w = 10, h = 3.1, res = 300, lwd = 1.8, mar = c(3.6, 3.8, 2.4, 0.8), mgp = c(2.2, 0.7, 0))
jc  <- setNames(hcl.colors(length(Js), "Zissou 1"), Js)
pick <- function(est, se = NULL, eps = 0, cond = "none", lay) {
  d <- s[s$estimator == est & s$eps == eps & s$cond == cond & s$layout == lay, ]
  if (!is.null(se)) d <- d[d$se_type == se, ]
  d
}
curves <- function(d, y, lty = 1, label = "") {
  if (!nrow(d) || !(y %in% names(d)) || all(is.na(d[[y]]))) {
    cat("  WARNING: no data for", y, if (nzchar(label)) paste0("(", label, ")"), "\n")
    text(0, mean(par("usr")[3:4]), "no data in file", col = "grey40"); return(invisible(FALSE)) }
  for (J in Js) { k <- d[d$J == J, ]; k <- k[order(k$theta), ]
    if (nrow(k)) lines(k$theta, k[[y]], col = jc[as.character(J)], lwd = th$lwd, lty = lty) }
  invisible(TRUE)
}
frame <- function(main, ylab, ylim, hline = NULL) {
  plot(NA, xlim = c(-3, 3), ylim = ylim, xlab = expression(theta), ylab = ylab, main = main,
       cex.main = 0.9, cex.axis = 0.85, cex.lab = 0.95, las = 1)
  if (!is.null(hline)) abline(h = hline, lty = 3, col = "grey45")
}
open_png <- function(f, ncol, nrow = 1) {
  png(f, width = th$w, height = nrow * th$h + 0.6, units = "in", res = th$res)
  layout(rbind(matrix(seq_len(nrow * ncol), nrow, ncol, byrow = TRUE), rep(nrow * ncol + 1, ncol)), heights = c(rep(1, nrow), 0.22))
  par(mar = th$mar, mgp = th$mgp)
}
legend_row <- function(extra = NULL, lty_extra = NULL) {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = c(paste("J =", Js), extra), col = c(jc, rep("black", length(extra))),
         lty = c(rep(1, length(Js)), lty_extra), lwd = th$lwd, horiz = TRUE, bty = "n", cex = 0.85, seg.len = 2.5)
}
written <- character(0)
save_as <- function(f) { dev.off(); written <<- c(written, basename(f)) }

## ---- figures per bank ---------------------------------------------------------
for (lay in banks) {
  tag <- gsub("[^A-Za-z0-9]+", "_", lay)
  ## G1: wrong decisions of two-outcome rules
  f <- file.path(out_dir, paste0("G1_wrong_two_outcome_", tag, ".png")); open_png(f, length(cuts))
  for (ct in cuts) { v <- paste0("wrong_point_", ct)
    frame(paste0(lay, ": wrong decision, cut ", ct), "probability", c(0, 1))
    curves(pick("RAW_hi", lay = lay), v, 1); curves(pick("WLE", "classical", lay = lay), v, 2); curves(pick("RUDNER", lay = lay), v, 3) }
  legend_row(c("raw score (k+)", "WLE point", "Rudner approx."), c(1, 2, 3)); save_as(f)
  ## G2: three-way rules
  f <- file.path(out_dir, paste0("G2_three_way_", tag, ".png")); open_png(f, length(cuts), 2)
  for (ct in cuts) { v <- paste0("undet_", ct); frame(paste0(lay, ": undetermined, cut ", ct), "probability", c(0, 1))
    curves(pick("EXACT_CI", lay = lay), v, 1); curves(pick("WLE", "classical", lay = lay), v, 2); curves(pick("EAP", lay = lay), v, 3) }
  for (ct in cuts) { v <- paste0("wrong3_", ct); frame(paste0(lay, ": wrong 3-way decision, cut ", ct), "probability", c(0, 0.3), alpha / 2)
    curves(pick("EXACT_CI", lay = lay), v, 1); curves(pick("WLE", "classical", lay = lay), v, 2); curves(pick("EAP", lay = lay), v, 3) }
  legend_row(c("exact", "WLE Wald", "EAP"), c(1, 2, 3)); save_as(f)
  ## G3: coverage under contamination
  f <- file.path(out_dir, paste0("G3_coverage_", cond_show, eps_show, "_", tag, ".png")); open_png(f, 4)
  for (e in list(c("WLE", "classical", "WLE Wald"), c("EAP", "psd", "EAP"), c("HU", "sandwich", "Huber"), c("EXACT_CI", "exact", "exact"))) {
    frame(paste0(lay, ": ", e[3], " (", cond_show, ", eps = ", eps_show, ")"), "coverage", c(0, 1), 1 - alpha)
    curves(pick(e[1], e[2], eps_show, cond_show, lay), "coverage") }
  legend_row(); save_as(f)
  ## G4: screening, one figure per screen size present in the file
  for (sz in sizes) {
    pc <- paste0("pass_s", sz); cc <- paste0("cov_pass_s", sz)
    d  <- pick("EXACT_CI", NULL, eps_show, cond_show, lay); d0 <- pick("EXACT_CI", lay = lay)
    if (!nrow(d) || !all(c(pc, cc) %in% names(d))) { cat("  skip G4 for size", sz, "in", lay, "(columns or rows missing)\n"); next }
    d$cgp <- d[[cc]] / d[[pc]]
    cat(sprintf("  G4 %s size %s: %d rows, %d non-missing pass values\n", lay, sz, nrow(d), sum(!is.na(d[[pc]]))))
    f <- file.path(out_dir, paste0("G4_screening_s", sz, "_", cond_show, eps_show, "_", tag, ".png")); open_png(f, 3)
    frame(paste0(lay, ": exact interval, screen ", sz, ": all (dashed) vs pass (solid)"), "coverage", c(0, 1), 1 - alpha)
    curves(d, "cgp", 1, paste(lay, "G4 cov|pass")); curves(d, "coverage", 2, paste(lay, "G4 coverage"))
    frame(paste0(lay, ": pass rate under contamination"), "P(pass)", c(0, 1)); curves(d, pc, 1, paste(lay, "G4 pass contam"))
    frame(paste0(lay, ": pass rate under the model"), "P(pass)", c(0.5, 1), 1 - as.numeric(sz)); curves(d0, pc, 1, paste(lay, "G4 pass model"))
    legend_row(); save_as(f)
  }
}
cat("Written to", out_dir, ":\n"); print(written)
