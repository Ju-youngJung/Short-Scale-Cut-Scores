## =============================================================================
## 04_summary.R -- exact performance measures for one cell
## =============================================================================

## three-way decisions from an interval [lower, upper] at cut ct
dec3 <- function(th, pr, lower, upper, cuts) {
  out <- c()
  for (ct in cuts) {
    pa <- sum(pr[which(lower > ct)]); pb <- sum(pr[which(upper < ct)])
    out[paste0("wrong3_", ct)] <- if (th >= ct) pb else pa
    out[paste0("undet_", ct)]  <- 1 - pa - pb
  }
  out
}
## point classification (estimate >= cut); missing estimates count as misclassified
dec_point <- function(th, pr, est, cuts) {
  out <- c()
  pna <- sum(pr[!is.finite(est)])
  for (ct in cuts) {
    pa <- sum(pr[which(est >= ct)]); pb <- sum(pr[which(est < ct)])
    out[paste0("wrong_point_", ct)] <- if (th >= ct) pb + pna else pa + pna
  }
  out
}

summarize_theta <- function(th, Y, a, b, sc, ci, cfg) {
  q <- item_q(th, a, b); z <- qnorm(1 - cfg$alpha / 2); wl <- cfg$width_lim; rows <- list()
  cut_cols <- c(paste0("wrong3_", cfg$cuts), paste0("undet_", cfg$cuts), paste0("wrong_point_", cfg$cuts))
  blank <- setNames(rep(NA_real_, length(cut_cols)), cut_cols)
  for (eps in c(0, cfg$eps_levels)) for (cond in if (eps == 0) "none" else cfg$conds) {
    cs <- contam_spec(cond, eps, q$P, q$Q)
    pr <- exp(pattern_logprob(Y, (1 - cs$e) * q$P + cs$e * cs$g)); pr <- pr / sum(pr)
    ## exact CI (does not depend on the estimator)
    inci <- ci$L <= th & th <= ci$U
    wid  <- pmax(0, pmin(ci$U, wl) - pmax(ci$L, -wl))
    unb  <- !is.finite(ci$L) | !is.finite(ci$U) | ci$L < -wl | ci$U > wl
    d <- blank; d[names(dec3(th, pr, ci$L, ci$U, cfg$cuts))] <- dec3(th, pr, ci$L, ci$U, cfg$cuts)
    rows[[length(rows) + 1]] <- data.frame(theta = th, eps = eps, cond = cond, estimator = "EXACT_CI", se_type = "exact",
      mass_bound = NA, mass_multi = NA, mass_na = 0, bias = NA, sd = NA, rmse = NA, se_ratio = NA,
      coverage = sum(pr * inci), width = sum(pr * wid), mass_unbounded = sum(pr * unb), t(d), check.names = FALSE)
    for (nm in names(sc)) {
      est <- sc[[nm]]$theta; fin <- is.finite(est); w <- pr[fin] / sum(pr[fin])
      mu <- sum(w * est[fin]); sdv <- sqrt(sum(w * (est[fin] - mu)^2))
      mb <- sum(pr[fin & abs(est) >= cfg$L - 1e-8])
      mm <- if (all(is.na(sc[[nm]]$nroots))) NA_real_ else sum(pr[which(sc[[nm]]$nroots > 1)])
      dp <- dec_point(th, pr, est, cfg$cuts)
      for (ty in names(sc[[nm]]$se)) {
        se <- sc[[nm]]$se[[ty]]; ok <- fin & is.finite(se)
        lo <- est - z * se; up <- est + z * se
        cov <- sum(pr[ok & lo <= th & th <= up])
        d <- blank
        d3 <- dec3(th, pr, ifelse(ok, lo, -Inf), ifelse(ok, up, Inf), cfg$cuts)   # missing SE: no decision
        d[names(d3)] <- d3; d[names(dp)] <- dp
        rows[[length(rows) + 1]] <- data.frame(theta = th, eps = eps, cond = cond, estimator = nm, se_type = ty,
          mass_bound = mb, mass_multi = mm, mass_na = sum(pr[!fin]), bias = mu - th, sd = sdv,
          rmse = sqrt((mu - th)^2 + sdv^2), se_ratio = sum(pr[ok] * se[ok]) / sum(pr[ok]) / sdv,
          coverage = cov, width = sum(pr[ok] * 2 * z * se[ok]) / sum(pr[ok]), mass_unbounded = NA, t(d), check.names = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

cell_key <- function(spec) sprintf("%s%s_J%02d", spec$layout, if (is.null(spec$a_const)) "" else paste0("_a", spec$a_const), spec$J)
cell_label <- function(spec) if (is.null(spec$a_const)) spec$layout else paste0(spec$layout, " (a = ", spec$a_const, ")")

run_cell <- function(spec, cfg, out_root) {
  key <- cell_key(spec); dir <- file.path(out_root, key); done <- file.path(dir, "cell.rds")
  if (file.exists(done)) { cat("== cell", key, "already complete; loading\n"); return(readRDS(done)) }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE); t0 <- Sys.time()
  cat(sprintf("== cell %s (%d patterns)\n", key, 2^spec$J))
  bank <- make_bank(spec$J, spec$layout, spec$a_const); a <- bank$a; b <- bank$b
  Y <- all_patterns(spec$J)
  sc <- score_all(Y, a, b, cfg); ci <- exact_ci_setup(Y, a, b, cfg$alpha)
  gh <- gauss_hermite(cfg$gh_nodes)
  thetas <- sort(unique(round(c(cfg$theta_grid, gh$x), 10)))
  summ <- do.call(rbind, par_lapply(thetas, summarize_theta, .cores = cfg$cores, Y = Y, a = a, b = b, sc = sc, ci = ci, cfg = cfg))
  summ <- cbind(cell = key, layout = cell_label(spec), J = spec$J, summ)
  summ$on_grid <- round(summ$theta, 10) %in% round(cfg$theta_grid, 10)
  write.csv(summ, file.path(dir, "summary.csv"), row.names = FALSE)
  roots <- do.call(rbind, lapply(names(sc), function(nm) data.frame(cell = key, layout = cell_label(spec), J = spec$J, estimator = nm,
    n_multi = if (all(is.na(sc[[nm]]$nroots))) NA else sum(sc[[nm]]$nroots > 1),
    n_bound = sum(abs(sc[[nm]]$theta) >= cfg$L - 1e-8, na.rm = TRUE), n_na = sum(!is.finite(sc[[nm]]$theta)))))
  write.csv(roots, file.path(dir, "roots.csv"), row.names = FALSE)
  res <- list(summary = summ, roots = roots, gh = list(x = round(gh$x, 10), w = gh$w),
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(spec = spec, bank = bank, scores = sc, ci = ci[c("us", "Ls", "Us")]), file.path(dir, "scores.rds"))
  saveRDS(res, done); cat(sprintf("   done in %.1f min\n", res$minutes)); res
}

## N(0,1)-marginal summaries by Gauss-Hermite quadrature (15 nodes; numerical integration)
marginal_summary <- function(summ, gh) {
  s <- summ[round(summ$theta, 10) %in% gh$x, ]; s$w <- gh$w[match(round(s$theta, 10), gh$x)]
  dcols <- grep("^(wrong3_|undet_|wrong_point_)", names(s), value = TRUE)
  out <- lapply(split(s, list(s$cell, s$estimator, s$se_type, s$eps, s$cond), drop = TRUE), function(k) {
    ww <- k$w / sum(k$w)
    r <- data.frame(cell = k$cell[1], layout = k$layout[1], J = k$J[1], estimator = k$estimator[1], se_type = k$se_type[1],
                    eps = k$eps[1], cond = k$cond[1], rmse = sqrt(sum(ww * k$rmse^2)), abs_bias = sum(ww * abs(k$bias)),
                    coverage = sum(ww * k$coverage), width = sum(ww * k$width), mass_bound = sum(ww * k$mass_bound),
                    mass_unbounded = sum(ww * k$mass_unbounded))
    for (dc in dcols) r[[dc]] <- sum(ww * k[[dc]])
    r
  })
  do.call(rbind, out)
}

## Descriptive convergence summary: slope of log|SE ratio - 1| on log J under the model.
## There is no sampling error: the regression SE measures departure from a power law.
rate_slopes <- function(summ, thetas = c(-2, -1, 0, 1, 2), floor_err = 0.005) {
  m <- summ[summ$eps == 0 & summ$theta %in% thetas & summ$estimator != "EXACT_CI", ]
  m$err <- abs(m$se_ratio - 1)
  out <- lapply(split(m, list(m$layout, m$estimator, m$se_type, m$theta), drop = TRUE), function(k) {
    k <- k[order(k$J), ]
    if (length(unique(k$J)) < 3 || any(!is.finite(k$err)) || any(k$err <= floor_err)) { sl <- NA; se <- NA }
    else { co <- summary(lm(log(k$err) ~ log(k$J)))$coefficients; sl <- co[2, 1]; se <- co[2, 2] }
    data.frame(layout = k$layout[1], estimator = k$estimator[1], se_type = k$se_type[1], theta = k$theta[1], n_J = nrow(k),
               err_minJ = k$err[1], err_maxJ = k$err[nrow(k)], slope = sl, lack_of_fit_se = se)
  })
  do.call(rbind, out)
}
