## =============================================================================
## 10_final_summary.R (v2) -- exact performance measures per theta, final design
## Rows per (theta, eps, cond): estimator in {WLE, EAP, HU, EXACT_CI, RAW_lo, RAW_hi, RUDNER}
## Decision columns per cut ct: wrong_point_ct, cons_point_ct (two-outcome rules);
##   wrong3_ct, undet_ct (three-way interval rules).
## Screening columns per screen size a (suffix _sA, e.g. _s0.05): pass_sA, cov_pass_sA (joint),
##   wrong3_pass_sA_ct, undet_pass_sA_ct, wrong_point_pass_sA_ct (joint with passing).
## =============================================================================

summarize_theta_final <- function(th, Y, a, b, sc, ci, cfg, Xraw, rawk, pass_mat) {
  q <- item_q(th, a, b); z <- qnorm(1 - cfg$alpha / 2); wl <- cfg$width_lim
  cuts <- cfg$cuts; alphas <- cfg$screen_alphas; sfx <- paste0("_s", alphas)
  rows <- list()
  base_cols <- c(paste0("wrong_point_", cuts), paste0("cons_point_", cuts), paste0("wrong3_", cuts), paste0("undet_", cuts))
  pass_cols <- unlist(lapply(sfx, function(sx) c(paste0("pass", sx), paste0("cov_pass", sx),
                                                 paste0("wrong3_pass", sx, "_", cuts), paste0("undet_pass", sx, "_", cuts), paste0("wrong_point_pass", sx, "_", cuts))))
  blank <- setNames(rep(NA_real_, length(base_cols) + length(pass_cols)), c(base_cols, pass_cols))

  ## two-outcome rule from an 'above' indicator list (one per cut); joint-with-pass versions
  add_point <- function(d, pr, above) {
    for (i in seq_along(cuts)) {
      ct <- cuts[i]; ab <- above[[i]]; ab[is.na(ab)] <- FALSE
      t2 <- two_outcome(th, ct, sum(pr[ab])); d[paste0("wrong_point_", ct)] <- t2[["wrong"]]; d[paste0("cons_point_", ct)] <- t2[["cons"]]
      for (k in seq_along(alphas)) { ps <- pass_mat[, k]; pp <- sum(pr[ps]); pa <- sum(pr[ab & ps])
        d[paste0("wrong_point_pass", sfx[k], "_", ct)] <- if (th >= ct) pp - pa else pa }
    }
    d
  }
  ## three-way rule from interval limits; joint-with-pass versions
  add_three <- function(d, pr, lower, upper) {
    for (ct in cuts) {
      ab <- lower > ct; be <- upper < ct; ab[is.na(ab)] <- FALSE; be[is.na(be)] <- FALSE
      pa <- sum(pr[ab]); pb <- sum(pr[be]); d[paste0("wrong3_", ct)] <- if (th >= ct) pb else pa; d[paste0("undet_", ct)] <- 1 - pa - pb
      for (k in seq_along(alphas)) { ps <- pass_mat[, k]; pp <- sum(pr[ps]); pa2 <- sum(pr[ab & ps]); pb2 <- sum(pr[be & ps])
        d[paste0("wrong3_pass", sfx[k], "_", ct)] <- if (th >= ct) pb2 else pa2; d[paste0("undet_pass", sfx[k], "_", ct)] <- pp - pa2 - pb2 }
    }
    d
  }
  add_pass <- function(d, pr, covered = NULL) {
    for (k in seq_along(alphas)) { ps <- pass_mat[, k]; d[paste0("pass", sfx[k])] <- sum(pr[ps])
      if (!is.null(covered)) d[paste0("cov_pass", sfx[k])] <- sum(pr[covered & ps]) }
    d
  }
  mk <- function(est, se_type, eps, cond, d, cov = NA, width = NA, unb = NA, bias = NA, sd = NA, rmse = NA, se_ratio = NA, mb = NA, mm = NA)
    data.frame(theta = th, eps = eps, cond = cond, estimator = est, se_type = se_type, bias = bias, sd = sd, rmse = rmse, se_ratio = se_ratio,
               mass_bound = mb, mass_multi = mm, coverage = cov, width = width, mass_unbounded = unb, t(d), check.names = FALSE, stringsAsFactors = FALSE)

  for (eps in c(0, cfg$eps_levels)) for (cond in if (eps == 0) "none" else cfg$conds) {
    cs <- contam_spec(cond, eps, q$P, q$Q)
    pr <- exp(pattern_logprob(Y, (1 - cs$e) * q$P + cs$e * cs$g)); pr <- pr / sum(pr)
    ## exact interval (estimator-free)
    inci <- ci$L <= th & th <= ci$U
    wid <- pmax(0, pmin(ci$U, wl) - pmax(ci$L, -wl)); unb <- !is.finite(ci$L) | !is.finite(ci$U) | ci$L < -wl | ci$U > wl
    d <- add_pass(add_three(blank, pr, ci$L, ci$U), pr, inci)
    rows[[length(rows) + 1]] <- mk("EXACT_CI", "exact", eps, cond, d, cov = sum(pr * inci), width = sum(pr * wid), unb = sum(pr * unb))
    ## raw-score rules (TCC-mapped integer cuts)
    for (side in c("lo", "hi")) {
      d <- add_pass(add_point(blank, pr, lapply(rawk, function(k) Xraw >= k[[side]])), pr)
      rows[[length(rows) + 1]] <- mk(paste0("RAW_", side), "raw", eps, cond, d)
    }
    ## Rudner's normal approximation for the WLE point rule
    d <- blank
    for (ct in cuts) { t2 <- two_outcome(th, ct, rudner_p_above(th, ct, a, b)); d[paste0("wrong_point_", ct)] <- t2[["wrong"]]; d[paste0("cons_point_", ct)] <- t2[["cons"]] }
    rows[[length(rows) + 1]] <- mk("RUDNER", "approx", eps, cond, add_pass(d, pr))
    ## estimators
    for (nm in names(sc)) {
      est <- sc[[nm]]$theta; fin <- is.finite(est); w <- pr[fin] / sum(pr[fin])
      mu <- sum(w * est[fin]); sdv <- sqrt(sum(w * (est[fin] - mu)^2))
      mb <- sum(pr[fin & abs(est) >= cfg$L - 1e-8]); mm <- if (all(is.na(sc[[nm]]$nroots))) NA_real_ else sum(pr[which(sc[[nm]]$nroots > 1)])
      above <- lapply(cuts, function(ct) est >= ct)
      for (ty in names(sc[[nm]]$se)) {
        se <- sc[[nm]]$se[[ty]]; ok <- fin & is.finite(se); lo <- est - z * se; up <- est + z * se
        inw <- ok & lo <= th & th <= up; inw[is.na(inw)] <- FALSE
        d <- add_pass(add_point(add_three(blank, pr, ifelse(ok, lo, -Inf), ifelse(ok, up, Inf)), pr, above), pr, inw)
        rows[[length(rows) + 1]] <- mk(nm, ty, eps, cond, d, cov = sum(pr[inw]), width = sum(pr[ok] * 2 * z * se[ok]) / sum(pr[ok]),
                                       bias = mu - th, sd = sdv, rmse = sqrt((mu - th)^2 + sdv^2), se_ratio = sum(pr[ok] * se[ok]) / sum(pr[ok]) / sdv, mb = mb, mm = mm)
      }
    }
  }
  do.call(rbind, rows)
}

run_cell_final <- function(spec, cfg, out_root) {
  key <- cell_key(spec); dir <- file.path(out_root, key); done <- file.path(dir, "cell.rds")
  if (file.exists(done)) { cat("== cell", key, "already complete; loading\n"); return(readRDS(done)) }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE); t0 <- Sys.time()
  cat(sprintf("== cell %s (%d patterns)\n", key, 2^(if (!is.null(spec$a)) length(spec$a) else spec$J)))
  bank <- if (!is.null(spec$a)) list(a = spec$a, b = spec$b, J = length(spec$a), layout = spec$layout)
          else make_bank(spec$J, spec$layout, spec$a_const)
  a <- bank$a; b <- bank$b; Y <- all_patterns(bank$J)
  sc <- score_all(Y, a, b, cfg); ci <- exact_ci_setup(Y, a, b, cfg$alpha)
  X <- rowSums(Y); rawk <- raw_cuts_for(cfg$cuts, a, b)
  L <- lz_star_wle(Y, sc$WLE$theta, a, b)
  scr <- lapply(cfg$screen_alphas, function(al) screen_pass(L, sc$WLE$theta, Y, a, b, al))
  pass_mat <- do.call(cbind, lapply(scr, `[[`, "pass")); colnames(pass_mat) <- paste0("s", cfg$screen_alphas)
  gh <- gauss_hermite(cfg$gh_nodes); thetas <- sort(unique(round(c(cfg$theta_grid, gh$x), 10)))
  summ <- do.call(rbind, par_lapply(thetas, summarize_theta_final, .cores = cfg$cores, Y = Y, a = a, b = b, sc = sc, ci = ci, cfg = cfg,
                                    Xraw = X, rawk = rawk, pass_mat = pass_mat))
  summ <- cbind(cell = key, layout = cell_label(spec), J = bank$J, summ)
  summ$on_grid <- round(summ$theta, 10) %in% round(cfg$theta_grid, 10)
  write.csv(summ, file.path(dir, "summary.csv"), row.names = FALSE)
  lw <- max(vapply(c(-2, 0, 1.5), function(t) { pr <- exp(pattern_logprob(Y, p2pl(t, a, b))); f <- tapply(pr, factor(X, levels = 0:bank$J), sum)
    max(abs(as.numeric(f) - lord_wingersky(t, a, b))) }, numeric(1)))
  size <- do.call(rbind, lapply(seq_along(cfg$screen_alphas), function(k) data.frame(cell = key, layout = cell_label(spec), J = bank$J,
    screen_alpha = cfg$screen_alphas[k], theta = cfg$theta_grid,
    size = vapply(cfg$theta_grid, function(t) { pr <- exp(pattern_logprob(Y, p2pl(t, a, b))); sum(pr[!pass_mat[, k]]) / sum(pr) }, numeric(1)))))
  info <- data.frame(cell = key, layout = cell_label(spec), J = bank$J, lw_max_diff = lw,
                     raw_cuts = paste(sapply(rawk, function(k) sprintf("%.2f->[%.0f,%.0f]", k["tcc"], k["lo"], k["hi"])), collapse = "; "),
                     n_multi_HU = if (!is.null(sc$HU)) sum(sc$HU$nroots > 1) else NA)
  res <- list(summary = summ, info = info, screen_size = size, gh = list(x = round(gh$x, 10), w = gh$w),
              minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(spec = spec, bank = bank, scores = sc, lz = L, pass = pass_mat, cutoffs = lapply(scr, `[[`, "cutoff")), file.path(dir, "scores.rds"))
  saveRDS(res, done); cat(sprintf("   done in %.1f min (LW check %.1e)\n", res$minutes, lw)); res
}

## Marginal summary: GH averages of all numeric columns (joint quantities stay joint)
marginal_final <- function(summ, gh) {
  s <- summ[round(summ$theta, 10) %in% gh$x, ]; s$w <- gh$w[match(round(s$theta, 10), gh$x)]
  num <- setdiff(names(s)[vapply(s, is.numeric, logical(1))], c("theta", "eps", "J", "w"))
  out <- lapply(split(s, list(s$cell, s$estimator, s$se_type, s$eps, s$cond), drop = TRUE), function(k) {
    ww <- k$w / sum(k$w)
    r <- data.frame(cell = k$cell[1], layout = k$layout[1], J = k$J[1], estimator = k$estimator[1], se_type = k$se_type[1], eps = k$eps[1], cond = k$cond[1])
    for (v in num) r[[v]] <- if (v == "rmse") sqrt(sum(ww * k[[v]]^2)) else sum(ww * k[[v]])
    r
  })
  do.call(rbind, out)
}
