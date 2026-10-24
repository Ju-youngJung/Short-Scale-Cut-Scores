## =============================================================================
## 05_checks.R -- implementation checks based on established results
## Each returns list(pass, value, note); pass = NA means "reported only".
## =============================================================================

## (1) Exponential-family representation: P_t(S = s) from H(s) equals enumeration
check_expfam <- function(Y, a, b, ci, thetas = c(-2, 0, 1.3)) {
  worst <- 0
  for (t in thetas) {
    pr <- exp(pattern_logprob(Y, p2pl(t, a, b)))
    en <- as.vector(tapply(pr, factor(ci$idx, levels = seq_along(ci$us)), sum))
    worst <- max(worst, max(abs(en - ci$ps(t))))
  }
  list(pass = worst < 1e-10, value = worst, note = "P(S=s) via sufficient statistic = enumeration")
}

## (2) Monotone likelihood ratio: P_t(S >= s) increasing in t (needed for interval form)
check_mlr <- function(ci, grid = seq(-6, 6, by = 0.05)) {
  K <- length(ci$us); worst <- 0
  for (k in 2:K) { v <- vapply(grid, function(t) sum(ci$ps(t)[k:K]), numeric(1)); worst <- max(worst, -min(diff(v))) }
  list(pass = worst < 1e-12, value = worst, note = "P_t(S >= s) non-decreasing in t for all s")
}

## (3) Exact CI coverage >= 1 - alpha under the model (Neyman construction)
check_ci_coverage <- function(Y, a, b, ci, alpha, thetas = seq(-3, 3, by = 0.25)) {
  cv <- vapply(thetas, function(t) { pr <- exp(pattern_logprob(Y, p2pl(t, a, b))); sum(pr[ci$L <= t & t <= ci$U]) }, numeric(1))
  list(pass = min(cv) >= 1 - alpha - 1e-10, value = min(cv), note = "min exact-CI coverage under the model (>= 1 - alpha)")
}

## (4) Robust estimating functions: E[Psi] = 0 and Var[Psi] = sum w^2 a^2 P Q over patterns
check_robust_moments <- function(Y, a, b, cfg, thetas = c(-2, 0, 1.5)) {
  worst <- 0
  for (nm in c("BW", "HU")) for (t in thetas) {
    co <- psi_coef(nm, t, a, b, cfg); Tv <- as.vector(Y %*% co$c1) + co$c0
    pr <- exp(pattern_logprob(Y, p2pl(t, a, b))); q <- item_q(t, a, b); w <- item_weights(nm, t, a, b, cfg)
    E <- sum(pr * Tv); V <- sum(pr * Tv^2) - E^2
    worst <- max(worst, abs(E), abs(V - sum(w^2 * q$s)))
  }
  list(pass = worst < 1e-9, value = worst, note = "robust Psi: E = 0 and Var = sum w^2 s (enumeration)")
}

## (5) Solver check: Psi(theta_hat) = 0 for interior estimates of the scored estimators
check_roots <- function(Y, a, b, cfg, sc) {
  worst <- 0
  for (nm in intersect(c("ML", "WLE", "MAP", "BW", "HU"), names(sc))) {
    th <- sc[[nm]]$theta; ok <- is.finite(th) & abs(th) < cfg$L - 1e-6
    if (any(ok)) worst <- max(worst, max(abs(psi_rows(nm, th[ok], Y[ok, , drop = FALSE], a, b, cfg))))
  }
  list(pass = worst < 1e-7, value = worst, note = "max |Psi(theta_hat)| at interior estimates")
}

## (6) External: WLE and EAP against mirt (optional)
check_mirt <- function(Y, a, b, cfg, sc, tol = 1e-4) {
  if (!requireNamespace("mirt", quietly = TRUE)) return(list(pass = NA, value = NA, note = "mirt not installed; skipped"))
  out <- tryCatch({
    dat <- data.frame(Y); names(dat) <- paste0("i", seq_len(ncol(Y)))
    pars <- mirt::mirt(dat, 1, itemtype = "2PL", pars = "values", verbose = FALSE)
    pars$value[pars$name == "a1"] <- a; pars$value[pars$name == "d"] <- -a * b; pars$est <- FALSE
    mod <- mirt::mirt(dat, 1, itemtype = "2PL", pars = pars, verbose = FALSE)
    wle <- as.vector(mirt::fscores(mod, method = "WLE"))
    eap <- as.vector(mirt::fscores(mod, method = "EAP", quadpts = cfg$eap_nodes, theta_lim = c(-cfg$eap_lim, cfg$eap_lim)))
    ok <- is.finite(wle) & abs(sc$WLE$theta) < cfg$L - 1e-6
    c(max(abs(wle[ok] - sc$WLE$theta[ok])), max(abs(eap - sc$EAP$theta)))
  }, error = function(e) c(NA, NA))
  list(pass = if (any(is.na(out))) NA else all(out < tol), value = out, note = "max |mirt - own| for WLE, EAP")
}

## (7) Summary integrity: no missing values where a value is expected
check_summary_integrity <- function(summ) {
  est <- summ[summ$estimator != "EXACT_CI", ]
  dc <- grep("^(wrong3_|undet_|wrong_point_)", names(summ), value = TRUE)
  bad <- sum(is.na(est[, c("bias", "sd", "rmse", "coverage", "width", dc)])) + sum(is.na(summ[summ$estimator == "EXACT_CI", c("coverage", "width", grep("^(wrong3_|undet_)", names(summ), value = TRUE))]))
  list(pass = bad == 0, value = bad, note = "number of unexpected NA cells in the summary")
}

run_checks <- function(checks) {
  for (nm in names(checks)) { r <- checks[[nm]]
    st <- if (is.na(r$pass)) "REPORT" else if (r$pass) "PASS" else "FAIL"
    cat(sprintf("[%-6s] %-18s %s | value = %s\n", st, nm, r$note, paste(format(r$value, digits = 4), collapse = ", "))) }
  invisible(checks)
}
