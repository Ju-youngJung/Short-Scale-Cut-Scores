## =============================================================================
## 03_exact_ci.R (v2) -- exact confidence intervals for theta (estimator-free)
## With known items the 2PL is an exponential family in theta with sufficient
## statistic S = sum a_j y_j and monotone likelihood ratio, so inverting the
## equal-tailed exact tests gives [L(s), U(s)]:
##   L(s): P_L(S >= s) = alpha/2 (L = -Inf at the smallest s);  U(s): P_U(S <= s) = alpha/2.
## Coverage >= 1 - alpha for every theta (Klauer, 1991; non-randomized, conservative).
##
## Endpoints: when S has few distinct values (lattice banks, K <= max_root) each endpoint is
## found by root-finding. Otherwise (real banks: S differs for almost every pattern) the
## endpoints are located on a theta grid of step h and taken at the OUTER grid point, so the
## interval is wider than the exact equal-tailed interval by at most h and its coverage
## guarantee is preserved. The method used is recorded in $method.
## =============================================================================

exact_ci_setup <- function(Y, a, b, alpha, digits = 9, max_root = 500L, h = 0.001, lim = 8) {
  S <- round(as.vector(Y %*% a), digits); us <- sort(unique(S)); idx <- match(S, us); K <- length(us)
  lw <- -as.vector(Y %*% (a * b))
  logH <- vapply(split(lw, factor(idx, levels = seq_len(K))), function(v) { m <- max(v); m + log(sum(exp(v - m))) }, numeric(1))
  ps <- function(t) { lp <- logH + t * us; p <- exp(lp - max(lp)); p / sum(p) }
  Ls <- Us <- numeric(K)
  if (K <= max_root) {
    for (k in seq_len(K)) {
      Ls[k] <- if (k == 1) -Inf else uniroot(function(t) sum(ps(t)[k:K]) - alpha / 2, c(-25, 25), tol = 1e-11)$root
      Us[k] <- if (k == K)  Inf else uniroot(function(t) sum(ps(t)[1:k]) - alpha / 2, c(-25, 25), tol = 1e-11)$root
    }
    method <- "root"
  } else {
    grid <- seq(-lim, lim, by = h); Ls[] <- -Inf; Us[] <- Inf
    doneL <- rep(FALSE, K); doneU <- rep(FALSE, K); prevL <- -Inf
    ## scanning t upward: P_t(S >= s_k) increases in t. L(s_k) = first t with tail_ge > alpha/2;
    ## take the previous grid point (outer). U(s_k) = last t with tail_le > alpha/2 ... scanned
    ## downward below, taking the next grid point (outer).
    for (t in grid) { p <- ps(t); tail_ge <- rev(cumsum(rev(p)))
      hit <- !doneL & tail_ge > alpha / 2; Ls[hit] <- prevL; doneL[hit] <- TRUE; prevL <- t
      if (all(doneL)) break }
    Ls[1] <- -Inf
    prevU <- Inf
    for (t in rev(grid)) { p <- ps(t); tail_le <- cumsum(p)
      hit <- !doneU & tail_le > alpha / 2; Us[hit] <- prevU; doneU[hit] <- TRUE; prevU <- t
      if (all(doneU)) break }
    Us[K] <- Inf
    method <- sprintf("grid(h=%g, outer)", h)
  }
  list(S = S, us = us, idx = idx, logH = logH, ps = ps, Ls = Ls, Us = Us, L = Ls[idx], U = Us[idx], method = method)
}
