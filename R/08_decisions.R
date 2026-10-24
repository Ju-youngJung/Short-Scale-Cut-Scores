## =============================================================================
## 08_decisions.R -- raw-score rule, classification accuracy/consistency,
## Rudner's normal approximation, Lord-Wingersky check
## =============================================================================

## Test characteristic curve and TCC-mapped integer raw cuts for a latent cut
tcc <- function(theta, a, b) sum(p2pl(theta, a, b))
raw_cuts_for <- function(cuts, a, b) {
  lapply(cuts, function(ct) { v <- tcc(ct, a, b); c(lo = floor(v), hi = ceiling(v), tcc = v) })
}

## Lord-Wingersky (1984) recursion: exact distribution of the raw score X at theta
lord_wingersky <- function(theta, a, b) {
  P <- p2pl(theta, a, b); f <- 1
  for (j in seq_along(P)) f <- c(f * (1 - P[j]), 0) + c(0, f * P[j])
  f   # f[k + 1] = P(X = k)
}

## Rudner (2001, 2005): P(theta_hat >= cut | theta) under a normal approximation
## with SD equal to the asymptotic SE 1/sqrt(I(theta)) at the TRUE theta.
rudner_p_above <- function(theta, cut, a, b) 1 - pnorm((cut - theta) * sqrt(test_info(theta, a, b)))

## Two-outcome rule summaries from P(above | theta): wrong-decision probability and
## the consistency integrand p^2 + (1-p)^2 (Hambleton & Novick, 1973; Lee, 2010)
two_outcome <- function(th, ct, p_above) c(wrong = if (th >= ct) 1 - p_above else p_above, cons = p_above^2 + (1 - p_above)^2)
