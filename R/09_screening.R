## =============================================================================
## 09_screening.R (v2) -- person-fit screening with Snijders' (2001) l_z*, calibrated
## by its exact null distribution; vectorized statistic; grid calibration for real banks
##   w_j = log(P_j/Q_j) = a_j(theta - b_j),  c = sum w_j P_j' / sum r_j P_j',  w~_j = w_j - c r_j,
##   l_z* = sum w~_j(theta_hat)(y_j - P_j(theta_hat)) / sqrt(sum w~_j^2 P_j Q_j),  r_j = a_j (WLE).
## VERIFY against Snijders (2001) / Magis, Raiche & Beland (2012).
## =============================================================================

lz_star_wle <- function(Y, theta_hat, a, b) {
  n <- nrow(Y); J <- ncol(Y); A <- matrix(a, n, J, byrow = TRUE)
  U <- sweep(outer(theta_hat, b, "-"), 2, a, "*"); P <- plogis(U); P <- pmin(pmax(P, 1e-15), 1 - 1e-15); Q <- 1 - P
  dP <- A * P * Q
  cc <- rowSums(U * dP) / rowSums(A * dP)                 # c(theta_hat), w_j = U_j
  Wt <- U - cc * A
  num <- rowSums(Wt * (Y - P)); den <- sqrt(rowSums(Wt^2 * P * Q))
  out <- num / den; out[!is.finite(theta_hat) | den <= 0] <- NA_real_; out
}

## exact lower alpha-quantile of a pattern statistic L under probabilities pr
exact_lower_quantile <- function(L, pr, alpha) { o <- order(L); cum <- cumsum(pr[o]); L[o][which(cum >= alpha)[1]] }

## Cutoffs at a set of theta values: q(theta) = exact alpha-quantile of L under the model at theta
screen_cutoffs <- function(theta_vals, Y, L_enum, a, b, alpha) {
  ok <- is.finite(L_enum)
  vapply(theta_vals, function(t) { pr <- exp(pattern_logprob(Y, p2pl(t, a, b))); pr <- pr / sum(pr)
    exact_lower_quantile(L_enum[ok], pr[ok] / sum(pr[ok]), alpha) }, numeric(1))
}

## Plug-in screening of the enumerated patterns. If the WLE takes few distinct values
## (lattice banks) the cutoff is calibrated at each distinct theta_hat; otherwise at the
## grid point (step grid_step) nearest theta_hat. Returns pass indicator and cutoff used.
screen_pass <- function(L, theta_hat, Y, a, b, alpha, max_exact = 1000L, grid_step = 0.02, lim = 4) {
  key <- round(theta_hat, 10); ut <- unique(key[is.finite(key)])
  if (length(ut) <= max_exact) { cal_pts <- ut; assign_idx <- match(key, ut); method <- "per theta_hat" }
  else { cal_pts <- seq(-lim, lim, by = grid_step); assign_idx <- pmin(pmax(round((key + lim) / grid_step) + 1, 1), length(cal_pts))
         method <- sprintf("grid(step=%g)", grid_step) }
  q <- screen_cutoffs(cal_pts, Y, L, a, b, alpha)
  cut <- q[assign_idx]; cut[!is.finite(key)] <- NA_real_
  list(pass = is.finite(L) & is.finite(cut) & L >= cut, cutoff = cut, method = method)
}
