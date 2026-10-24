## =============================================================================
## 01_model.R -- 2PL primitives, item banks, enumeration, contamination, quadrature
## =============================================================================

p2pl <- function(theta, a, b) plogis(a * (theta - b))

item_q <- function(theta, a, b) {
  P <- pmin(pmax(p2pl(theta, a, b), 1e-15), 1 - 1e-15); Q <- 1 - P
  list(P = P, Q = Q, s = a^2 * P * Q)
}
test_info <- function(theta, a, b) sum(item_q(theta, a, b)$s)
## Warm's J(theta) = sum a^3 P Q (Q - P)
warm_J <- function(theta, a, b) { q <- item_q(theta, a, b); sum(a^3 * q$P * q$Q * (q$Q - q$P)) }

## Structured banks: wide b on [-2,2], peaked b on [-0.75,0.75]; a cycles 0.8/1.2/1.6/2.0
make_bank <- function(J, layout = c("wide", "peaked"), a_const = NULL) {
  layout <- match.arg(layout)
  b <- switch(layout, wide = seq(-2, 2, length.out = J), peaked = seq(-0.75, 0.75, length.out = J))
  a <- if (is.null(a_const)) rep(c(0.8, 1.2, 1.6, 2.0), length.out = J) else rep(a_const, J)
  list(a = a, b = b, J = J, layout = layout)
}

## All 2^J response patterns (rows)
all_patterns <- function(J) {
  if (J > 20) stop("Exact enumeration is intended for J <= 20.")
  n <- 2^J; Y <- matrix(0L, n, J); idx <- 0:(n - 1)
  for (j in seq_len(J)) Y[, j] <- as.integer(bitwAnd(idx, as.integer(2^(j - 1))) > 0)
  Y
}

## log P(y | p) for every pattern, p = item success probabilities
pattern_logprob <- function(Y, p) {
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  as.vector(Y %*% log(p) + (1 - Y) %*% log(1 - p))
}

## Contamination (evaluation conditions only): P_j^eps = (1 - e_j) P_j + e_j g_j
contam_spec <- function(cond, eps, P, Q) {
  J <- length(P)
  switch(cond,
         none  = list(e = rep(0, J),   g = rep(0, J)),
         up    = list(e = rep(eps, J), g = rep(1, J)),      # all contaminating responses correct
         down  = list(e = rep(eps, J), g = rep(0, J)),      # all contaminating responses incorrect
         guess = list(e = eps * Q,     g = rep(0.25, J)),   # guessing (4 options) more likely on hard items
         slip  = list(e = eps * P,     g = rep(0, J)),      # slips more likely on easy items
         stop("unknown condition: ", cond))
}

## Gauss-Hermite nodes/weights for N(0,1) (Golub & Welsch, 1969)
gauss_hermite <- function(n) {
  i <- seq_len(n - 1); Jm <- matrix(0, n, n)
  Jm[cbind(i, i + 1)] <- sqrt(i); Jm[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(Jm, symmetric = TRUE)
  o <- order(e$values); list(x = e$values[o], w = (e$vectors[1, ]^2)[o])
}

## Parallel lapply: fork-based on Unix, socket cluster on Windows (same results either way).
## Formals are dotted so that arguments forwarded through ... cannot collide with them.
par_lapply <- function(.X, .FUN, .cores = 1L, ...) {
  if (.cores <= 1L) return(lapply(.X, .FUN, ...))
  if (.Platform$OS.type != "windows") return(parallel::mclapply(.X, .FUN, ..., mc.cores = .cores))
  cl <- parallel::makeCluster(.cores); on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, ls(globalenv()), envir = globalenv())
  parallel::parLapply(cl, .X, .FUN, ...)
}
