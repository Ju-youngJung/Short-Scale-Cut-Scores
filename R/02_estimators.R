## =============================================================================
## 02_estimators.R -- established ability estimators and their standard errors
##   ML   : Lord (1980)
##   WLE  : Warm (1989)
##   MAP  : Bayes modal, N(mu, sd^2) prior
##   EAP  : Bock & Mislevy (1982), posterior mean and posterior SD by quadrature
##   BW   : biweight estimator (Mislevy & Bock, 1982)
##   HU   : Huber-weighted estimator (Schuster & Yuan, 2011)
## Robust weights use the standardized distance u_j = a_j (theta - b_j).
## VERIFY tuning constants (cfg$bw_c, cfg$hu_k) against the original articles.
## =============================================================================

w_biweight <- function(u, c) ifelse(abs(u) <= c, (1 - (u / c)^2)^2, 0)
w_huber    <- function(u, k) pmin(1, k / pmax(abs(u), 1e-12))

item_weights <- function(name, theta, a, b, cfg) {
  u <- a * (theta - b)
  switch(name, BW = w_biweight(u, cfg$bw_c), HU = w_huber(u, cfg$hu_k), rep(1, length(a)))
}

## At fixed theta every estimating function used here is linear in y:
##   Psi(theta; y) = sum_j c1_j y_j + c0
psi_coef <- function(name, theta, a, b, cfg) {
  q <- item_q(theta, a, b); w <- item_weights(name, theta, a, b, cfg)
  c1 <- w * a; c0 <- -sum(w * a * q$P)
  if (name == "WLE") c0 <- c0 + warm_J(theta, a, b) / (2 * sum(q$s))             # Warm (1989)
  if (name == "MAP") c0 <- c0 - (theta - cfg$prior_mean) / cfg$prior_sd^2
  list(c1 = c1, c0 = c0)
}

## Psi evaluated row-wise with a row-specific theta (vectorized)
psi_rows <- function(name, theta, Y, a, b, cfg) {
  n <- nrow(Y); J <- ncol(Y)
  U <- sweep(outer(theta, b, "-"), 2, a, "*"); P <- plogis(U)
  W <- switch(name, BW = w_biweight(U, cfg$bw_c), HU = w_huber(U, cfg$hu_k), matrix(1, n, J))
  A <- matrix(a, n, J, byrow = TRUE)
  out <- rowSums(W * A * (Y - P))
  if (name == "WLE") { PQ <- P * (1 - P); out <- out + rowSums(A^3 * PQ * (1 - 2 * P)) / (2 * rowSums(A^2 * PQ)) }
  if (name == "MAP") out <- out - (theta - cfg$prior_mean) / cfg$prior_sd^2
  out
}

## Solve Psi = 0 for every pattern.
##  * roots are located by sign changes on a grid over [-L, L] (step cfg$grid_step) and
##    refined by vectorized bisection (40 steps);
##  * no sign change: estimate truncated at +L or -L (convention, e.g. catR range c(-4, 4));
##  * several roots (possible for BW/HU): root nearest the reference estimate (WLE) is
##    reported and the multiplicity is recorded (nroots), so the convention is transparent.
solve_estimator <- function(name, Y, a, b, cfg, ref = NULL, chunk = 4096L) {
  grid <- seq(-cfg$L, cfg$L, by = cfg$grid_step); K <- length(grid); J <- length(a)
  co <- lapply(grid, function(t) psi_coef(name, t, a, b, cfg))
  C1 <- do.call(rbind, lapply(co, `[[`, "c1")); C0 <- vapply(co, `[[`, numeric(1), "c0")
  n <- nrow(Y); theta <- rep(NA_real_, n); nroots <- integer(n); lo <- hi <- rep(NA_real_, n)
  for (st in seq(1L, n, by = chunk)) {
    id <- st:min(n, st + chunk - 1L)
    Fm <- Y[id, , drop = FALSE] %*% t(C1) + matrix(C0, length(id), K, byrow = TRUE)
    Sg <- sign(Fm); ch <- (Sg[, -K, drop = FALSE] * Sg[, -1, drop = FALSE]) < 0
    nroots[id] <- rowSums(ch)
    for (r in seq_along(id)) {
      i <- id[r]; ks <- which(ch[r, ])
      if (!length(ks)) {
        s <- Sg[r, Sg[r, ] != 0]
        theta[i] <- if (!length(s)) NA_real_ else if (all(s > 0)) cfg$L else if (all(s < 0)) -cfg$L else NA_real_
        next
      }
      if (length(ks) > 1L) {
        apx <- grid[ks] - Fm[r, ks] * cfg$grid_step / (Fm[r, ks + 1] - Fm[r, ks])
        k <- ks[which.min(abs(apx - (if (is.null(ref)) 0 else ref[i])))]
      } else k <- ks
      lo[i] <- grid[k]; hi[i] <- grid[k + 1]
    }
  }
  need <- which(!is.na(lo))
  if (length(need)) {
    Yn <- Y[need, , drop = FALSE]; l <- lo[need]; h <- hi[need]
    fl <- psi_rows(name, l, Yn, a, b, cfg)
    for (it in 1:40) {
      m <- (l + h) / 2; fm <- psi_rows(name, m, Yn, a, b, cfg)
      left <- sign(fm) == sign(fl) & fm != 0
      l[left] <- m[left]; fl[left] <- fm[left]; h[!left] <- m[!left]
    }
    theta[need] <- (l + h) / 2
  }
  list(theta = theta, nroots = nroots)
}

## EAP and posterior SD by equally spaced quadrature (Bock & Mislevy, 1982)
eap_scores <- function(Y, a, b, cfg, chunk = 8192L) {
  nodes <- seq(-cfg$eap_lim, cfg$eap_lim, length.out = cfg$eap_nodes)
  lp <- dnorm(nodes, cfg$prior_mean, cfg$prior_sd, log = TRUE)
  Pn <- pmin(pmax(plogis(sweep(outer(nodes, b, "-"), 2, a, "*")), 1e-15), 1 - 1e-15)
  lP <- t(log(Pn)); lQ <- t(log(1 - Pn)); K <- length(nodes)
  n <- nrow(Y); e <- v <- numeric(n)
  for (st in seq(1L, n, by = chunk)) {
    id <- st:min(n, st + chunk - 1L); Yc <- Y[id, , drop = FALSE]
    LL <- Yc %*% lP + (1 - Yc) %*% lQ + matrix(lp, length(id), K, byrow = TRUE)
    W <- exp(LL - apply(LL, 1, max)); W <- W / rowSums(W)
    e[id] <- as.vector(W %*% nodes); v[id] <- as.vector(W %*% nodes^2) - e[id]^2
  }
  list(theta = e, psd = sqrt(pmax(v, 0)))
}

## Standard errors evaluated at the estimate
##  classical : 1/sqrt(I) (ML, WLE); 1/sqrt(I + 1/sd^2) (MAP)
##  sandwich  : M-estimation ASE sqrt(E[Psi^2]) / (-E[Psi']) (Magis, 2014b, 2016)
##              VERIFY the WLE and MAP forms against Magis (2016) before reporting.
##  psd       : EAP posterior SD (Bayesian; frequentist coverage not guaranteed)
se_one <- function(name, type, th, a, b, cfg) {
  if (!is.finite(th)) return(NA_real_)
  q <- item_q(th, a, b); I <- sum(q$s)
  if (type == "classical") return(if (name == "MAP") 1 / sqrt(I + 1 / cfg$prior_sd^2) else 1 / sqrt(I))
  if (name == "WLE") {
    g <- function(t) warm_J(t, a, b) / (2 * test_info(t, a, b))
    gd <- (g(th + 1e-5) - g(th - 1e-5)) / 2e-5
    return(sqrt(I + g(th)^2) / (I - gd))
  }
  if (name == "MAP") return(sqrt(I + ((th - cfg$prior_mean) / cfg$prior_sd^2)^2) / (I + 1 / cfg$prior_sd^2))
  w <- item_weights(name, th, a, b, cfg)
  den <- sum(w * q$s); if (den <= 0) return(NA_real_)
  sqrt(sum(w^2 * q$s)) / den
}

se_vector <- function(name, type, theta, a, b, cfg) {
  key <- round(theta, 10); u <- unique(key[is.finite(key)])
  val <- vapply(u, function(t) se_one(name, type, t, a, b, cfg), numeric(1))
  out <- rep(NA_real_, length(theta)); ok <- is.finite(key); out[ok] <- val[match(key[ok], u)]; out
}

## Score all patterns with the estimators listed in cfg$estimators (WLE is always
## computed: it is the reference used when a robust equation has several roots).
score_all <- function(Y, a, b, cfg) {
  res <- list(); w <- solve_estimator("WLE", Y, a, b, cfg)
  res$WLE <- list(theta = w$theta, nroots = w$nroots)
  for (nm in setdiff(cfg$estimators, "WLE")) {
    if (nm == "EAP") { e <- eap_scores(Y, a, b, cfg); res$EAP <- list(theta = e$theta, nroots = rep(NA_integer_, nrow(Y)), se = list(psd = e$psd)) }
    else { r <- solve_estimator(nm, Y, a, b, cfg, ref = if (nm %in% c("BW", "HU")) w$theta else NULL); res[[nm]] <- list(theta = r$theta, nroots = r$nroots) }
  }
  for (nm in names(res)) for (ty in cfg$se_types[[nm]]) if (ty != "psd")
    res[[nm]]$se[[ty]] <- se_vector(nm, ty, res[[nm]]$theta, a, b, cfg)
  res[intersect(cfg$estimators, names(res))]
}
