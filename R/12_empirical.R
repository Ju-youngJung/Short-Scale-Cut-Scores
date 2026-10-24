## =============================================================================
## 12_empirical.R -- empirical illustration on a real item bank and real respondents
## Data: ICAR cognitive items, psych::ability (16 dichotomous items, N = 1525, public domain)
## The item parameters are estimated once and then treated as known; the exact pipeline
## is applied to that bank. Real data illustrate decision rates for actual respondents;
## they cannot validate coverage (the true theta is unknown).
## =============================================================================

icar_load <- function() {
  ## The ICAR 'ability' data are distributed in psychTools (older psych versions had them in psych)
  find_in <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(NULL)
    env <- new.env()
    ok <- tryCatch({ suppressWarnings(data("ability", package = pkg, envir = env)); exists("ability", envir = env) },
                   error = function(e) FALSE)
    if (ok) get("ability", envir = env) else NULL
  }
  x <- find_in("psychTools"); if (is.null(x)) x <- find_in("psych")
  if (is.null(x)) stop("ICAR 'ability' data not found: run install.packages('psychTools') and retry.")
  dat <- as.data.frame(x); dat[] <- lapply(dat, as.integer); dat
}

## 2PL calibration by marginal ML (mirt); returns a, b (IRT parameterization) and the model
## mirt is attached (not only loaded) so that the S4 methods coef() and residuals() dispatch
attach_mirt <- function() {
  if (!requireNamespace("mirt", quietly = TRUE)) stop("install.packages('mirt')")
  if (!"package:mirt" %in% search()) suppressPackageStartupMessages(library(mirt))
}

calibrate_2pl <- function(dat, SE = TRUE) {
  attach_mirt()
  mod <- mirt::mirt(dat, 1, itemtype = "2PL", SE = SE, verbose = FALSE)
  it <- coef(mod, simplify = TRUE, IRTpars = TRUE)$items   # S4 method from mirt
  list(mod = mod, a = unname(it[, "a"]), b = unname(it[, "b"]), items = rownames(it))
}

fit_summary <- function(mod) {
  attach_mirt(); out <- list()
  out$M2 <- tryCatch(mirt::M2(mod, na.rm = TRUE), error = function(e) NULL)
  out$itemfit <- tryCatch(mirt::itemfit(mod, fit_stats = "S_X2", na.rm = TRUE), error = function(e) NULL)
  out$Q3 <- tryCatch({ q3 <- residuals(mod, type = "Q3", verbose = FALSE); q3 <- as.matrix(q3); q3[lower.tri(q3)] }, error = function(e) NULL)
  out
}

## Exact null cutoffs of l_z* for a set of theta_hat values (used for respondents)
screen_cutoffs <- function(theta_vals, Y, L_enum, a, b, alpha) {
  ok <- is.finite(L_enum)
  vapply(theta_vals, function(t) { pr <- exp(pattern_logprob(Y, p2pl(t, a, b))); pr <- pr / sum(pr)
    exact_lower_quantile(L_enum[ok], pr[ok] / sum(pr[ok]), alpha) }, numeric(1))
}

## Score observed response vectors (complete cases) with every procedure of the design
score_respondents <- function(R, a, b, cfg, Y, ci, L_enum) {
  R <- as.matrix(R); storage.mode(R) <- "integer"; n <- nrow(R); z <- qnorm(1 - cfg$alpha / 2)
  w <- solve_estimator("WLE", R, a, b, cfg); th <- w$theta
  se <- vapply(th, function(t) if (is.finite(t)) 1 / sqrt(test_info(t, a, b)) else NA_real_, numeric(1))
  e <- eap_scores(R, a, b, cfg)
  S <- round(as.vector(R %*% a), 9); idx <- match(S, ci$us)
  L <- ci$Ls[idx]; U <- ci$Us[idx]
  X <- rowSums(R); rawk <- raw_cuts_for(cfg$cuts, a, b)
  lz <- lz_star_wle(R, th, a, b)
  key <- round(th, 10); grid_pts <- seq(-cfg$L, cfg$L, by = 0.02)
  pass <- sapply(cfg$screen_alphas, function(al) {
    q <- screen_cutoffs(grid_pts, Y, L_enum, a, b, al)
    gi <- pmin(pmax(round((key + cfg$L) / 0.02) + 1, 1), length(grid_pts))
    is.finite(lz) & is.finite(key) & lz >= q[gi] })
  pass <- matrix(pass, nrow = n); colnames(pass) <- paste0("pass_s", cfg$screen_alphas)
  out <- data.frame(wle = th, se_wle = se, eap = e$theta, psd = e$psd, raw = X, S = S, exact_L = L, exact_U = U, lz_star = lz, pass)
  dec3 <- function(lo, up, ct) ifelse(lo > ct, "above", ifelse(up < ct, "below", "undetermined"))
  for (i in seq_along(cfg$cuts)) { ct <- cfg$cuts[i]
    out[[paste0("raw_lo_", ct)]]  <- ifelse(X >= rawk[[i]][["lo"]], "above", "below")
    out[[paste0("raw_hi_", ct)]]  <- ifelse(X >= rawk[[i]][["hi"]], "above", "below")
    out[[paste0("wle_point_", ct)]] <- ifelse(th >= ct, "above", "below")
    out[[paste0("wald_", ct)]]  <- dec3(th - z * se, th + z * se, ct)
    out[[paste0("eap_", ct)]]   <- dec3(e$theta - z * e$psd, e$theta + z * e$psd, ct)
    out[[paste0("exact_", ct)]] <- dec3(L, U, ct)
  }
  out
}

## Proportions of decisions among respondents, by rule and cut
respondent_table <- function(sc, cfg) {
  rows <- list()
  for (ct in cfg$cuts) for (rule in c("raw_lo", "raw_hi", "wle_point", "wald", "eap", "exact")) {
    v <- sc[[paste0(rule, "_", ct)]]
    rows[[length(rows) + 1]] <- data.frame(cut = ct, rule = rule, above = mean(v == "above"), below = mean(v == "below"), undetermined = mean(v == "undetermined"))
  }
  tab <- do.call(rbind, rows)
  dis <- do.call(rbind, lapply(cfg$cuts, function(ct) data.frame(cut = ct,
    raw_hi_vs_wle = mean(sc[[paste0("raw_hi_", ct)]] != sc[[paste0("wle_point_", ct)]]),
    raw_lo_vs_wle = mean(sc[[paste0("raw_lo_", ct)]] != sc[[paste0("wle_point_", ct)]]),
    exact_vs_wald_differ = mean(sc[[paste0("exact_", ct)]] != sc[[paste0("wald_", ct)]]))))
  flags <- colMeans(sc[, grep("^pass_s", names(sc)), drop = FALSE])
  list(decisions = tab, disagreement = dis, pass_rates = flags)
}

## Nonparametric bootstrap of the calibration: refit, rerun the exact decision rates (model only,
## GH nodes), and rescore respondents; returns the bootstrap distribution of key quantities
bootstrap_calibration <- function(dat, R_complete, cfg, B = 50, seed = 1) {
  set.seed(seed); n <- nrow(dat); out <- list()
  cfg_b <- cfg; cfg_b$eps_levels <- numeric(0); cfg_b$theta_grid <- numeric(0); cfg_b$screen_alphas <- cfg$screen_alphas[1]
  for (bb in seq_len(B)) {
    cal <- tryCatch(calibrate_2pl(dat[sample.int(n, n, replace = TRUE), ], SE = FALSE), error = function(e) NULL)
    if (is.null(cal) || any(!is.finite(c(cal$a, cal$b))) || any(cal$a <= 0)) next
    a <- cal$a; b <- cal$b; Y <- all_patterns(length(a))
    sc <- score_all(Y, a, b, cfg_b); ci <- exact_ci_setup(Y, a, b, cfg$alpha)
    X <- rowSums(Y); rawk <- raw_cuts_for(cfg$cuts, a, b)
    L <- lz_star_wle(Y, sc$WLE$theta, a, b); pass <- as.matrix(screen_pass(L, sc$WLE$theta, Y, a, b, cfg_b$screen_alphas)$pass)
    gh <- gauss_hermite(cfg$gh_nodes)
    summ <- do.call(rbind, lapply(round(gh$x, 10), summarize_theta_final, Y = Y, a = a, b = b, sc = sc, ci = ci, cfg = cfg_b, Xraw = X, rawk = rawk, pass_mat = pass))
    summ <- cbind(cell = "boot", layout = "boot", J = length(a), summ)
    ms <- marginal_final(summ, list(x = round(gh$x, 10), w = gh$w))
    rs <- score_respondents(R_complete, a, b, cfg, Y, ci, L); rt <- respondent_table(rs, cfg)
    rec <- data.frame(draw = bb)
    for (ct in cfg$cuts) {
      g <- function(est, col) ms[[col]][ms$estimator == est]
      rec[[paste0("wrong_wle_", ct)]] <- g("WLE", paste0("wrong_point_", ct)); rec[[paste0("wrong_rawhi_", ct)]] <- g("RAW_hi", paste0("wrong_point_", ct))
      rec[[paste0("undet_exact_", ct)]] <- g("EXACT_CI", paste0("undet_", ct)); rec[[paste0("undet_wald_", ct)]] <- g("WLE", paste0("undet_", ct))
      rec[[paste0("resp_undet_exact_", ct)]] <- rt$decisions$undetermined[rt$decisions$cut == ct & rt$decisions$rule == "exact"]
      rec[[paste0("resp_rawhi_vs_wle_", ct)]] <- rt$disagreement$raw_hi_vs_wle[rt$disagreement$cut == ct]
    }
    out[[length(out) + 1]] <- rec
    cat(sprintf("  bootstrap draw %d / %d\n", bb, B))
  }
  do.call(rbind, out)
}
