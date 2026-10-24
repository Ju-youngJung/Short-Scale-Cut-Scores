## =============================================================================
## make_supplement_xlsx_v2.R -- supplementary materials workbook
##
## Every number quoted in the article can be checked in this workbook. Sheet
## "A_Article_Values" lists each quoted value, the value recomputed here from the
## output files, the sheet that contains it, and the filter that reproduces it.
##
## Inputs are located automatically anywhere under `root` (see the console log).
## Usage: setwd(<project>); source("scripts/make_supplement_xlsx.R")   (or via scripts/run_everything.R)
## =============================================================================

for (p in c("openxlsx")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(openxlsx)

## ---- settings ------------------------------------------------------------------
if (!exists("root"))       root       <- getwd()                       # project folder (contains results/)
if (!exists("out_file"))   out_file   <- file.path(root, "results", "Supplementary_Materials.xlsx")
if (!exists("refit_icar")) refit_icar <- TRUE   # recompute ICAR fit statistics (M2, S-X2, Q3) with mirt + psychTools

## ---- locate inputs -----------------------------------------------------------------
find_all <- function(pattern) list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
read1 <- function(p) {
  d <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
  if ("on_grid" %in% names(d)) d$on_grid <- as.logical(d$on_grid)
  d
}
classify <- function(p) {
  d <- tryCatch(read.csv(p, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !("layout" %in% names(d))) return("other")
  if (any(d$layout == "ICAR")) return("icar")
  if (setequal(unique(d$J), c(8, 12, 16, 20)) && all(c("pass_s0.05", "pass_s0.2") %in% names(d))) return("final")
  "other"
}
pick_file <- function(pattern, kind = NULL) {
  f <- find_all(pattern)
  f <- f[!grepl("/(v1_backup|old|archive)/", f, ignore.case = TRUE)]
  if (!is.null(kind) && length(f)) f <- f[vapply(f, classify, character(1)) == kind]
  if (!length(f)) { message("NOT FOUND (sheets skipped): ", pattern, if (!is.null(kind)) paste0(" [", kind, "]")); return(NULL) }
  f <- f[order(file.mtime(f), decreasing = TRUE)]
  message("  ", sprintf("%-38s", paste0(pattern, if (!is.null(kind)) paste0(" [", kind, "]"))), " -> ", f[1])
  f[1]
}
rd <- function(p) if (is.null(p)) NULL else read1(p)

cat("Locating input files under", root, "\n")
s     <- rd(pick_file("^ALL_summary.*\\.csv$", "final"))
ms    <- rd(pick_file("^ALL_marginal_summary.*\\.csv$", "final"))
ssz   <- rd(pick_file("^ALL_screen_size.*\\.csv$"))
cinfo <- rd(pick_file("^ALL_cell_info.*\\.csv$"))
ibank <- rd(pick_file("^icar_bank_2pl.*\\.csv$"))
is    <- rd(pick_file("^ALL_summary.*\\.csv$", "icar"))
ims   <- rd(pick_file("^ALL_marginal_summary.*\\.csv$", "icar"))
irs   <- rd(pick_file("^respondent_scores.*\\.csv$"))
iboot <- rd(pick_file("^bootstrap_calibration.*\\.csv$"))
logs  <- find_all("^run_log.*\\.txt$")
if (!is.null(ssz))   ssz   <- ssz[ssz$layout %in% c("wide", "peaked") & ssz$J %in% c(8, 12, 16, 20), ]
if (!is.null(cinfo)) cinfo <- cinfo[cinfo$layout %in% c("wide", "peaked") & cinfo$J %in% c(8, 12, 16, 20), ]

## ---- labels and filtering ---------------------------------------------------------------
proc_lab <- c(WLE = "WLE (Wald, 1/sqrt(I))", EAP = "EAP (posterior SD)", HU = "Huber (M-estimation SE)",
              EXACT_CI = "Exact interval", RAW_lo = "Raw score (k-)", RAW_hi = "Raw score (k+)",
              RUDNER = "Rudner approximation (WLE)")
cond_lab <- c(none = "Model", up = "All-correct", down = "All-incorrect", guess = "Guessing")
prep <- function(d) {
  if (is.null(d)) return(NULL)
  d <- d[d$estimator %in% names(proc_lab), ]
  d <- d[!(d$estimator == "WLE" & d$se_type != "classical"), ]
  d$Procedure <- unname(proc_lab[d$estimator]); d$Condition <- unname(cond_lab[d$cond]); d
}
s <- prep(s); ms <- prep(ms); is <- prep(is); ims <- prep(ims)
g  <- if (!is.null(s))  s[s$on_grid, ]  else NULL
ig <- if (!is.null(is)) is[is$on_grid, ] else NULL

sel <- function(d, ...) { a <- list(...); k <- rep(TRUE, nrow(d)); for (n in names(a)) k <- k & d[[n]] %in% a[[n]]; d[k, , drop = FALSE] }
v1  <- function(d, col, ...) { x <- sel(d, ...)[[col]]; if (length(x) == 1) x else NA_real_ }
cuts_of  <- function(d) sort(as.numeric(sub("^undet_", "", grep("^undet_[0-9.]+$", names(d), value = TRUE))))
sizes_of <- function(d) sub("^pass_s", "", grep("^pass_s[0-9.]+$", names(d), value = TRUE))
pick <- function(d, cols) d[, intersect(cols, names(d)), drop = FALSE]
ord  <- function(d, keys) d[do.call(order, d[intersect(keys, names(d))]), , drop = FALSE]
dec_cols <- function(d) grep("^(wrong_point|wrong3|undet)_[0-9.]+$", names(d), value = TRUE)

two_long <- function(d, keys) {
  d <- d[d$estimator %in% c("RAW_lo", "RAW_hi", "WLE", "EAP", "RUDNER"), ]
  do.call(rbind, lapply(cuts_of(d), function(ct) { w <- d[[paste0("wrong_point_", ct)]]
    data.frame(d[, keys, drop = FALSE], Cut = ct, Wrong = w, Accuracy = 1 - w, Consistency = d[[paste0("cons_point_", ct)]], check.names = FALSE) }))
}
three_long <- function(d, keys) {
  d <- d[d$estimator %in% c("WLE", "EAP", "EXACT_CI"), ]
  do.call(rbind, lapply(cuts_of(d), function(ct)
    data.frame(d[, keys, drop = FALSE], Cut = ct, Wrong = d[[paste0("wrong3_", ct)]], Undetermined = d[[paste0("undet_", ct)]], check.names = FALSE)))
}
add_ratios_to <- function(d, ref = "WLE", keys = c("layout", "J", "eps", "cond")) {   # ratio to WLE (width, rmse)
  r <- d[d$estimator == ref, c(keys, "width", "rmse")]; names(r)[-seq_along(keys)] <- c("width_WLE", "rmse_WLE")
  x <- merge(d, r, by = keys, all.x = TRUE)
  x$width_ratio_to_WLE <- x$width / x$width_WLE; x$rmse_ratio_to_WLE <- x$rmse / x$rmse_WLE
  x[, setdiff(names(x), c("width_WLE", "rmse_WLE"))]
}

## ---- workbook ---------------------------------------------------------------------------
wb <- createWorkbook()
st_title <- createStyle(textDecoration = "bold", fontSize = 12)
st_note  <- createStyle(textDecoration = "italic")
st_head  <- createStyle(textDecoration = "bold", border = "Bottom", fgFill = "#E7E6E6", halign = "center", wrapText = TRUE)
st_num   <- createStyle(numFmt = "0.000")
index <- data.frame(Sheet = character(0), Content = character(0), stringsAsFactors = FALSE)
add_sheet <- function(name, d, title, note = NULL, int_cols = c("J", "theta", "eps", "Cut", "n", "draw", "Count", "N")) {
  if (is.null(d) || !nrow(d)) { message("  skipped empty sheet: ", name); return(invisible()) }
  addWorksheet(wb, name)
  writeData(wb, name, title, startRow = 1); addStyle(wb, name, st_title, rows = 1, cols = 1)
  if (!is.null(note)) { writeData(wb, name, note, startRow = 2); addStyle(wb, name, st_note, rows = 2, cols = 1) }
  writeData(wb, name, d, startRow = 4, headerStyle = st_head)
  num <- which(vapply(d, is.numeric, logical(1)) & !(names(d) %in% int_cols))
  if (length(num)) addStyle(wb, name, st_num, rows = 5:(nrow(d) + 4), cols = num, gridExpand = TRUE)
  freezePane(wb, name, firstActiveRow = 5)
  setColWidths(wb, name, cols = seq_along(d), widths = "auto")
  index[nrow(index) + 1, ] <<- c(name, title)
}
keysC <- c("layout", "J", "theta", "Condition", "eps", "Procedure")
keysM <- c("layout", "J", "Condition", "eps", "Procedure")

## ---- V: verification ----------------------------------------------------------------------
if (length(logs)) {
  rx <- "^\\[(PASS|FAIL|REPORT)\\s*\\]\\s+(\\S+)\\s+(.*?)\\s*\\|\\s*value\\s*=\\s*(.*)$"
  ver <- do.call(rbind, lapply(logs, function(f) {
    L <- grep(rx, readLines(f, warn = FALSE), value = TRUE, perl = TRUE)
    if (!length(L)) return(NULL)
    m <- regmatches(L, regexec(rx, L, perl = TRUE))
    data.frame(Log = sub(paste0("^", root, "/?"), "", f), Status = sapply(m, `[`, 2), Check = sapply(m, `[`, 3),
               Description = sapply(m, `[`, 4), Value = sapply(m, `[`, 5), stringsAsFactors = FALSE)
  }))
  add_sheet("V_Verification", ver, "Implementation checks (copied from the run logs)",
            "Each check is an identity that must hold to numerical precision; a failure would indicate an implementation error.")
}
if (!is.null(cinfo)) add_sheet("S1_Cells", cinfo, "Cells of the exact evaluation study",
  "lw_max_diff = maximum difference between the enumerated and Lord-Wingersky raw-score distributions; raw_cuts = TCC(cut) -> [k-, k+].")

## ---- S2-S11: exact study --------------------------------------------------------------------
if (!is.null(ms)) {
  m0 <- add_ratios_to(ms[ms$eps == 0 & ms$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI"), ])
  add_sheet("S2_Interval_Marginal", ord(pick(m0, c("layout", "J", "Procedure", "coverage", "width", "rmse", "width_ratio_to_WLE", "rmse_ratio_to_WLE", "mass_unbounded", "mass_bound")), c("layout", "J", "Procedure")),
            "Marginal coverage, expected width, and RMSE of interval procedures under the model",
            "Marginal = average over theta ~ N(0,1) by 15-node Gauss-Hermite quadrature. Ratios are relative to the WLE Wald interval / WLE.")
  t4 <- two_long(ms[ms$eps == 0, ], c("layout", "J", "Procedure", "estimator"))
  wl <- t4[t4$estimator == "WLE", c("layout", "J", "Cut", "Wrong")]; names(wl)[4] <- "Wrong_WLE"
  t4 <- merge(t4, wl, by = c("layout", "J", "Cut")); t4$Diff_from_WLE <- t4$Wrong - t4$Wrong_WLE
  t4$Relative_to_WLE <- t4$Wrong / t4$Wrong_WLE - 1; t4$estimator <- NULL; t4$Wrong_WLE <- NULL
  add_sheet("S4_TwoOutcome_Marginal", ord(t4, c("layout", "J", "Cut", "Procedure")),
            "Marginal wrong-decision rate, accuracy, and consistency of two-outcome rules under the model",
            "Accuracy = 1 - wrong. Diff_from_WLE and Relative_to_WLE compare each rule with the WLE point rule (Relative = ratio - 1).")
  add_sheet("S6_ThreeWay_Marginal", ord(three_long(ms[ms$eps == 0, ], c("layout", "J", "Procedure")), c("layout", "J", "Cut", "Procedure")),
            "Marginal wrong and undetermined probabilities of three-way interval rules under the model",
            "Three-way rule: above if lower limit > cut, below if upper limit < cut, otherwise undetermined.")
  mc <- add_ratios_to(ms[ms$eps > 0 & ms$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI", "RAW_hi", "RAW_lo"), ])
  cov0 <- ms[ms$eps == 0, c("layout", "J", "estimator", "coverage")]; names(cov0)[4] <- "coverage_model"
  mc <- merge(mc, cov0, by = c("layout", "J", "estimator"), all.x = TRUE); mc$coverage_loss <- mc$coverage_model - mc$coverage
  add_sheet("S8_Contam_Marginal", ord(pick(mc, c(keysM, "coverage", "coverage_model", "coverage_loss", "width", "rmse", "rmse_ratio_to_WLE", dec_cols(mc))), keysM),
            "Marginal coverage, width, RMSE, and decision rates under contamination",
            "Contamination: P_eps = (1-e_j) P_j + e_j g_j. All-correct g = 1; all-incorrect g = 0; guessing e_j = eps*Q_j, g = 1/4.")
  if (length(sizes_of(ms))) {
    sc <- ms[ms$estimator %in% c("WLE", "EXACT_CI"), ]
    for (sz in sizes_of(sc)) { sx <- paste0("_s", sz)
      p0 <- sc[sc$eps == 0, c("layout", "J", "estimator", paste0("pass", sx))]; names(p0)[4] <- paste0("pass_model", sx)
      sc <- merge(sc, p0, by = c("layout", "J", "estimator"), all.x = TRUE)
      sc[[paste0("implied_detection", sx)]] <- 1 - sc[[paste0("pass", sx)]] / sc[[paste0("pass_model", sx)]]
      sc[[paste0("cov_given_pass", sx)]] <- sc[[paste0("cov_pass", sx)]] / sc[[paste0("pass", sx)]]
      sc[[paste0("coverage_gain", sx)]] <- sc[[paste0("cov_given_pass", sx)]] - sc$coverage
      for (ct in cuts_of(sc)) {
        sc[[paste0("wrong3_given_pass", sx, "_", ct)]] <- sc[[paste0("wrong3_pass", sx, "_", ct)]] / sc[[paste0("pass", sx)]]
        sc[[paste0("undet_given_pass", sx, "_", ct)]]  <- sc[[paste0("undet_pass", sx, "_", ct)]]  / sc[[paste0("pass", sx)]] } }
    cols <- c(keysM, "coverage", grep("^(pass_s|pass_model_s|implied_detection_s|cov_given_pass_s|coverage_gain_s)", names(sc), value = TRUE),
              grep("^(wrong3|undet)_[0-9.]+$|given_pass_s[0-9.]+_", names(sc), value = TRUE))
    add_sheet("S10_Screening", ord(pick(sc, unique(cols)), keysM),
              "Person-fit screening with l_z* (exact-null plug-in cutoff): pass rates, coverage, and decisions",
              "implied_detection = 1 - pass / pass under the model; *_given_pass = value among examinees who pass; coverage_gain = cov_given_pass - coverage.")
  }
}
if (!is.null(g)) {
  add_sheet("S3_Interval_Conditional", ord(pick(g[g$eps == 0 & g$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI"), ], c("layout", "J", "theta", "Procedure", "bias", "sd", "rmse", "se_ratio", "coverage", "width", "mass_unbounded", "mass_bound", "mass_multi")), c("layout", "J", "Procedure", "theta")),
            "Conditional bias, SD, RMSE, SE ratio, coverage, and width under the model",
            "se_ratio = mean reported SE / exact conditional SD; mass_bound = probability of a truncated estimate; mass_multi = probability of a multi-rooted Huber equation.")
  add_sheet("S5_TwoOutcome_Conditional", ord(two_long(g[g$eps == 0, ], c("layout", "J", "theta", "Procedure")), c("layout", "J", "Cut", "Procedure", "theta")),
            "Conditional wrong-decision probability of two-outcome rules under the model",
            "Wrong = P(decide below | theta) if theta >= cut, P(decide above | theta) otherwise.")
  add_sheet("S7_ThreeWay_Conditional", ord(three_long(g[g$eps == 0, ], c("layout", "J", "theta", "Procedure")), c("layout", "J", "Cut", "Procedure", "theta")),
            "Conditional wrong and undetermined probabilities of three-way interval rules under the model")
  cc <- g[g$eps > 0 & g$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI"), ]
  add_sheet("S9_Contam_Conditional", ord(pick(cc, c(keysC, "bias", "rmse", "coverage", "width", "pass_s0.05", "cov_pass_s0.05", "pass_s0.2", "cov_pass_s0.2", dec_cols(cc))), c("layout", "J", "Condition", "eps", "Procedure", "theta")),
            "Conditional coverage, decision rates, and screening joint probabilities under contamination",
            "cov_pass_sA / pass_sA = coverage among examinees who pass the screen of size A.")
}
if (!is.null(ssz)) add_sheet("S11_Screen_Size", ord(pick(ssz, c("layout", "J", "screen_alpha", "theta", "size")), c("layout", "J", "screen_alpha", "theta")),
  "Exact conditional size of the screens under the model", "size = P(fail | theta) with the cutoff calibrated at the examinee's WLE.")

## ---- S12-S19: empirical study -------------------------------------------------------------
if (!is.null(ibank)) {
  tcc <- function(t) sum(plogis(ibank$a * (t - ibank$b)))
  add_sheet("S12_ICAR_Items", ibank, "ICAR item parameters (2PL, marginal maximum likelihood)",
            sprintf("TCC(0.67) = %.2f; TCC(1.28) = %.2f. Parameters treated as known in the exact computations.", tcc(0.67), tcc(1.28)))
}
icar_fit <- NULL
if (refit_icar && requireNamespace("mirt", quietly = TRUE) && requireNamespace("psychTools", quietly = TRUE)) {
  icar_fit <- tryCatch({
    suppressPackageStartupMessages(library(mirt))
    env <- new.env(); data("ability", package = "psychTools", envir = env)
    dat <- as.data.frame(get("ability", envir = env)); dat[] <- lapply(dat, as.integer)
    mod <- mirt(dat, 1, itemtype = "2PL", verbose = FALSE)
    m2 <- M2(mod, na.rm = TRUE); itf <- itemfit(mod, fit_stats = "S_X2", na.rm = TRUE)
    q3 <- as.matrix(residuals(mod, type = "Q3", verbose = FALSE)); q3 <- q3[lower.tri(q3)]
    list(m2 = data.frame(Statistic = names(m2), Value = as.numeric(m2[1, ])), itemfit = itf,
         q3 = data.frame(Statistic = c("max |Q3|", "pairs with |Q3| > 0.20", "pairs"), Value = c(max(abs(q3), na.rm = TRUE), sum(abs(q3) > 0.2, na.rm = TRUE), sum(is.finite(q3)))))
  }, error = function(e) { message("  ICAR refit failed: ", conditionMessage(e)); NULL })
}
if (!is.null(icar_fit)) {
  add_sheet("S13_ICAR_Fit", rbind(icar_fit$m2, icar_fit$q3), "ICAR 2PL model fit (M2 limited-information statistics and Q3 local dependence)",
            "M2 computed on complete cases (mirt removes incomplete rows for M2).")
  add_sheet("S14_ICAR_ItemFit", icar_fit$itemfit, "ICAR item fit (S-X2)")
}
if (!is.null(ims)) {
  i0 <- ims[ims$eps == 0, ]
  add_sheet("S15_ICAR_Marginal", rbind(
      transform(two_long(i0, "Procedure"), Rule_type = "two-outcome", Undetermined = NA),
      transform(three_long(i0, "Procedure"), Rule_type = "three-way", Accuracy = NA, Consistency = NA)[, c("Procedure", "Cut", "Wrong", "Accuracy", "Consistency", "Rule_type", "Undetermined")]),
    "ICAR bank: exact marginal decision rates under the model")
  add_sheet("S16_ICAR_Intervals", ord(pick(add_ratios_to(ims[ims$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI"), ]), c("Condition", "eps", "Procedure", "coverage", "width", "rmse", "width_ratio_to_WLE")), c("eps", "Condition", "Procedure")),
            "ICAR bank: exact marginal coverage, width, and RMSE under the model and contamination")
}
if (!is.null(ig)) add_sheet("S17_ICAR_Conditional", ord(pick(ig[ig$estimator %in% names(proc_lab), ], c(keysC, "bias", "sd", "se_ratio", "coverage", "width", dec_cols(ig), "pass_s0.05", "cov_pass_s0.05", "pass_s0.2", "cov_pass_s0.2")), c("Condition", "eps", "Procedure", "theta")),
  "ICAR bank: exact conditional results under the model and under contamination")
if (!is.null(irs)) {
  rules <- c(raw_lo = "Raw score (k-)", raw_hi = "Raw score (k+)", wle_point = "WLE point", wald = "WLE Wald three-way", eap = "EAP three-way", exact = "Exact three-way")
  cuts_r <- sort(unique(as.numeric(sub("^exact_", "", grep("^exact_[0-9.]+$", names(irs), value = TRUE)))))
  tab <- do.call(rbind, lapply(cuts_r, function(ct) do.call(rbind, lapply(names(rules), function(r) { v <- irs[[paste0(r, "_", ct)]]
    data.frame(Cut = ct, Rule = rules[[r]], Above = mean(v == "above"), Below = mean(v == "below"), Undetermined = mean(v == "undetermined"), N = length(v)) }))))
  add_sheet("S18_ICAR_Respondents", tab, "ICAR respondents (complete cases): classification by rule",
            sprintf("N = %d; WLE mean %.2f, SD %.2f; screens flagged %.1f%% (size 0.05) and %.1f%% (size 0.20).", nrow(irs), mean(irs$wle), sd(irs$wle),
                    100 * (1 - mean(irs$pass_s0.05)), 100 * (1 - mean(irs$pass_s0.2))))
  xt <- do.call(rbind, lapply(cuts_r, function(ct) { x <- as.data.frame(table(Exact = irs[[paste0("exact_", ct)]], Wald = irs[[paste0("wald_", ct)]]), stringsAsFactors = FALSE)
    names(x)[3] <- "Count"; cbind(Cut = ct, x[x$Count > 0, ]) }))
  Jmax <- max(irs$raw); ps <- irs[irs$raw == Jmax, ]
  und <- irs[irs[[paste0("exact_", cuts_r[1])]] == "undetermined", ]
  agree <- do.call(rbind, lapply(cuts_r, function(ct) data.frame(Cut = ct,
    Agree_raw_kplus_WLE = mean(irs[[paste0("raw_hi_", ct)]] == irs[[paste0("wle_point_", ct)]]),
    Agree_raw_kminus_WLE = mean(irs[[paste0("raw_lo_", ct)]] == irs[[paste0("wle_point_", ct)]]))))
  detail <- data.frame(Quantity = c("Perfect scorers (n)", "Perfect score: WLE", "Perfect score: CSEM", "Perfect score: Wald lower limit",
                                    "Perfect score: exact lower limit", "Perfect score: information at WLE", "Perfect score: EAP", "Perfect score: EAP lower limit",
                                    sprintf("Undetermined by exact rule at %.2f: WLE min", cuts_r[1]), "  WLE max", "  raw score min", "  raw score max"),
                       Value = c(nrow(ps), ps$wle[1], ps$se_wle[1], ps$wle[1] - qnorm(0.975) * ps$se_wle[1], ps$exact_L[1], 1 / ps$se_wle[1]^2,
                                 ps$eap[1], ps$eap[1] - qnorm(0.975) * ps$psd[1], min(und$wle), max(und$wle), min(und$raw), max(und$raw)))
  add_sheet("S19_ICAR_Resp_Detail", detail, "ICAR respondents: perfect scorers and respondents left undetermined", "Information at the WLE = 1 / CSEM^2.")
  add_sheet("S20_ICAR_CrossTab", xt, "ICAR respondents: exact versus Wald three-way classification (counts)")
  add_sheet("S21_ICAR_Agreement", agree, "ICAR respondents: agreement of raw-score rules with the WLE point rule")
}
if (!is.null(iboot)) {
  q <- t(apply(iboot[, setdiff(names(iboot), "draw"), drop = FALSE], 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  add_sheet("S22_ICAR_Bootstrap", data.frame(Quantity = rownames(q), `2.5%` = q[, 1], Median = q[, 2], `97.5%` = q[, 3], Span = q[, 3] - q[, 1], check.names = FALSE),
            sprintf("ICAR calibration bootstrap (B = %d refits)", nrow(iboot)), "Span = width of the 95% bootstrap interval.")
}

## ---- A: article values ------------------------------------------------------------------------
AV <- list()
f3 <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 3))
rg <- function(x) if (all(is.na(x))) "NA" else paste(f3(min(x, na.rm = TRUE)), "to", f3(max(x, na.rm = TRUE)))
add <- function(section, statement, article, computed, sheet, how)
  AV[[length(AV) + 1]] <<- data.frame(Section = section, Statement = statement, Article = article, Recomputed = computed, Sheet = sheet, Filter = how, stringsAsFactors = FALSE)
safe <- function(expr) tryCatch(expr, error = function(e) "NA")

if (!is.null(ms) && !is.null(g)) {
  m0 <- ms[ms$eps == 0, ]; g0 <- g[g$eps == 0, ]
  cmin <- function(est, lay, Js) sapply(Js, function(J) min(sel(g0, estimator = est, layout = lay, J = J)$coverage))
  sec <- "Interval validity"
  add(sec, "Exact interval: minimum conditional coverage", "0.950", safe(f3(min(sel(g0, estimator = "EXACT_CI")$coverage))), "S3", "Procedure = Exact; min of coverage")
  add(sec, "Exact interval: marginal coverage", "0.960 to 0.975", safe(rg(sel(m0, estimator = "EXACT_CI")$coverage)), "S2", "Procedure = Exact")
  add(sec, "WLE Wald: minimum conditional coverage, peaked", "0.86 to 0.90", safe(rg(cmin("WLE", "peaked", c(8, 12, 16, 20)))), "S3", "WLE; peaked; min over theta per J")
  add(sec, "WLE Wald: minimum conditional coverage, wide", "0.93 to 0.94", safe(rg(cmin("WLE", "wide", c(8, 12, 16, 20)))), "S3", "WLE; wide; min over theta per J")
  add(sec, "EAP: max |marginal coverage - 0.95|", "0.003", safe(f3(max(abs(sel(m0, estimator = "EAP")$coverage - 0.95)))), "S2", "Procedure = EAP")
  add(sec, "EAP: minimum conditional coverage, wide", "0.43 to 0.62", safe(rg(cmin("EAP", "wide", c(8, 12, 16, 20)))), "S3", "EAP; wide; min over theta per J")
  add(sec, "EAP: minimum conditional coverage, peaked, J >= 16", "0.48 to 0.55", safe(rg(cmin("EAP", "peaked", c(16, 20)))), "S3", "EAP; peaked; J = 16, 20")
  add(sec, "EAP: coverage at theta = +/-3, peaked, J <= 12", "0", safe(rg(sel(g0, estimator = "EAP", layout = "peaked", J = c(8, 12), theta = c(-3, 3))$coverage)), "S3", "EAP; peaked; J = 8, 12; theta = -3, 3")
  add(sec, "Huber: marginal coverage", "0.959 to 0.981", safe(rg(sel(m0, estimator = "HU")$coverage)), "S2", "Procedure = Huber")
  r2 <- add_ratios_to(m0[m0$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI"), ])
  add(sec, "Huber RMSE / WLE RMSE", "1.11 to 1.73 (11% to 73%)", safe(rg(sel(r2, estimator = "HU")$rmse_ratio_to_WLE)), "S2", "Huber; rmse_ratio_to_WLE")
  add(sec, "EAP RMSE / WLE RMSE", "0.82 to 0.91", safe(rg(sel(r2, estimator = "EAP")$rmse_ratio_to_WLE)), "S2", "EAP; rmse_ratio_to_WLE")
  add(sec, "Exact width / WLE Wald width", "1.07 to 1.21 (7% to 21%)", safe(rg(sel(r2, estimator = "EXACT_CI")$width_ratio_to_WLE)), "S2", "Exact; width_ratio_to_WLE")

  sec <- "Two-outcome rules"
  dif <- unlist(lapply(c(0.25, 1.25, 2.25), function(ct) { c <- paste0("wrong_point_", ct)
    m <- merge(sel(m0, estimator = "RAW_hi")[, c("layout", "J", c)], sel(m0, estimator = "WLE")[, c("layout", "J", c)], by = c("layout", "J")); m[[2 + 1]] - m[[2 + 2]] }))
  add(sec, "Raw (k+) wrong rate - WLE wrong rate, all cells and cuts", "-0.008 to 0.034", safe(rg(dif)), "S4", "Raw (k+); Diff_from_WLE")
  add(sec, "Raw (k-) vs WLE, wide, J = 8, cut 1.25", "0.280 vs 0.090", safe(paste(f3(v1(m0, "wrong_point_1.25", estimator = "RAW_lo", layout = "wide", J = 8)), "vs", f3(v1(m0, "wrong_point_1.25", estimator = "WLE", layout = "wide", J = 8)))), "S4", "wide; J = 8; cut 1.25")
  for (th in c(1.0, 1.5)) add(sec, sprintf("Conditional wrong, wide, J = 20, cut 1.25, theta = %.1f (k-, k+, WLE)", th), if (th == 1) "0.59, 0.36, 0.28" else "0.09, 0.23, 0.28",
    safe(paste(sapply(c("RAW_lo", "RAW_hi", "WLE"), function(e) f3(v1(g0, "wrong_point_1.25", estimator = e, layout = "wide", J = 20, theta = th))), collapse = ", ")), "S5", sprintf("wide; J = 20; cut 1.25; theta = %.1f", th))
  rel <- function(ct) { c <- paste0("wrong_point_", ct); m <- merge(sel(m0, estimator = "RUDNER")[, c("layout", "J", c)], sel(m0, estimator = "WLE")[, c("layout", "J", c)], by = c("layout", "J")); m[[3]] / m[[4]] - 1 }
  add(sec, "Rudner relative error, cut 0.25", "up to +16%", safe(rg(rel(0.25))), "S4", "Rudner; cut 0.25; Relative_to_WLE")
  add(sec, "Rudner relative error, cut 2.25", "-22% to +48%", safe(rg(rel(2.25))), "S4", "Rudner; cut 2.25; Relative_to_WLE")
  add(sec, "Rudner vs exact WLE wrong, theta = 2.5, cut 1.25, peaked, J = 20", "0.086 vs 0.013",
      safe(paste(f3(v1(g0, "wrong_point_1.25", estimator = "RUDNER", layout = "peaked", J = 20, theta = 2.5)), "vs", f3(v1(g0, "wrong_point_1.25", estimator = "WLE", layout = "peaked", J = 20, theta = 2.5)))), "S5", "peaked; J = 20; cut 1.25; theta = 2.5")
  add(sec, "EAP point rule at tail cut, peaked: accuracy (all J)", "0.981", safe(rg(1 - sel(m0, estimator = "EAP", layout = "peaked")$wrong_point_2.25)), "S4", "EAP; peaked; cut 2.25")
  add(sec, "EAP point rule at tail cut, peaked: consistency (all J)", "1.000", safe(rg(sel(m0, estimator = "EAP", layout = "peaked")$cons_point_2.25)), "S4", "EAP; peaked; cut 2.25")
  add(sec, "WLE point rule at tail cut, peaked, J = 8: consistency", "1.000", safe(f3(v1(m0, "cons_point_2.25", estimator = "WLE", layout = "peaked", J = 8))), "S4", "WLE; peaked; J = 8; cut 2.25")
  add(sec, "EAP conditional wrong at theta = 2.5 and 3, peaked, cut 2.25", "1.0", safe(rg(sel(g0, estimator = "EAP", layout = "peaked", theta = c(2.5, 3))$wrong_point_2.25)), "S5", "EAP; peaked; cut 2.25; theta 2.5, 3")

  sec <- "Three-way rules"; I3 <- m0[m0$estimator %in% c("WLE", "EAP", "EXACT_CI"), ]
  add(sec, "Maximum marginal wrong-decision probability", "0.003", safe(f3(max(unlist(I3[, paste0("wrong3_", c(0.25, 1.25, 2.25))])))), "S6", "Wrong; all rules and cuts")
  for (ct in c(0.25, 1.25, 2.25)) add(sec, sprintf("Marginal undetermined probability, cut %.2f", ct), c("0.47 to 0.95", "0.28 to 0.57", "0.06 to 0.27")[match(ct, c(0.25, 1.25, 2.25))],
    safe(rg(I3[[paste0("undet_", ct)]])), "S6", sprintf("Undetermined; cut %.2f", ct))
  add(sec, "Exact rule decided, theta >= 2, cut 1.25, peaked, J = 20", "0.17 to 0.53", safe(rg(1 - sel(g0, estimator = "EXACT_CI", layout = "peaked", J = 20, theta = c(2, 2.5, 3))$undet_1.25)), "S7", "Exact; peaked; J = 20; cut 1.25; 1 - Undetermined")
  add(sec, "WLE Wald undetermined, theta >= 2, cut 1.25, peaked, J = 20", "1.000", safe(rg(sel(g0, estimator = "WLE", layout = "peaked", J = 20, theta = c(2, 2.5, 3))$undet_1.25)), "S7", "WLE; peaked; J = 20; cut 1.25")
  add(sec, "WLE Wald wrong, theta = 2.5, cut 2.25, peaked, J = 20", "0.039", safe(f3(v1(g0, "wrong3_2.25", estimator = "WLE", layout = "peaked", J = 20, theta = 2.5))), "S7", "WLE; peaked; J = 20; cut 2.25; theta = 2.5")
  add(sec, "EAP wrong, theta = 2.5, cut 2.25, maximum over cells", "0.12 (wide, J = 8)", safe(f3(max(sel(g0, estimator = "EAP", theta = 2.5)$wrong3_2.25))), "S7", "EAP; cut 2.25; theta = 2.5")
  add(sec, "Exact rule: maximum conditional wrong, any theta and cut", "0.011", safe(f3(max(unlist(sel(g0, estimator = "EXACT_CI")[, paste0("wrong3_", c(0.25, 1.25, 2.25))])))), "S7", "Exact; Wrong; all")
  add(sec, "EAP undetermined at tail cut (marginal)", "0.06 to 0.11", safe(rg(sel(m0, estimator = "EAP")$undet_2.25)), "S6", "EAP; cut 2.25")

  sec <- "Aberrant responding"; C <- ms[ms$eps > 0 & ms$estimator %in% c("WLE", "EAP", "HU", "EXACT_CI"), ]
  add(sec, "Marginal coverage, all-correct, eps = 0.10", "0.87 to 0.95", safe(rg(sel(C, cond = "up", eps = 0.1)$coverage)), "S8", "All-correct; eps 0.10")
  add(sec, "Marginal coverage, all-correct, eps = 0.20", "0.73 to 0.91", safe(rg(sel(C, cond = "up", eps = 0.2)$coverage)), "S8", "All-correct; eps 0.20")
  mc0 <- merge(C, ms[ms$eps == 0, c("layout", "J", "estimator", "coverage")], by = c("layout", "J", "estimator"), suffixes = c("", "_model"))
  add(sec, "Maximum coverage loss under guessing", "0.02", safe(f3(max(sel(mc0, cond = "guess")$coverage_model - sel(mc0, cond = "guess")$coverage))), "S8", "Guessing; coverage_loss")
  for (e in c("WLE", "EXACT_CI", "HU", "EAP")) add(sec, paste("Conditional minimum coverage, eps = 0.20, all-correct/all-incorrect:", e), c(WLE = "0.04", EXACT_CI = "0.12", HU = "0.12", EAP = "0.00")[[e]],
    safe(f3(min(sel(g, estimator = e, eps = 0.2, cond = c("up", "down"))$coverage))), "S9", paste("eps 0.20; all-correct, all-incorrect;", e))
  add(sec, "WLE point wrong, cut 0.25, wide, J = 20: model -> all-correct 0.20", "0.106 -> 0.274",
      safe(paste(f3(v1(ms, "wrong_point_0.25", estimator = "WLE", layout = "wide", J = 20, eps = 0)), "->", f3(v1(ms, "wrong_point_0.25", estimator = "WLE", layout = "wide", J = 20, eps = 0.2, cond = "up")))), "S4, S8", "wide; J = 20; cut 0.25")
  add(sec, "Raw (k+) wrong, cut 0.25, wide, J = 20: model -> all-correct 0.20", "0.134 -> 0.318",
      safe(paste(f3(v1(ms, "wrong_point_0.25", estimator = "RAW_hi", layout = "wide", J = 20, eps = 0)), "->", f3(v1(ms, "wrong_point_0.25", estimator = "RAW_hi", layout = "wide", J = 20, eps = 0.2, cond = "up")))), "S4, S8", "wide; J = 20; cut 0.25")
  for (e in c("EXACT_CI", "WLE", "EAP")) add(sec, paste("Maximum three-way wrong under contamination:", e), c(EXACT_CI = "0.035", WLE = "0.035", EAP = "0.040")[[e]],
    safe(f3(max(unlist(sel(ms, estimator = e)[ms$eps[ms$estimator == e] > 0, paste0("wrong3_", c(0.25, 1.25, 2.25))])))), "S8", paste(e, "; all cuts"))
  hw <- merge(sel(C, estimator = "HU"), sel(C, estimator = "WLE"), by = c("layout", "J", "eps", "cond"), suffixes = c("_hu", "_wle"))
  add(sec, "Huber - WLE coverage under contamination (max)", "0.05", safe(f3(max(hw$coverage_hu - hw$coverage_wle))), "S8", "Huber vs WLE")
  add(sec, "Huber RMSE / WLE RMSE under contamination", "1.02 to 1.67", safe(rg(hw$rmse_hu / hw$rmse_wle)), "S8", "Huber; rmse_ratio_to_WLE")
  for (ct in c(0.25, 2.25)) add(sec, sprintf("|Huber - WLE| point wrong, cut %.2f (max)", ct), if (ct == 0.25) "0.02" else "0.06",
    safe(f3(max(abs(hw[[paste0("wrong_point_", ct, "_hu")]] - hw[[paste0("wrong_point_", ct, "_wle")]])))), "S8", sprintf("cut %.2f", ct))

  sec <- "Person-fit screening"; E <- ms[ms$estimator == "EXACT_CI", ]; W <- ms[ms$estimator == "WLE", ]
  for (lay in c("peaked", "wide")) for (sz in c("0.05", "0.2")) {
    p_c <- v1(E, paste0("pass_s", sz), layout = lay, J = 20, eps = 0.2, cond = "up"); p_m <- v1(E, paste0("pass_s", sz), layout = lay, J = 20, eps = 0)
    add(sec, sprintf("Pass rate (contamination; model) and implied detection, %s, size %s, all-correct 0.20, J = 20", lay, sz),
        c(peaked0.05 = "0.894; 0.951; 0.06", wide0.05 = "0.751; 0.950; 0.21", peaked0.2 = "0.695; 0.801; 0.13", wide0.2 = "0.522; 0.801; 0.35")[[paste0(lay, sz)]],
        safe(paste0(f3(p_c), "; ", f3(p_m), "; ", f3(1 - p_c / p_m))), "S10", sprintf("%s; J = 20; all-correct 0.20; size %s", lay, sz))
    for (nm in c("E", "W")) { D <- get(nm)
      add(sec, sprintf("%s coverage all -> given pass, %s, size %s, all-correct 0.20, J = 20", if (nm == "E") "Exact" else "WLE Wald", lay, sz),
          c(Epeaked0.05 = "0.796 -> 0.812", Epeaked0.2 = "0.796 -> 0.829", Ewide0.05 = "0.749 -> 0.793", Ewide0.2 = "0.749 -> 0.818",
            Wpeaked0.05 = "0.786 -> 0.806", Wpeaked0.2 = "0.786 -> 0.828", Wwide0.05 = "0.731 -> 0.784", Wwide0.2 = "0.731 -> 0.815")[[paste0(nm, lay, sz)]],
          safe(paste(f3(v1(D, "coverage", layout = lay, J = 20, eps = 0.2, cond = "up")), "->",
                     f3(v1(D, paste0("cov_pass_s", sz), layout = lay, J = 20, eps = 0.2, cond = "up") / v1(D, paste0("pass_s", sz), layout = lay, J = 20, eps = 0.2, cond = "up")))),
          "S10", sprintf("%s; J = 20; all-correct 0.20; cov_given_pass_s%s", lay, sz)) }
  }
  Ec <- E[E$eps > 0, ]
  for (sz in c("0.05", "0.2")) add(sec, sprintf("Largest coverage gain from screening, size %s", sz), if (sz == "0.05") "0.044" else "0.069",
    safe(f3(max(Ec[[paste0("cov_pass_s", sz)]] / Ec[[paste0("pass_s", sz)]] - Ec$coverage))), "S10", paste0("coverage_gain_s", sz, "; max"))
  for (sz in c("0.05", "0.2")) { gm <- merge(sel(E, cond = "guess"), E[E$eps == 0, c("layout", "J", paste0("pass_s", sz))], by = c("layout", "J"), suffixes = c("", "_m"))
    add(sec, sprintf("Pass-rate reduction under guessing, size %s", sz), "0.01 to 0.10", safe(rg(gm[[paste0("pass_s", sz, "_m")]] - gm[[paste0("pass_s", sz)]])), "S10", "Guessing; pass_model - pass") }
  add(sec, "Exact wrong at cut 0.25, wide, J = 20, all-correct 0.20: all -> given pass (size 0.20)", "0.035 -> 0.031",
      safe(paste(f3(v1(E, "wrong3_0.25", layout = "wide", J = 20, eps = 0.2, cond = "up")), "->", f3(v1(E, "wrong3_pass_s0.2_0.25", layout = "wide", J = 20, eps = 0.2, cond = "up") / v1(E, "pass_s0.2", layout = "wide", J = 20, eps = 0.2, cond = "up")))), "S10", "wide; J = 20; all-correct 0.20")
  add(sec, "theta = -3, wide, J = 20, all-correct 0.20: pass (size 0.20); coverage given pass", "0.12; 0.41",
      safe(paste0(f3(v1(g, "pass_s0.2", estimator = "EXACT_CI", layout = "wide", J = 20, eps = 0.2, cond = "up", theta = -3)), "; ",
                  f3(v1(g, "cov_pass_s0.2", estimator = "EXACT_CI", layout = "wide", J = 20, eps = 0.2, cond = "up", theta = -3) / v1(g, "pass_s0.2", estimator = "EXACT_CI", layout = "wide", J = 20, eps = 0.2, cond = "up", theta = -3)))), "S9", "Exact; wide; J = 20; theta = -3")
}

if (!is.null(ims) && !is.null(ig)) {
  sec <- "Empirical study (ICAR)"; i0 <- ims[ims$eps == 0, ]; ig0 <- ig[ig$eps == 0, ]
  if (!is.null(ibank)) add(sec, "TCC(0.67); TCC(1.28)", "11.00; 13.06", safe(paste(formatC(tcc(0.67), format = "f", digits = 2), formatC(tcc(1.28), format = "f", digits = 2), sep = "; ")), "S12", "note of S12")
  add(sec, "Marginal wrong, cut 0.67: raw k-, WLE, raw k+", "0.095, 0.108, 0.128", safe(paste(sapply(c("RAW_lo", "WLE", "RAW_hi"), function(e) f3(v1(i0, "wrong_point_0.67", estimator = e))), collapse = ", ")), "S15", "two-outcome; cut 0.67")
  add(sec, "Marginal wrong, cut 1.28: raw k+, WLE, raw k-", "0.054, 0.054, 0.075", safe(paste(sapply(c("RAW_hi", "WLE", "RAW_lo"), function(e) f3(v1(i0, "wrong_point_1.28", estimator = e))), collapse = ", ")), "S15", "two-outcome; cut 1.28")
  add(sec, "Rudner vs exact WLE, marginal: cut 0.67; cut 1.28", "0.108 vs 0.108; 0.056 vs 0.054",
      safe(paste0(f3(v1(i0, "wrong_point_0.67", estimator = "RUDNER")), " vs ", f3(v1(i0, "wrong_point_0.67", estimator = "WLE")), "; ", f3(v1(i0, "wrong_point_1.28", estimator = "RUDNER")), " vs ", f3(v1(i0, "wrong_point_1.28", estimator = "WLE")))), "S15", "Rudner and WLE")
  add(sec, "Rudner vs exact WLE at theta = 2.5 and 3, cut 1.28", "0.086, 0.085 vs 0.016, 0.002",
      safe(paste0(f3(v1(ig0, "wrong_point_1.28", estimator = "RUDNER", theta = 2.5)), ", ", f3(v1(ig0, "wrong_point_1.28", estimator = "RUDNER", theta = 3)), " vs ",
                  f3(v1(ig0, "wrong_point_1.28", estimator = "WLE", theta = 2.5)), ", ", f3(v1(ig0, "wrong_point_1.28", estimator = "WLE", theta = 3)))), "S17", "cut 1.28; theta 2.5, 3")
  add(sec, "EAP vs WLE point wrong, theta = 1.5, cut 1.28", "0.51 vs 0.35", safe(paste(f3(v1(ig0, "wrong_point_1.28", estimator = "EAP", theta = 1.5)), "vs", f3(v1(ig0, "wrong_point_1.28", estimator = "WLE", theta = 1.5)))), "S17", "theta = 1.5; cut 1.28")
  add(sec, "Minimum conditional coverage: exact; WLE Wald", "0.950; 0.91", safe(paste(f3(min(sel(ig0, estimator = "EXACT_CI")$coverage)), f3(min(sel(ig0, estimator = "WLE")$coverage)), sep = "; ")), "S17", "Model; min over theta")
  add(sec, "EAP coverage at theta = -3; 3", "0.62; 0.66", safe(paste(f3(v1(ig0, "coverage", estimator = "EAP", theta = -3)), f3(v1(ig0, "coverage", estimator = "EAP", theta = 3)), sep = "; ")), "S17", "EAP; theta -3, 3")
  add(sec, "WLE SE ratio at theta = -3; 3", "2.0 to 2.3", safe(paste(f3(v1(ig0, "se_ratio", estimator = "WLE", theta = -3)), f3(v1(ig0, "se_ratio", estimator = "WLE", theta = 3)), sep = "; ")), "S17", "WLE; se_ratio")
  add(sec, "Undetermined at cut 0.67, theta = 1.5, 2.0, 2.5: exact | WLE Wald", "0.58, 0.24, 0.07 | 0.78, 0.47, 0.22",
      safe(paste(paste(sapply(c(1.5, 2, 2.5), function(t) f3(v1(ig0, "undet_0.67", estimator = "EXACT_CI", theta = t))), collapse = ", "), "|",
                 paste(sapply(c(1.5, 2, 2.5), function(t) f3(v1(ig0, "undet_0.67", estimator = "WLE", theta = t))), collapse = ", "))), "S17", "cut 0.67")
  add(sec, "Undetermined at cut 1.28 for theta >= 1.5 (min over the three rules)", "> 0.98", safe(f3(min(sel(ig0, estimator = c("WLE", "EAP", "EXACT_CI"), theta = c(1.5, 2, 2.5, 3))$undet_1.28))), "S17", "cut 1.28; theta >= 1.5")
  cu <- sel(ims, eps = 0.2, cond = "up")
  add(sec, "Marginal coverage, all-correct 0.20: WLE, EAP, Huber, exact", "0.820, 0.793, 0.842, 0.799", safe(paste(sapply(c("WLE", "EAP", "HU", "EXACT_CI"), function(e) f3(v1(cu, "coverage", estimator = e))), collapse = ", ")), "S16", "All-correct; eps 0.20")
  add(sec, "Marginal coverage under guessing (range)", "0.93 to 0.96", safe(rg(sel(ims, cond = "guess", estimator = c("WLE", "EAP", "HU", "EXACT_CI"))$coverage)), "S16", "Guessing")
  add(sec, "Coverage at theta = -3, all-correct 0.20: WLE, EAP, Huber, exact", "0.09, 0.02, 0.23, 0.17",
      safe(paste(sapply(c("WLE", "EAP", "HU", "EXACT_CI"), function(e) f3(v1(ig, "coverage", estimator = e, eps = 0.2, cond = "up", theta = -3))), collapse = ", ")), "S17", "All-correct 0.20; theta = -3")
  add(sec, "Width, all-correct 0.20: Huber vs WLE Wald", "2.07 vs 1.79", safe(paste(f3(v1(cu, "width", estimator = "HU")), "vs", f3(v1(cu, "width", estimator = "WLE")))), "S16", "All-correct 0.20; width")
  add(sec, "Exact coverage all -> given pass, all-correct 0.20 (sizes 0.05, 0.20)", "0.799 -> 0.827, 0.847",
      safe(paste(f3(v1(cu, "coverage", estimator = "EXACT_CI")), "->", f3(v1(cu, "cov_pass_s0.05", estimator = "EXACT_CI") / v1(cu, "pass_s0.05", estimator = "EXACT_CI")), ",",
                 f3(v1(cu, "cov_pass_s0.2", estimator = "EXACT_CI") / v1(cu, "pass_s0.2", estimator = "EXACT_CI")))), "S17 (marginal from ALL_marginal_summary)", "Exact; all-correct 0.20")
}
if (!is.null(irs)) {
  sec <- "Empirical study (ICAR respondents)"
  add(sec, "N; WLE mean; WLE SD", "1,248; 0.04; 1.07", safe(paste(nrow(irs), f3(mean(irs$wle)), f3(sd(irs$wle)), sep = "; ")), "S18", "note of S18")
  add(sec, "Exact vs Wald 'above' at cut 0.67", "0.092 vs 0.057", safe(paste(f3(mean(irs$`exact_0.67` == "above")), "vs", f3(mean(irs$`wald_0.67` == "above")))), "S18", "cut 0.67")
  add(sec, "Above by exact, undetermined by Wald (count); below by Wald, undetermined by exact (count)", "44; 2",
      safe(paste(sum(irs$`exact_0.67` == "above" & irs$`wald_0.67` == "undetermined"), sum(irs$`exact_0.67` == "undetermined" & irs$`wald_0.67` == "below"), sep = "; ")), "S20", "cut 0.67")
  add(sec, "Undetermined by exact at 0.67: WLE range; raw-score range", "-0.13 to 1.51; 7 to 14", safe(paste0(rg(und$wle), "; ", min(und$raw), " to ", max(und$raw))), "S19", "detail rows")
  add(sec, "Perfect scorers: n; exact lower limit; Wald lower limit; information", "30; 1.251; 0.68; 0.9",
      safe(paste(nrow(ps), f3(ps$exact_L[1]), f3(ps$wle[1] - qnorm(0.975) * ps$se_wle[1]), f3(1 / ps$se_wle[1]^2), sep = "; ")), "S19", "detail rows")
  add(sec, "Agreement with WLE: k+ and k- at 0.67; at 1.28", "0.978, 0.933; 0.983, 0.959",
      safe(paste(paste(f3(agree$Agree_raw_kplus_WLE[1]), f3(agree$Agree_raw_kminus_WLE[1]), sep = ", "), paste(f3(agree$Agree_raw_kplus_WLE[2]), f3(agree$Agree_raw_kminus_WLE[2]), sep = ", "), sep = "; ")), "S21", "")
  add(sec, "Flagged by screens (sizes 0.05, 0.20)", "0.095, 0.225", safe(paste(f3(1 - mean(irs$pass_s0.05)), f3(1 - mean(irs$pass_s0.2)), sep = ", ")), "S18", "note of S18")
}
if (!is.null(iboot)) {
  sp <- apply(iboot[, setdiff(names(iboot), "draw"), drop = FALSE], 2, function(x) diff(quantile(x, c(0.025, 0.975), na.rm = TRUE)))
  add("Empirical study (ICAR bootstrap)", "Widest bootstrap interval (raw-WLE disagreement at 0.67); range of the others", "0.053; 0.008 to 0.043",
      safe(paste0(f3(sp[["resp_rawhi_vs_wle_0.67"]]), "; ", rg(sp[names(sp) != "resp_rawhi_vs_wle_0.67"]))), "S22", "Span column")
}
if (length(AV)) {
  A <- do.call(rbind, AV)
  add_sheet("A_Article_Values", A, "Values quoted in the article and their recomputation from the output files",
            "Article = value as written in the manuscript; Recomputed = value computed here; Sheet and Filter show where to find it. Small differences are rounding.")
}

## ---- README and order ----------------------------------------------------------------------------
addWorksheet(wb, "README")
writeData(wb, "README", "Supplementary materials: Enumeration-Based Evaluation of Cut Score Decisions Using Reported Item Response Theory Precision in Short Tests", startRow = 1)
addStyle(wb, "README", st_title, rows = 1, cols = 1)
writeData(wb, "README", c("All quantities are computed exactly by enumeration of response patterns (no simulation error).",
                          "Exact study: wide (b on [-2, 2]) and peaked (b on [-0.75, 0.75]) banks; a cycles 0.8/1.2/1.6/2.0; J = 8, 12, 16, 20; cuts 0.25, 1.25, 2.25.",
                          "Empirical study: 16 ICAR items (psychTools::ability), 1,248 complete-case respondents; cuts 0.67, 1.28.",
                          "Start with sheet A_Article_Values, which maps each number in the article to the sheet and filter that reproduce it."), startRow = 2)
writeData(wb, "README", index, startRow = 7, headerStyle = st_head)
setColWidths(wb, "README", cols = 1:2, widths = c(26, 110))
nm <- names(wb); first <- c("README", "A_Article_Values"); worksheetOrder(wb) <- match(c(intersect(first, nm), setdiff(nm, first)), nm)

saveWorkbook(wb, out_file, overwrite = TRUE)
cat("\nWritten:", out_file, "\n")
if (length(AV)) { A <- do.call(rbind, AV); cat("Article values recomputed:", nrow(A), "(see sheet A_Article_Values)\n") }
