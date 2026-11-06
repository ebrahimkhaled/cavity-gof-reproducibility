# E51 -- confirmation of the tied recalibration (tau = kappa b, rho^2 from adjusted R^2). See README.
#
# PRE-REGISTRATION (written 2026-11-06 after E50, before any E51 result).
# E50 result that motivates this: at n = 400 the tied fit and tau = 0 gave almost the same decisions, so n = 400 cannot show what the tie is for. The tie is
#   meant to stop the drift of the two-parameter test at large n (size .100 at n = 20000 in the bent design, E25). E50's I2 + R2 pairing was picked after seeing
#   E50, so it needs fresh seeds. Nothing in the method changes: same estimates, same clamp, same G = 10 and chi2_{G-2}.
# Part 1 (drift): p = 50, AR .9, lambda = 2n, n = 2000, 5000, 20000; logistic ridge index and I1.
# Part 2 (growth): p/n = .25, AR .7, lambda/n = 1.04, n = 800, 1600, 3200; I1 and I2.
# Part 3 (fresh seeds): E47c cells, seed + 7e6.
# PASS needs ALL of:
#   (i) Part 1, n = 20000: T0 size above .0695 (the design still discriminates) and R2 size at most .0695, for the logistic index and for I1;
#   (ii) Part 2: R2 size in [.0305, .0695] at every n, for I1 and I2;
#   (iii) Part 3, I2 + R2: size in band in A_p5, B_p40, B_p100, B_p160, H5_k25_ar03, at most .0695 at binary_k25 and t3_k25; quadratic power not below d3 by more
#         than 2 paired SE; sd of the implied tau at B_p100 below one third of its median.
# If (i) fails: the tie has no use over tau = 0 and the paper keeps d = 3. No tuning.
suppressPackageStartupMessages({ library(parallel) })
args <- commandArgs(trailingOnly = TRUE); W <- if (length(args) >= 1) as.integer(args[1]) else 14; PHASE <- if (length(args) >= 2) args[2] else "run"
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUTF <- file.path(ROOT, "results", "E51_confirm_perrep.csv")
src <- readLines(file.path(ROOT, "sims", "E47c_standardised.R"))
eval(parse(text = src[grep("^SEED0", src):(grep("^# default statistic", src) - 1)]))
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
eval(parse(text = src[grep("^loo_linear", src):grep("^std <-", src)]))
s50 <- readLines(file.path(ROOT, "sims", "E50_tau_separate.R")); eval(parse(text = s50[grep("^cells_e", s50):(grep("^one <- ", s50) - 1)]))

D <- list()
for (n in c(2000, 5000, 20000)) { d <- mk(paste0("bent_n", n), n, 50, 2 * n, ar = .9, seed = 20261205 + n * 1000); d$part <- 1; d$useFull <- TRUE; d$xf <- TRUE; D[[length(D) + 1]] <- d }
for (n in c(800, 1600, 3200)) { d <- mk(paste0("k25_n", n), n, n / 4, 1.04 * n, seed = 20261305 + n * 1000); d$part <- 2; d$useFull <- FALSE; d$xf <- FALSE; D[[length(D) + 1]] <- d }
for (nm in c("A_p5", "B_p40", "B_p100", "B_p160", "H5_k25_ar03", "binary_k25", "t3_k25", "B_p40_quad1", "B_p100_quad1", "A_cubic1")) {
  d <- Filter(function(c) c$name == nm, CELLS)[[1]]; d$seed <- d$seed + 7000000; d$name <- paste0("fresh_", nm); d$part <- 3; d$useFull <- TRUE; d$xf <- FALSE; D[[length(D) + 1]] <- d }

one <- function(r, d) {
  z <- gen(r, d); X <- z$X; y <- z$y; n <- d$n; sx <- d$seed + r + 5000000; a <- adjR2(X, y); pen <- 4 * d$lamI
  L1 <- loo_linear(X, y, pen)
  out <- score("I1", L1, y, if (d$xf) rho_XF(X, y, sx, function(Xa, ya, s) lin_b(Xa, ya, pen * s)) else NA_real_, rho_R2(L1, y, a), rho_OR(lin_b(X, y, pen), d))
  if (d$p < n - 2 && d$part != 1) { L2 <- loo_linear(X, y, 0); out <- c(out, score("I2", L2, y, NA_real_, rho_R2(L2, y, a), rho_OR(lin_b(X, y, 0), d))) }
  if (d$useFull) { b <- cf_ridge(cbind(1, X), y, d$lambda, diag(c(0, rep(1, d$p)))); xi <- cav_gof_fast(X, y, d$lambda, beta = b)$xi
    fl <- function(Xa, ya, s) cf_ridge(cbind(1, Xa), ya, d$lambda * s, diag(c(0, rep(1, ncol(Xa)))))[-1]
    out <- c(out, score("full", xi, y, if (d$xf) rho_XF(X, y, sx, fl) else NA_real_, rho_R2(xi, y, a), rho_OR(b[-1], d))) }
  cols <- as.vector(outer(c("I1", "I2", "full"), c("d3", "T0", "OR", "XF", "R2", "tau_d3", "tau_OR", "tau_XF", "tau_R2", "rho_OR", "rho_XF", "rho_R2"), paste, sep = "_"))
  res <- setNames(rep(NA_real_, length(cols)), cols); res[names(out)] <- out
  data.frame(rep = r, cell = d$name, part = d$part, t(res))
}
if (PHASE == "smoke") { for (k in c(1, 3, 6, 9)) print(system.time(o <- one(1, D[[k]]))); print(round(unlist(o[, c("I1_d3", "I1_R2", "I1_tau_R2", "I1_tau_OR")]), 3)); quit("no") }
done <- if (file.exists(OUTF)) { x <- read.csv(OUTF); paste(x$cell, x$rep) } else character()
cl <- makeCluster(W); clusterExport(cl, c("std", "gen", "one", "score", "d3", "tied", "kap", "cells_e", "adjR2", "rho_R2", "rho_XF", "rho_OR", "lin_b", "loo_linear", "a_for", "cf_expit", "statmod_gh", "cf_ridge", "cav_gof_fast", "cav_recal_fit", "cav_pvalue"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
for (d in D) for (st in seq(1, 1000, by = 200)) {
  rs <- st:(st + 199); rs <- rs[!paste(d$name, rs) %in% done]; if (!length(rs)) next
  out <- do.call(rbind, parLapply(cl, rs, one, d = d))     # every batch written to disk: resume by re-running
  write.table(out, OUTF, sep = ",", row.names = FALSE, col.names = !file.exists(OUTF), append = file.exists(OUTF))
  r5 <- function(v) round(mean(v < .05, na.rm = TRUE), 3)
  cat(d$name, max(rs), "I1 d3", r5(out$I1_d3), "T0", r5(out$I1_T0), "R2", r5(out$I1_R2), "full T0", r5(out$full_T0), "R2", r5(out$full_R2), "\n")
}
stopCluster(cl)
