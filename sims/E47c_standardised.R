# E47c -- CORRECTION RUN, same rule as E47. E47b showed the raw linear index (sd about 0.07) distorts the d = 3 optimiser (cubic .27 raw, .094 standardised),
# so every index is standardised to mean 0, sd 1 before the recalibration fit. Only that changes, plus the comparator F1 (logistic ridge, lambda = 1, p < n).
# E47 -- the cavity test on an EXPLICIT (linear) leave-one-out index. Trial 3 of the study (see README).
#
# PRE-REGISTRATION (written 2026-11-06 before any E47 result).
# Why: Prop 2 and the exact covariance identity hold for ANY estimator of the other observations, so the fit that
#   builds the index need not be the logistic fit. Under a Gaussian design least squares estimates the index direction
#   (Brillinger). A linear index makes the leave-one-out step exact and every influence quantity an explicit formula.
# Candidates (statistic, cells, recalibration and reference exactly the default: G = 10, d = 3, chi2_{G-2}):
#   I1 one-step ridge: b = (X1'X1 + 4 lambda D)^{-1} X1'y, the first Newton step from the null at v = 1/4, same lambda as
#      the fit (lasso cell: lambda = n). Index = exact leave-one-out (yhat_i - h_ii y_i) / (1 - h_ii).
#   I2 leave-one-out OLS (p < n only): no penalty, so the index is unchanged by any invertible map of the predictors.
# PASS for I1 needs ALL of:
#   (a) size in [0.0305, 0.0695] at 5% in A_p5, B_p40, B_p100, B_p160, H5_k25_ar03, lasso_sparse5_p2n (R = 1000);
#   (b) quadratic power at B_p40_quad1 and B_p100_quad1 not below the default full test by more than 2 paired SE
#       (reps 1-500, the reps E43 stored);
#   (c) outside Gaussian designs: size at most 0.0695 at binary_k25 and t3_k25. Sparse strong binary cells are reported
#       and the paper's usage rule applies if they fail; a failure there narrows the scope, it does not fail the trial.
# If (a) or (b) fails: stop, no tuning.
suppressPackageStartupMessages({ library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) as.integer(args[1]) else 16; PHASE <- if (length(args) >= 2) args[2] else "run"
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
OUTF <- file.path(ROOT, "results", "E47c_std_perrep.csv")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
SEED0 <- 20261111; SEED1 <- 20261147

a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda = NULL, ar = 0.7, norm = 2.31, prev = .5, indep = FALSE, sparse = 0, dep = "null", rel = 0,
               law = "gaussian", seed = SEED0, full = FALSE, lamI = NULL, b0 = NULL, R = 1000) {
  if (law %in% c("binary_sparse", "binary_rare")) { Sig <- diag(p) } else if (indep) { Sig <- diag(p); b0 <- c(.5, -.4, .3, -.3, .2) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-"))
    b0 <- if (sparse > 0) c(rep_len(c(.35, -.3, .25, -.2, .15), sparse), rep(0, p - sparse)) else rep_len(c(.35, -.3, .25, -.2, .15), p)
    b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  if (length(b0) < p) b0 <- c(b0, rep(0, p - length(b0)))
  list(name = name, n = n, p = p, lambda = lambda, lamI = if (is.null(lamI)) lambda else lamI, b0 = b0, Sc = chol(Sig), prev = prev,
       dep = dep, rel = rel, law = law, seed = seed, full = full, R = R, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0))) }
CELLS <- list(
  mk("A_p5", 500, 5, 100, indep = TRUE), mk("B_p40", 400, 40, 416), mk("B_p100", 400, 100, 416), mk("B_p160", 400, 160, 416),
  mk("H5_k25_ar03", 400, 100, 416, ar = .3), mk("prev08_k25", 400, 100, 416, prev = .08),
  mk("lasso_sparse5_p2n", 300, 600, sparse = 5, lamI = 300),
  mk("A_cubic1", 500, 5, 100, indep = TRUE, dep = "cubic", rel = 1), mk("B_p40_quad1", 400, 40, 416, dep = "quadratic", rel = 1),
  mk("B_p100_quad1", 400, 100, 416, dep = "quadratic", rel = 1),
  mk("k75", 400, 300, 416, seed = SEED1, full = TRUE), mk("ridge_p2n", 300, 600, 312, seed = SEED1, full = TRUE),
  mk("binary_k25", 400, 100, 416, law = "binary", seed = SEED1, full = TRUE), mk("t3_k25", 400, 100, 416, law = "t3", seed = SEED1, full = TRUE),
  mk("bin_sparse_p100", 400, 100, 100, law = "binary_sparse", b0 = c(1.2, -1.2, .9, -.9, .6), prev = .3, seed = SEED1, full = TRUE),
  mk("bin_rare_p100", 400, 100, 100, law = "binary_rare", b0 = c(2.5, .6, -.6, .6, -.6), prev = .2, seed = SEED1, full = TRUE))

gen <- function(r, d) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(d$seed + r)
  n <- d$n; p <- d$p
  X <- switch(d$law, gaussian = matrix(rnorm(n * p), n) %*% d$Sc, t3 = (matrix(rt(n * p, 3) / sqrt(3), n)) %*% d$Sc,
              binary = (matrix((rbinom(n * p, 1, .3) - .3) / sqrt(.21), n)) %*% d$Sc,
              binary_sparse = matrix(rbinom(n * p, 1, .3), n), binary_rare = cbind(rbinom(n, 1, .1), matrix(rbinom(n * (p - 1), 1, .3), n)))
  eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(d$dep, null = 0, cubic = d$sd_eta * d$rel * (Zs^3 - 3 * Zs) / sqrt(6), quadratic = d$sd_eta * d$rel * (Zs^2 - 1) / sqrt(2))
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta, d$prev); list(X = X, y = rbinom(n, 1, cf_expit(a + eta)))
}
# default statistic, E43's q3_S copied line for line
q3_S <- function(xi, y, G = 10, nodes = 24) {
  n <- length(y); gh <- statmod_gh(nodes)
  lam <- function(th) { P <- 0; for (k in seq_len(nodes)) P <- P + gh$w[k] * cf_expit(th[1] + th[2] * xi + exp(th[3]) * gh$x[k]); pmin(pmax(P, 1e-9), 1 - 1e-9) }
  nll <- function(th) { m <- lam(th); -sum(y * log(m) + (1 - y) * log(1 - m)) }
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients)
  op <- optim(c(rc, log(0.5)), nll, method = "L-BFGS-B", lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)))
  m <- lam(op$par); v <- m * (1 - m); g <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G); idx <- split(seq_len(n), g)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0) / sqrt(vapply(idx, function(I) sum(v[I]), 0))
  pchisq(sum(e^2), length(idx) - 2, lower.tail = FALSE)
}
loo_linear <- function(X, y, pen) {   # exact leave-one-out index of a linear smoother; intercept unpenalised
  X1 <- cbind(1, X); D <- diag(c(0, rep(pen, ncol(X))))
  A <- solve(crossprod(X1) + D, t(X1)); h <- rowSums(X1 * t(A)); fit <- drop(X1 %*% (A %*% y))
  (fit - h * y) / (1 - h)
}
std <- function(x) (x - mean(x)) / sd(x)
one <- function(r, d) {
  z <- gen(r, d); X <- z$X; y <- z$y; T <- function(xi) tryCatch(q3_S(std(xi), y), error = function(e) NA_real_)
  pI1 <- T(loo_linear(X, y, 4 * d$lamI))
  pI2 <- if (d$p < d$n - 2) T(loo_linear(X, y, 0)) else NA_real_
  pF <- if (!is.null(d$lambda)) T(cav_pieces(X, y, d$lambda)$xi) else NA_real_
  pF1 <- if (d$p < d$n) tryCatch(T(cav_pieces(X, y, 1)$xi), error = function(e) NA_real_) else NA_real_
  data.frame(rep = r, cell = d$name, I1 = pI1, I2 = pI2, full = pF, F1 = pF1)
}
if (PHASE == "smoke") {
  old <- read.csv(file.path(ROOT, "results", "E43_split_q3_perrep.csv"))
  d <- CELLS[[2]]; for (r in 1:3) { o <- one(r, d); print(c(o$I1, o$I2, o$full, old$full_q3[old$cell == d$name & old$rep == r])) }
  print(system.time(one(4, CELLS[[7]]))); quit("no")
}
done <- if (file.exists(OUTF)) { x <- read.csv(OUTF); paste(x$cell, x$rep) } else character()
cl <- makeCluster(W); clusterExport(cl, c("std", "gen", "one", "q3_S", "loo_linear", "a_for", "PHASE", "cf_expit", "statmod_gh", "cav_pieces", "cav_recal_fit", "cf_ridge"))
for (d in CELLS) for (st in seq(1, d$R, by = 200)) {
  rs <- st:min(st + 199, d$R); rs <- rs[!paste(d$name, rs) %in% done]; if (!length(rs)) next
  out <- do.call(rbind, parLapply(cl, rs, one, d = d))     # every batch written to disk: resume by re-running
  write.table(out, OUTF, sep = ",", row.names = FALSE, col.names = !file.exists(OUTF), append = file.exists(OUTF))
  cat(d$name, max(rs), "size I1", round(mean(out$I1 < .05, na.rm = TRUE), 3), "\n")
}
stopCluster(cl)
