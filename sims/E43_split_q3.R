# E43 -- the sample-split and cross-fitted cavity tests with the THREE-PARAMETER recalibration (q = 3) and the chi2_{G-2}
# reference -- the version Theorem 2 covers -- next to the two-parameter version (q = 2) that E31 reported, on E31's
# IDENTICAL data and splits (referee item m7 / S11). Optional phase "ads": the same for the Internet-Ads split placebo
# (E35c's data and splits).
# Data: E31's generator, Mersenne-Twister / Inversion / Rejection reset before set.seed(20261111 + r) (E31 ran with that
# default); split seed 20261111 + 7e5 + r, exactly as cav_split_gof() draws it. Each half's index is computed once and both
# recalibrations are applied to it. q = 3 uses cav_lognormal()'s fit (same start, optimiser and bounds) and computes only S
# (the E42 route), so a singular information matrix cannot stop it; chi2_{G-2} needs no information matrix.
# Fingerprint: the q = 2 full / split / cross P-values are checked against results/E31_split_crossfit_perrep.csv.
# Appended per batch, resumable. Usage: Rscript E43_split_q3.R [workers] [phase: e31 | ads | smoke] [R_ads]
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) as.integer(args[1]) else 12; PHASE <- if (length(args) >= 2) args[2] else "e31"
R_ADS <- if (length(args) >= 3) as.integer(args[3]) else 120
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- if (PHASE == "smoke") tempdir() else file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
SEED <- 20261111; R <- if (PHASE == "smoke") 2 else 500

# ---- E31's designs, copied line for line ---------------------------------------------------------------------------
a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda = NULL, ar = 0.7, norm = 2.31, prev = .5, indep = FALSE, sparse = 0, lambda_g = NULL, dep = "null", rel = 0) {
  if (indep) { Sig <- diag(p); b0 <- c(.5, -.4, .3, -.3, .2) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-"))
    b0 <- if (sparse > 0) c(rep_len(c(.35, -.3, .25, -.2, .15), sparse), rep(0, p - sparse)) else rep_len(c(.35, -.3, .25, -.2, .15), p)
    b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), prev = prev, lambda_g = lambda_g, dep = dep, rel = rel,
       sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0))) }
CELLS <- list(
  mk("A_p5", 500, 5, 100, indep = TRUE), mk("B_p40", 400, 40, 416), mk("B_p100", 400, 100, 416), mk("B_p160", 400, 160, 416),
  mk("H5_k25_ar03", 400, 100, 416, ar = .3), mk("prev08_k25", 400, 100, 416, prev = .08),
  mk("lasso_sparse5_p2n", 300, 600, sparse = 5, lambda_g = 0.04),
  mk("A_cubic1", 500, 5, 100, indep = TRUE, dep = "cubic", rel = 1), mk("B_p40_quad1", 400, 40, 416, dep = "quadratic", rel = 1),
  mk("B_p100_quad1", 400, 100, 416, dep = "quadratic", rel = 1))

# ---- q = 3 statistic only (cav_lognormal's fit, no information matrix) -------------------------------------------------
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
cauchy <- function(p) 0.5 - atan(mean(tan((0.5 - p) * pi))) / pi
# split and cross-fit, both recalibrations, on one split (as cav_split_gof draws it)
split_both <- function(X, y, lambda, lambda_g, seed, alpha = 1, G = 10) {
  n <- nrow(X)
  fit_fun <- if (is.null(lambda_g)) function(Xa, ya) cf_ridge(cbind(1, Xa), ya, lambda * nrow(Xa) / n, diag(c(0, rep(1, ncol(Xa))))) else
    function(Xa, ya) as.numeric(coef(glmnet::glmnet(Xa, ya, family = "binomial", alpha = alpha, lambda = lambda_g, standardize = FALSE,
                                                    thresh = 1e-12, maxit = 1e6), s = lambda_g))
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed)
  idxA <- sort(sample.int(n, floor(n / 2))); idxB <- setdiff(seq_len(n), idxA)
  half <- function(train, test) { xi <- cav_split_index(X, y, train, fit_fun)
    c(q2 = cav_from_index(xi, y[test], G)$p_chi2, q3 = q3_S(xi, y[test], G)) }
  h1 <- half(idxA, idxB); h2 <- half(idxB, idxA)
  c(split_q2 = h1[["q2"]], cross_q2 = cauchy(c(h1[["q2"]], h2[["q2"]])), split_q3 = h1[["q3"]], cross_q3 = cauchy(c(h1[["q3"]], h2[["q3"]])))
}
glmnet_xi <- function(X, y, lambda_g) {   # cav_gof_glmnet's index (lasso)
  n <- nrow(X); fit <- glmnet::glmnet(X, y, family = "binomial", alpha = 1, lambda = lambda_g, standardize = FALSE, thresh = 1e-12, maxit = 1e6)
  b <- as.numeric(coef(fit, s = lambda_g)); E <- which(b[-1] != 0)
  X1 <- cbind(1, X); eta <- drop(X1 %*% b); pi <- pmin(pmax(cf_expit(eta), 1e-9), 1 - 1e-9); w <- pi * (1 - pi)
  XE <- cbind(1, X[, E, drop = FALSE]); q <- rowSums((XE %*% solve(crossprod(XE, w * XE))) * XE); H <- w * q
  eta - q / (1 - H) * (y - pi)
}
one <- function(r, d) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED + r)
  X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(d$dep, null = 0, cubic = d$sd_eta * d$rel * (Zs^3 - 3 * Zs) / sqrt(6), quadratic = d$sd_eta * d$rel * (Zs^2 - 1) / sqrt(2))
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta, d$prev); y <- rbinom(d$n, 1, cf_expit(a + eta))
  xi <- if (is.null(d$lambda_g)) cav_gof_fast(X, y, d$lambda, G = 10)$xi else glmnet_xi(X, y, d$lambda_g)
  sp <- split_both(X, y, d$lambda, d$lambda_g, SEED + 7e5 + r)
  data.frame(rep = r, cell = d$name, full_q2 = cav_from_index(xi, y, 10)$p_chi2, full_q3 = q3_S(xi, y, 10), t(sp))
}
# ---- Internet Ads (E35c's data and splits) ---------------------------------------------------------------------------
if (PHASE == "ads") {
  SEED_A <- 20261133
  dd <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
  dd <- dd[!is.na(dd$V4), ]; Xa <- as.matrix(dd[, 4:(ncol(dd) - 1)]); ya <- as.numeric(dd[, ncol(dd)] == "ad.")
  Xa <- Xa[, colSums(Xa) >= 5]; na <- nrow(Xa); pa <- ncol(Xa)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED_A)
  cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10); lam_a <- na * cvf$lambda.1se
  ba <- cf_ridge(cbind(1, Xa), ya, lam_a, diag(c(0, rep(1, pa))))
  one_ads <- function(r) { RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED_A + 1000 + r)
    y <- rbinom(na, 1, cf_expit(drop(cbind(1, Xa) %*% ba))); data.frame(rep = r, t(split_both(Xa, y, lam_a, NULL, SEED_A + r))) }
}
cl <- makeCluster(W)
clusterExport(cl, c("one", "split_both", "q3_S", "cauchy", "glmnet_xi", "a_for", "SEED", ls(pattern = "^(cf_|cav_|naive_|statmod_)"),
                    if (PHASE == "ads") c("one_ads", "Xa", "ba", "lam_a", "na", "SEED_A")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(CompQuadForm); library(glmnet) })))
t0 <- Sys.time()
if (PHASE %in% c("e31", "smoke")) {
  f_out <- file.path(OUT, "E43_split_q3_perrep.csv")
  for (d in CELLS) {
    have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$cell == d$name] } else integer(0); todo <- setdiff(seq_len(R), have)
    for (batch in split(todo, ceiling(seq_along(todo) / (4 * W)))) {
      o <- do.call(rbind, parLapply(cl, batch, function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
      if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out)) }
    cat(sprintf("done %-18s (%.1f min)\n", d$name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  if (PHASE == "smoke") print(read.csv(f_out))
}
if (PHASE == "ads") {
  f_out <- file.path(OUT, "E43_ads_split_q3_perrep.csv")
  have <- if (file.exists(f_out)) read.csv(f_out)$rep else integer(0); todo <- setdiff(seq_len(R_ADS), have)
  for (batch in split(todo, ceiling(seq_along(todo) / W))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r) tryCatch(one_ads(r), error = function(e) NULL)))
    if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
    cat(sprintf("ads reps %d-%d (%.1f min)\n", min(batch), max(batch), as.numeric(difftime(Sys.time(), t0, units = "mins")))) }
}
stopCluster(cl); cat(sprintf("E43 %s finished (%.1f min)\n", PHASE, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
