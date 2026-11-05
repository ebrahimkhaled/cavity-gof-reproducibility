# E38 -- size of the THREE-PARAMETER (q = 3, logistic-normal) cavity test with the chi2_{G-2} reference, next to the
# two-parameter test (q = 2), in ALL 32 designs of Figure 3 (E18: 11 ridge designs, E19: 13 stress designs, E22: 8 lasso /
# elastic-net designs), on the IDENTICAL data sets (referee items M2 / S2, m8 / S12).
# Data are regenerated with each original script's generator and seed (Mersenne-Twister, set.seed(SEED + r), as those
# scripts ran: RNGkind is reset to that default before every set.seed). The q = 2 P-value is recomputed and checked against
# the stored one (fingerprint) in E38_summarise.R. Replications: as the originals (500; 300/200/300 in the three largest
# E19 designs), and extended to 2000 in the three key designs the referee named (B p/n 0.25, B p/n 0.40, lasso p = 2n).
# Every batch is appended to disk; a rerun skips (cell, rep) pairs already written.
# Usage: Rscript E38_q3_size_all.R [workers] [phase: base | extend | all]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) as.integer(args[1]) else 14; PHASE <- if (length(args) >= 2) args[2] else "all"
SMOKE <- PHASE == "smoke"     # smoke: 3 replicates of every cell, written to a scratch file
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
# the intercept for a target prevalence, with each original script's bracket (E18: -12..12, E19 / E22: -15..15), so that
# uniroot returns the identical root and the identical responses
a_for <- function(eta, target, src) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, if (src == "E18") c(-12, 12) else c(-15, 15))$root

# ---- designs, copied from E18 / E19 / E22 line for line --------------------------------------------------------------
mk18 <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, law = "gaussian", prev = 0.5, G = 10, indep = FALSE) {
  if (indep) { Sig <- diag(p); b0 <- c(.5, -.4, .3, -.3, .2) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  list(src = "E18", seed = 20261030, name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), law = law, prev = prev, G = G, R = 500) }
mk19 <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, law = "gaussian", prev = 0.5, G = 10, indep = FALSE, R = 500) {
  if (indep) { Sig <- diag(p); b0 <- rep_len(c(.5, -.4, .3, -.3, .2), p) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  list(src = "E19", seed = 20261031, name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), law = law, prev = prev, G = G, R = R) }
mk22 <- function(name, n, p, lambda_g, alpha = 1, ar = 0.7, norm = 2.31, sparse = 0, prev = .5) {
  Sig <- ar^abs(outer(1:p, 1:p, "-"))
  b0 <- if (sparse > 0) c(rep_len(c(.35, -.3, .25, -.2, .15), sparse), rep(0, p - sparse)) else rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0))
  list(src = "E22", seed = 20261101, name = name, n = n, p = p, lambda_g = lambda_g, alpha = alpha, b0 = b0, Sc = chol(Sig), law = "gaussian",
       prev = prev, G = 10, R = 500) }
CELLS <- list(
  mk18("A_p5", 500, 5, 100, indep = TRUE), mk18("B_p40", 400, 40, 416), mk18("B_p100", 400, 100, 416), mk18("B_p160", 400, 160, 416),
  mk18("H3_k40_lam50", 400, 160, 50), mk18("H5_k25_ar03", 400, 100, 416, ar = 0.3), mk18("H8_n800_k25", 800, 200, 832),
  mk18("k10_t3", 400, 40, 416, law = "t3"), mk18("prev15_k25", 400, 100, 416, prev = .15), mk18("prev08_k25", 400, 100, 416, prev = .08),
  mk18("prev08_A", 500, 5, 100, indep = TRUE, prev = .08),
  mk19("n1600_k25", 1600, 400, 1664, R = 300), mk19("n3200_k25", 3200, 800, 3328, R = 200),
  mk19("strong4_k25", 400, 100, 416, norm = 4), mk19("weak1_k25", 400, 100, 416, norm = 1),
  mk19("lam10_k25", 400, 100, 10), mk19("lam5000_k25", 400, 100, 5000), mk19("k50", 400, 200, 416),
  mk19("binary_k25", 400, 100, 416, law = "binary"), mk19("t3_k25", 400, 100, 416, law = "t3"),
  mk19("G5_k25", 400, 100, 416, G = 5), mk19("G20_k25", 400, 100, 416, G = 20),
  mk19("prev05_strong3_k25", 400, 100, 416, norm = 3, prev = .05), mk19("A_n5000", 5000, 5, 1000, indep = TRUE, R = 300),
  mk22("lasso_dense_k25_l02", 400, 100, 0.02), mk22("lasso_dense_k25_l05", 400, 100, 0.05),
  mk22("enet5_dense_k25_l05", 400, 100, 0.05, alpha = .5),
  mk22("lasso_sparse5_k25", 400, 100, 0.03, sparse = 5), mk22("lasso_sparse5_k50", 400, 200, 0.03, sparse = 5),
  mk22("lasso_sparse5_p2n", 300, 600, 0.04, sparse = 5), mk22("enet5_sparse10_p2n", 300, 600, 0.04, alpha = .5, sparse = 10),
  mk22("lasso_sparse5_prev10", 400, 100, 0.02, sparse = 5, prev = .10))
KEY <- c("B_p100", "B_p160", "lasso_sparse5_p2n"); R_KEY <- 2000

draw <- function(d, r) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(d$seed + r)
  if (d$src == "E22") X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc else {
    Zm <- switch(d$law, gaussian = matrix(rnorm(d$n * d$p), d$n), t3 = matrix(rt(d$n * d$p, 3) / sqrt(3), d$n),
                 binary = matrix((rbinom(d$n * d$p, 1, .3) - .3) / sqrt(.21), d$n))
    X <- Zm %*% d$Sc }
  eta <- drop(X %*% d$b0); a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta, d$prev, d$src)
  list(X = X, y = rbinom(d$n, 1, cf_expit(a + eta)))
}
glmnet_xi <- function(X, y, lambda_g, alpha) {        # the index of cav_gof_glmnet, returned
  n <- nrow(X); fit <- glmnet::glmnet(X, y, family = "binomial", alpha = alpha, lambda = lambda_g, standardize = FALSE, thresh = 1e-12, maxit = 1e6)
  b <- as.numeric(coef(fit, s = lambda_g)); E <- which(b[-1] != 0)
  X1 <- cbind(1, X); eta <- drop(X1 %*% b); pi <- pmin(pmax(cf_expit(eta), 1e-9), 1 - 1e-9); w <- pi * (1 - pi)
  XE <- cbind(1, X[, E, drop = FALSE]); DE <- diag(c(0, rep(1, length(E))), length(E) + 1)
  ME <- crossprod(XE, w * XE) + n * lambda_g * (1 - alpha) * DE
  q <- rowSums((XE %*% solve(ME)) * XE); H <- w * q
  eta - q / (1 - H) * (y - pi)
}
one <- function(r, d) {
  z <- draw(d, r)
  xi <- if (d$src == "E22") glmnet_xi(z$X, z$y, d$lambda_g, d$alpha) else cav_gof_fast(z$X, z$y, d$lambda, G = d$G)$xi
  c2 <- cav_from_index(xi, z$y, d$G)
  c3 <- tryCatch(cav_lognormal(xi, z$y, d$G), error = function(e) NULL)
  data.frame(rep = r, cell = d$name, src = d$src, G = d$G, q2_chi2 = c2$p_chi2, q2_exact = c2$p_a, S2 = c2$S,
             q3_chi2G2 = if (is.null(c3)) NA else pchisq(c3$S, d$G - 2, lower.tail = FALSE),
             q3_exact = if (is.null(c3)) NA else c3$p_a, S3 = if (is.null(c3)) NA else c3$S, tr3 = if (is.null(c3)) NA else c3$tr,
             tau3 = if (is.null(c3)) NA else c3$tau, b2 = c2$b)
}
f_out <- if (SMOKE) file.path(tempdir(), "E38_smoke.csv") else file.path(OUT, "E38_q3_size_all_perrep.csv")
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c("one", "draw", "glmnet_xi", "a_for", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(CompQuadForm); library(glmnet) })))
run_cell <- function(d, Rmax) {
  have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$cell == d$name] } else integer(0)
  todo <- setdiff(seq_len(Rmax), have); if (!length(todo)) return(invisible())
  for (batch in split(todo, ceiling(seq_along(todo) / (4 * W)))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
    if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  }
  cat(sprintf("done %-22s up to rep %d  (%.1f min)\n", d$name, Rmax, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (SMOKE) { for (d in CELLS) run_cell(d, 3); print(read.csv(f_out)[, 1:9]) }
if (PHASE %in% c("base", "all")) for (d in CELLS) run_cell(d, d$R)
if (PHASE %in% c("extend", "all")) for (d in CELLS) if (d$name %in% KEY) run_cell(d, R_KEY)
stopCluster(cl)
cat(sprintf("E38 finished (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
