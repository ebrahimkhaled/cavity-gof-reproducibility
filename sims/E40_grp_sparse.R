# E40 -- the residual-prediction test of Jankova, Shah, Buhlmann and Samworth (2020) (GRPtests::GRPtest, defaults:
# penalised version, random-forest residual prediction, 5 splits) against the cavity test in a SPARSE design, where
# GRP's theory applies (referee item M7(b) / S4(b)); in E33 (dense coefficients) GRP's size was 0.005 and 0.33.
# Designs (n = 400, s = 5 active coefficients scaled so that the true index has SD 2.31):
#   indep : p = 100 independent N(0, 1) predictors        AR07 : p = 100, autoregressive correlation 0.7 (as E22)
# Departures along the index (quadratic, size 1.5) and off the index (interaction x1 x2 of two ACTIVE predictors,
# size 1.5), defined as in E30 / the paper (true index sd s; Hermite quadratic r s (Z^2 - 1)/sqrt 2; x1 x2 scaled to sd r s).
# Tests on the same data: cavity after the lasso (glmnet lambda 0.03, the E22 sparse value; q = 2 chi2_8 and q = 3
# chi2_8), cavity after the ridge (lambda 416; q = 2, q = 3), GRPtest, textbook HL on the lasso fit.
# Phase "pilot": null only, 60 replicates, seeds 20261150 + 1e5 + r (to check GRP's size before the main run).
# Phase "main" : null, quadratic, interaction; R replicates (default 400); seeds 20261150 + r.
# RNGkind("L'Ecuyer-CMRG") is reset before every set.seed. Appended per batch; a rerun skips (cell, rep) pairs on disk.
# Usage: Rscript E40_grp_sparse.R [workers] [phase] [R]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) as.integer(args[1]) else 10; PHASE <- if (length(args) >= 2) args[2] else "main"
R <- if (length(args) >= 3) as.integer(args[3]) else 400
SEED <- 20261150; LAMG <- 0.03; LAMR <- 416
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
mk <- function(name, ar, dep = "null", rel = 0, n = 400, p = 100, s = 5) {
  Sig <- if (ar > 0) ar^abs(outer(1:p, 1:p, "-")) else diag(p)
  b0 <- c(rep_len(c(.35, -.3, .25, -.2, .15), s), rep(0, p - s)); b0 <- b0 * 2.31 / sqrt(drop(t(b0) %*% Sig %*% b0))
  list(name = name, n = n, p = p, Sc = chol(Sig), b0 = b0, sd_eta = 2.31, sd_x1x2 = sqrt(1 + Sig[1, 2]^2), dep = dep, rel = rel) }
CELLS <- list(mk("indep_null", 0), mk("indep_quad", 0, "quadratic", 1.5), mk("indep_int", 0, "interaction", 1.5),
              mk("AR07_null", .7), mk("AR07_quad", .7, "quadratic", 1.5), mk("AR07_int", .7, "interaction", 1.5))
# pilot 2 (added after pilot 1 gave GRP size 0/60 and median P near 1 at n = 400): larger n
PILOT2 <- list(mk("indep_n1000_null", 0, n = 1000), mk("indep_n1000_p20_null", 0, n = 1000, p = 20))
glmnet_fit_xi <- function(X, y, lambda_g) {
  n <- nrow(X); fit <- glmnet::glmnet(X, y, family = "binomial", alpha = 1, lambda = lambda_g, standardize = FALSE, thresh = 1e-12, maxit = 1e6)
  b <- as.numeric(coef(fit, s = lambda_g)); E <- which(b[-1] != 0)
  X1 <- cbind(1, X); eta <- drop(X1 %*% b); pi <- pmin(pmax(cf_expit(eta), 1e-9), 1 - 1e-9); w <- pi * (1 - pi)
  XE <- cbind(1, X[, E, drop = FALSE]); ME <- crossprod(XE, w * XE)
  q <- rowSums((XE %*% solve(ME)) * XE); H <- w * q
  list(xi = eta - q / (1 - H) * (y - pi), eta = eta, nnz = length(E))
}
both <- function(xi, y) { c2 <- cav_from_index(xi, y, 10); c3 <- tryCatch(cav_lognormal(xi, y, 10), error = function(e) NULL)
  c(q2 = c2$p_chi2, q3 = if (is.null(c3)) NA else pchisq(c3$S, 8, lower.tail = FALSE)) }
one <- function(r, d, seed0) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed0 + r)
  X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(d$dep, null = 0, quadratic = d$sd_eta * d$rel * (Zs^2 - 1) / sqrt(2),
                      interaction = (d$rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(d$n, 1, cf_expit(eta))
  la <- glmnet_fit_xi(X, y, LAMG); pl <- tryCatch(both(la$xi, y), error = function(e) c(q2 = NA, q3 = NA))
  pr <- tryCatch(both(cav_gof_fast(X, y, LAMR, G = 10)$xi, y), error = function(e) c(q2 = NA, q3 = NA))
  grp <- tryCatch(as.numeric(suppressWarnings(suppressMessages(GRPtests::GRPtest(X, y, fam = "binomial")))), error = function(e) NA_real_)
  data.frame(rep = r, cell = d$name, lasso_q2 = pl[["q2"]], lasso_q3 = pl[["q3"]], ridge_q2 = pr[["q2"]], ridge_q3 = pr[["q3"]],
             grp = grp, textbook = naive_hl(la$eta, y, 10), nnz = la$nnz)
}
cl <- makeCluster(W); clusterExport(cl, c("one", "both", "glmnet_fit_xi", "LAMG", "LAMR", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(CompQuadForm); library(glmnet); library(GRPtests) }); options(ranger.num.threads = 1) }))
t0 <- Sys.time()
run <- function(cells, Rn, seed0, f_out) for (d in cells) {
  have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$cell == d$name] } else integer(0)
  todo <- setdiff(seq_len(Rn), have)
  for (batch in split(todo, ceiling(seq_along(todo) / (2 * W)))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r, d, s0) tryCatch(one(r, d, s0), error = function(e) NULL), d = d, s0 = seed0))
    if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  }
  cat(sprintf("done %-11s (%.1f min)\n", d$name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (PHASE == "smoke") run(CELLS[c(1, 3)], 2, SEED + 2e5, file.path(tempdir(), "E40_smoke.csv"))
if (PHASE == "pilot") run(CELLS[c(1, 4)], 60, SEED + 1e5, file.path(OUT, "E40_grp_sparse_pilot_perrep.csv"))
if (PHASE == "pilot2") run(PILOT2, 40, SEED + 3e5, file.path(OUT, "E40_grp_sparse_pilot_perrep.csv"))
# main: GRPtest costs 13-25 s per call, so the main run uses the independent design (3 cells); "main_ar" adds AR07
if (PHASE == "main") run(CELLS[1:3], R, SEED, file.path(OUT, "E40_grp_sparse_perrep.csv"))
if (PHASE == "main_ar") run(CELLS[4:6], R, SEED, file.path(OUT, "E40_grp_sparse_perrep.csv"))
if (PHASE == "smoke") print(read.csv(file.path(tempdir(), "E40_smoke.csv")))
stopCluster(cl); cat(sprintf("E40 %s finished (%.1f min)\n", PHASE, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
