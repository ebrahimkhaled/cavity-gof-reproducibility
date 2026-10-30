# E31 -- the SAMPLE-SPLIT cavity test (fully proved: no dependence, exact curve) and its CROSS-FIT version (swap halves,
# Cauchy combination), against the full leave-one-out cavity test, on identical data. Size in the hardest designs and
# power against the departures where the cavity test had power (quadratic along the index; cubic at fixed p).
# Usage: Rscript E31_split_crossfit.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261111
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
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
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(d$dep, null = 0, cubic = d$sd_eta * d$rel * (Zs^3 - 3 * Zs) / sqrt(6), quadratic = d$sd_eta * d$rel * (Zs^2 - 1) / sqrt(2))
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta, d$prev); y <- rbinom(d$n, 1, cf_expit(a + eta))
  full <- if (is.null(d$lambda_g)) cav_gof_fast(X, y, d$lambda, G = 10)$p_chi2 else cav_gof_glmnet(X, y, d$lambda_g, G = 10)$p_chi2
  sp <- cav_split_gof(X, y, lambda = d$lambda, lambda_g = d$lambda_g, G = 10, seed = SEED + 7e5 + r)
  data.frame(rep = r, cell = d$name, full = full, split = sp$p_split, cross = sp$p_cross)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(glmnet)})))
res <- do.call(rbind, lapply(CELLS, function(d) { o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
  cat(sprintf("done %-18s n=%d (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl); write.csv(res, file.path(OUT, "E31_split_crossfit_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(full, split, cross) ~ cell, res, rej); agg$R <- aggregate(rep ~ cell, res, length)$rep
agg <- agg[match(sapply(CELLS, `[[`, "name"), agg$cell), ]
write.csv(agg, file.path(OUT, "E31_split_crossfit_summary.csv"), row.names = FALSE)
cat(sprintf("\n=== E31 | split / cross-fit / full leave-one-out cavity test | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
