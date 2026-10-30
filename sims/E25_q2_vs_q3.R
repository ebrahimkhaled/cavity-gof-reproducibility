# E25 -- two-parameter (logistic) vs three-parameter (logistic-normal) recalibration, under the null.
# (i) LARGE n with substantial index noise tau: p = 50, AR 0.9, heavy ridge (the penalty bends the coefficient direction away
#     from beta_0, so tau > 0 at small p) -- n = 2000, 5000, 20000. The q = 2 test should drift as n grows; q = 3 should not.
# (ii) small-sample cells (n = 400) where tau is weakly identified: does q = 3 keep its level?
# Also reports tau from Proposition 1 (full-fit beta in place of beta_(-i)) so the regime is documented.
# Usage: Rscript E25_q2_vs_q3.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 300; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20261105
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda, ar, norm = 2.31, prev = .5, alt = FALSE) {
  Sig <- ar^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0))
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sig = Sig, Sc = chol(Sig), prev = prev) }
CELLS <- list(mk("bent_n2000", 2000, 50, 2000 * 2, 0.9), mk("bent_n5000", 5000, 50, 5000 * 2, 0.9), mk("bent_n20000", 20000, 50, 20000 * 2, 0.9),
              mk("B_p100", 400, 100, 416, 0.7), mk("H5_k25_ar03", 400, 100, 416, 0.3), mk("prev08_k25", 400, 100, 416, 0.7, prev = .08))
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta0 <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta0, d$prev); y <- rbinom(d$n, 1, cf_expit(a + eta0))
  c2 <- cav_gof_fast(X, y, d$lambda, G = 10)
  c3 <- cav_lognormal(c2$xi, y, G = 10)
  X1 <- cbind(1, X); bh <- cf_ridge(X1, y, d$lambda, diag(c(0, rep(1, d$p))))
  s2 <- drop(t(d$b0) %*% d$Sig %*% d$b0); t2 <- drop(t(bh[-1]) %*% d$Sig %*% bh[-1]); cc <- drop(t(d$b0) %*% d$Sig %*% bh[-1])
  data.frame(rep = r, cell = d$name, q2_chi2 = c2$p_chi2, q2_exact = c2$p_a, q3_chi2 = c3$p_chi2, q3_exact = c3$p_a,
             S2 = c2$S, S3 = c3$S, tr3 = c3$tr, tau_fit = c3$tau, tau_prop1 = sqrt(max(s2 - cc^2 / t2, 0)), b_prop1 = cc / t2)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
f_out <- file.path(OUT, "E25_q2_vs_q3_perrep.csv"); if (file.exists(f_out)) file.remove(f_out)
for (d in CELLS) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
  if (is.null(o)) { cat("ALL FAILED", d$name, "\n"); next }
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("done %-14s n=%d  (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out); rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(q2_chi2, q2_exact, q3_chi2, q3_exact) ~ cell, res, rej)
agg$meanS2 <- aggregate(S2 ~ cell, res, mean)$S2; agg$meanS3 <- aggregate(S3 ~ cell, res, mean)$S3; agg$tr3 <- aggregate(tr3 ~ cell, res, mean)$tr3
agg$tau_prop1 <- aggregate(tau_prop1 ~ cell, res, mean)$tau_prop1; agg$tau_fit_med <- aggregate(tau_fit ~ cell, res, median)$tau_fit
agg$b_prop1 <- aggregate(b_prop1 ~ cell, res, mean)$b_prop1; agg$R <- aggregate(rep ~ cell, res, length)$rep
agg <- agg[match(sapply(CELLS, `[[`, "name"), agg$cell), ]
write.csv(agg, file.path(OUT, "E25_q2_vs_q3_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E25 | q = 2 vs q = 3 recalibration (null) | %.1f min ===\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
