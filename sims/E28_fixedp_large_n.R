# E28 -- is the excess mean of S at fixed p and large n (E19 A_n5000: 8.69 vs trace 8.03, R = 300) real?
# Fresh seeds, R = 1000, design A (p = 5 independent, lambda = n / 5), n = 500, 5000, 20000. If real, it grows with n.
# Also reports the mean of S for the TEXTBOOK statistic on the same data (as a check of the grouping machinery itself).
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20261108
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
b0 <- c(.5, -.4, .3, -.3, .2)
one <- function(r, n) {
  set.seed(SEED + r + n); X <- matrix(rnorm(n * 5), n); y <- rbinom(n, 1, cf_expit(drop(X %*% b0)))
  c2 <- cav_gof_fast(X, y, n / 5, G = 10); c3 <- cav_lognormal(c2$xi, y, G = 10)
  data.frame(rep = r, n = n, p2 = c2$p_a, S2 = c2$S, tr2 = c2$tr_a, b = c2$b, p3 = c3$p_a, S3 = c3$S, tr3 = c3$tr, tau = c3$tau)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(c(500, 5000, 20000), function(n) { o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, n) tryCatch(one(r, n), error = function(e) NULL), n = n))
  cat(sprintf("done n=%d (%.1f min)\n", n, as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl); write.csv(res, file.path(OUT, "E28_fixedp_large_n_perrep.csv"), row.names = FALSE)
for (n in unique(res$n)) { s <- res[res$n == n, ]
  cat(sprintf("n %6d R %d | q2: size %.3f mean S %.2f (SE %.2f) trace %.2f | q3: size %.3f mean S %.2f trace %.2f | median tau %.3f | slope %.2f\n",
              n, nrow(s), mean(s$p2 < .05), mean(s$S2), sd(s$S2) / sqrt(nrow(s)), mean(s$tr2), mean(s$p3 < .05), mean(s$S3), mean(s$tr3), median(s$tau), mean(s$b))) }
