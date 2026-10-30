# E29 -- why does the THREE-parameter cavity test keep a small excess at n = 20000 (bent design, tau ~ .75)?
# On the same data: S with the TRUE Proposition-2 curve (theta from the full-fit beta_hat with the true beta_0, Sigma; no
# estimation, reference = identity G), S with the fitted logistic-normal (q = 3), S with the fitted logistic (q = 2).
# If S_true is centred at G = 10 and S_q3 is not, the excess comes from ESTIMATING the curve; if S_true is also high, it comes
# from the grouping or from the dependence between cavity indices. Also stores the decile means of y - Lambda_true.
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 400; W <- if (length(args) >= 2) as.integer(args[2]) else 12
n <- if (length(args) >= 3) as.integer(args[3]) else 20000
SEED <- 20261109
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
p <- 50; Sg <- 0.9^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(drop(t(b0) %*% Sg %*% b0)); Cg <- chol(Sg)
gh <- statmod_gh(40)
Lam <- function(xi, a, b, tau) { P <- 0; for (k in seq_along(gh$x)) P <- P + gh$w[k] * cf_expit(a + b * xi + tau * gh$x[k]); P }
one <- function(r) {
  set.seed(SEED + r); X <- matrix(rnorm(n * p), n) %*% Cg; y <- rbinom(n, 1, cf_expit(drop(X %*% b0)))
  c2 <- cav_gof_fast(X, y, 2 * n, G = 10); c3 <- cav_lognormal(c2$xi, y, G = 10)
  X1 <- cbind(1, X); bh <- cf_ridge(X1, y, 2 * n, diag(c(0, rep(1, p))))
  s2 <- drop(t(b0) %*% Sg %*% b0); t2 <- drop(t(bh[-1]) %*% Sg %*% bh[-1]); cc <- drop(t(b0) %*% Sg %*% bh[-1])
  b <- cc / t2; tau <- sqrt(max(s2 - cc^2 / t2, 0)); a <- -b * bh[1]
  m <- Lam(c2$xi, a, b, tau); v <- m * (1 - m)
  idx <- split(seq_len(n), c2$groups)
  e_true <- vapply(idx, function(I) sum(y[I] - m[I]) / sqrt(sum(v[I])), 0.0)
  data.frame(rep = r, S_true = sum(e_true^2), S3 = c3$S, tr3 = c3$tr, p3 = c3$p_a, S2 = c2$S, tr2 = c2$tr_a, p2 = c2$p_a,
             tau_true = tau, tau_fit = c3$tau, b_true = b, b_fit3 = c3$b, t(setNames(e_true, paste0("e", 1:10))))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r), error = function(e) NULL)))
stopCluster(cl); write.csv(res, file.path(OUT, sprintf("E29_q3_excess_n%d_perrep.csv", n)), row.names = FALSE)
se <- function(v) sd(v) / sqrt(length(v))
cat(sprintf("=== E29 | n = %d | R = %d | %.1f min ===\n", n, nrow(res), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("S with the TRUE curve (no estimation; reference G = 10): mean %.3f (SE %.3f)\n", mean(res$S_true), se(res$S_true)))
cat(sprintf("S with fitted q = 3: mean %.3f (SE %.3f) vs trace %.3f  | size %.3f\n", mean(res$S3), se(res$S3), mean(res$tr3), mean(res$p3 < .05)))
cat(sprintf("S with fitted q = 2: mean %.3f (SE %.3f) vs trace %.3f  | size %.3f\n", mean(res$S2), se(res$S2), mean(res$tr2), mean(res$p2 < .05)))
cat(sprintf("tau true %.3f vs fitted %.3f (median) | b true %.3f vs fitted %.3f\n", mean(res$tau_true), median(res$tau_fit), mean(res$b_true), mean(res$b_fit3)))
ec <- grep("^e[0-9]+$", names(res)); mu <- colMeans(res[, ec])
cat("decile means of the TRUE-curve residual:", paste(sprintf("%+.3f", mu), collapse = " "), sprintf("| ||mean||^2 = %.3f\n", sum(mu^2)))
cat("covariance eigenvalues (true curve):", paste(sprintf("%.2f", sort(eigen(cov(res[, ec]), only.values = TRUE)$values, decreasing = TRUE)), collapse = " "), "\n")
