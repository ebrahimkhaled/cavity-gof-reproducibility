# E23 -- direct numerical verification of Proposition 1 (exact cavity calibration law) and of the ALO index.
# Design B (n 400, p 100, AR 0.7, lambda 416) and the hardest held-out cell (AR 0.3). Per replicate, for 25 random
# observations i: refit the ridge WITHOUT i (exact leave-one-out), record the exact cavity index, the ALO index,
# y_i, and the Proposition-1 parameters (a, b, tau) computed from beta_hat_(-i) with the TRUE beta_0 and Sigma.
# Checks: (1) ALO vs exact error; (2) mean of y_i - Lambda(xi_i; theta_(-i)) = 0; (3) binned calibration curve.
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 200; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- as.integer(if (length(args) >= 3) args[3] else 20261102); KI <- 25
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
mk <- function(name, n, p, lambda, ar) { Sig <- ar^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * 2.31 / sqrt(drop(t(b0) %*% Sig %*% b0)); list(name = name, n = n, p = p, lambda = lambda, Sig = Sig, b0 = b0, Sc = chol(Sig)) }
CELLS <- list(mk("B_p100", 400, 100, 416, 0.7), mk("H5_ar03", 400, 100, 416, 0.3))
gh <- statmod_gh(40)
Lam <- function(xi, a, b, tau) { P <- 0; for (k in seq_along(gh$x)) P <- P + gh$w[k] * cf_expit(a + b * xi + tau * gh$x[k]); P }
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; pt <- cf_expit(drop(X %*% d$b0)); y <- rbinom(d$n, 1, pt)
  X1 <- cbind(1, X); Dp <- diag(c(0, rep(1, d$p))); bh <- cf_ridge(X1, y, d$lambda, Dp)
  cv <- cav_gof_fast(X, y, d$lambda, beta = bh)
  s2 <- drop(t(d$b0) %*% d$Sig %*% d$b0)
  ii <- sample.int(d$n, KI)
  do.call(rbind, lapply(ii, function(i) {
    bi <- cf_ridge(X1[-i, ], y[-i], d$lambda, Dp)                  # exact leave-one-out refit
    t2 <- drop(t(bi[-1]) %*% d$Sig %*% bi[-1]); cc <- drop(t(d$b0) %*% d$Sig %*% bi[-1])
    b <- cc / t2; tau <- sqrt(max(s2 - cc^2 / t2, 0)); a <- -b * bi[1]
    xe <- drop(X1[i, ] %*% bi)
    data.frame(rep = r, cell = d$name, i = i, y = y[i], pi_true = pt[i], xi_exact = xe, xi_alo = cv$xi[i], eta_hat = drop(X1[i, ] %*% bh),
               a = a, b = b, tau = tau, Lam_exact = Lam(xe, a, b, tau), m_fit = cv$m[i], b_fit = cv$b)
  }))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))))
stopCluster(cl)
write.csv(res, file.path(OUT, sprintf("E23_verify_prop1_perobs_seed%d.csv", SEED)), row.names = FALSE)
cat(sprintf("=== E23 | Proposition 1 + ALO check | R = %d x %d obs | %.1f min ===\n", R, KI, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
for (cc in unique(res$cell)) {
  s <- res[res$cell == cc, ]
  err <- s$xi_alo - s$xi_exact; self <- s$eta_hat - s$xi_exact
  cat(sprintf("\n[%s]  n_obs = %d\n", cc, nrow(s)))
  cat(sprintf("  ALO error: sd %.4f, max |.| %.4f   vs self-influence eta_hat - xi_exact: sd %.4f   (sd xi %.3f)\n", sd(err), max(abs(err)), sd(self), sd(s$xi_exact)))
  res1 <- s$y - s$Lam_exact
  cat(sprintf("  mean(y - Lambda(xi_exact; theta_-i)) = %+.4f (SE %.4f)   mean(y - m_fit) = %+.4f (SE %.4f)\n",
              mean(res1), sd(res1) / sqrt(nrow(s)), mean(s$y - s$m_fit), sd(s$y - s$m_fit) / sqrt(nrow(s))))
  cat(sprintf("  Prop-1 parameters: b = %.2f (sd %.2f), tau = %.2f (sd %.2f); fitted logistic slope b_fit = %.2f\n", mean(s$b), sd(s$b), mean(s$tau), sd(s$tau), mean(s$b_fit)))
  br <- quantile(s$xi_exact, 0:10 / 10); g <- cut(s$xi_exact, br, include.lowest = TRUE)
  tab <- data.frame(decile = 1:10, obs = tapply(s$y, g, mean), prop1 = tapply(s$Lam_exact, g, mean), fit = tapply(s$m_fit, g, mean), n = as.vector(table(g)))
  tab$z <- (tab$obs - tab$prop1) / sqrt(tab$prop1 * (1 - tab$prop1) / tab$n)
  print(round(tab, 3), row.names = FALSE)
}
