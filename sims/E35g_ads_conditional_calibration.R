# E35g -- size of the design-conditional (simulation-calibrated) cavity test on the Internet-Ads predictors.
# E35 showed that the closed-form reference fails on this fixed design (0.57), and E35b-E35f that the cause is a
# systematic misfit of the logistic curve in xi for this design, not a variance problem. The remedy tested here:
# refer S to its distribution under the fitted model, by simulating y* from the ridge fit on the observed X, refitting
# and recomputing S (B times). Outer loop: placebo data sets as in E35 (same seeds), each tested this way.
# Usage: Rscript E35g_ads_conditional_calibration.R [R_outer] [B_inner] [workers]
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 60; B <- if (length(args) >= 2) as.integer(args[2]) else 39
W <- if (length(args) >= 3) as.integer(args[3]) else 12
SEED <- 20261133
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

d <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
d <- d[!is.na(d$V4), ]
Xa <- as.matrix(d[, 4:(ncol(d) - 1)]); ya <- as.numeric(d[, ncol(d)] == "ad.")
Xa <- Xa[, colSums(Xa) >= 5]; na <- nrow(Xa); pa <- ncol(Xa)
set.seed(SEED); cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
lam_a <- na * cvf$lambda.1se
X1a <- cbind(1, Xa); Dpen <- diag(c(0, rep(1, pa))); ba <- cf_ridge(X1a, ya, lam_a, Dpen)

newton <- function(y, beta, tol = 1e-9, maxit = 50) {      # warm-started ridge Newton on the Internet-Ads design
  for (it in seq_len(maxit)) {
    eta <- drop(X1a %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1a, w * X1a) + lam_a * Dpen, crossprod(X1a, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn)
    beta <- bn
  }
  beta
}
stat <- function(y, beta0) {                               # one fit -> cavity statistic S (two-parameter recalibration)
  bh <- newton(y, beta0); eta <- drop(X1a %*% bh); pi <- cf_expit(eta); w <- pi * (1 - pi)
  q <- rowSums((X1a %*% solve(crossprod(X1a, w * X1a) + lam_a * Dpen)) * X1a)
  list(S = cav_from_index(eta - q / (1 - w * q) * (y - pi), y)$S, beta = bh)
}
outer_rep <- function(r) {
  set.seed(SEED + 1000 + r); y <- rbinom(na, 1, cf_expit(drop(X1a %*% ba)))   # the E35 placebo data set r
  s0 <- stat(y, ba)
  set.seed(SEED + 5000 + r)
  Sb <- vapply(seq_len(B), function(b) { yb <- rbinom(na, 1, cf_expit(drop(X1a %*% s0$beta))); stat(yb, s0$beta)$S }, 0.0)
  data.frame(rep = r, S = s0$S, p_chi2 = pchisq(s0$S, 8, lower.tail = FALSE), p_boot = (1 + sum(Sb >= s0$S)) / (B + 1),
             Sb_mean = mean(Sb), Sb_sd = sd(Sb))
}

f_out <- file.path(OUT, "E35g_ads_conditional_perrep.csv")
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c("X1a", "ba", "lam_a", "Dpen", "na", "B", "SEED", "newton", "stat", "outer_rep", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
have <- if (file.exists(f_out)) read.csv(f_out)$rep else integer(0)
for (batch in split(setdiff(seq_len(R), have), ceiling(seq_along(setdiff(seq_len(R), have)) / W))) {
  o <- do.call(rbind, parLapply(cl, batch, function(r) tryCatch(outer_rep(r), error = function(e) NULL)))
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("reps %d-%d  (%.1f min)\n", min(batch), max(batch), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out); a <- mean(res$p_boot <= .05)
cat(sprintf("\n=== E35g | design-conditional cavity test, Internet Ads, %d placebo data sets, B = %d ===\n", nrow(res), B))
cat(sprintf("size: closed-form chi2 %.3f, simulation-calibrated %.3f (SE %.3f); mean simulated S %.2f (sd %.2f)\n",
            mean(res$p_chi2 < .05), a, sqrt(a * (1 - a) / nrow(res)), mean(res$Sb_mean), mean(res$Sb_sd)))
