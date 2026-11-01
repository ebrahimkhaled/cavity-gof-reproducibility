# E35b -- why does the cavity test reject 57 per cent of correct models in the Internet-Ads placebo (E35)?
# Same placebo data sets as E35 (same seeds). For each replicate: the one-step cavity index, the exact ten-fold
# cross-fitted index (no approximation), and per-observation leverage and feature counts. For replicate 1 the
# exact leave-one-out index is also computed on 240 random observations, to measure the one-step error directly.
# Usage: Rscript E35b_ads_placebo_diagnosis.R [reps] [workers]
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 12; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261133
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

d <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
d <- d[!is.na(d$V4), ]
Xa <- as.matrix(d[, 4:(ncol(d) - 1)]); ya <- as.numeric(d[, ncol(d)] == "ad.")
Xa <- Xa[, colSums(Xa) >= 5]
na <- nrow(Xa); pa <- ncol(Xa)
set.seed(SEED); cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
lam_a <- na * cvf$lambda.1se
X1a <- cbind(1, Xa); Dpen <- diag(c(0, rep(1, pa))); ba <- cf_ridge(X1a, ya, lam_a, Dpen)
cat(sprintf("lambda = %.4f\n", lam_a))

newton <- function(X1, y, lambda, Dpen, beta, tol = 1e-9, maxit = 50) {   # warm-started ridge Newton
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1, w * X1) + lambda * Dpen, crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn)
    beta <- bn
  }
  beta
}
sim_y <- function(r) { set.seed(SEED + 1000 * 1 + r); rbinom(na, 1, cf_expit(drop(X1a %*% ba))) }   # tag "ads" = 1, as in E35

one_rep <- function(r) {
  y <- sim_y(r); bh <- newton(X1a, y, lam_a, Dpen, ba)
  eta <- drop(X1a %*% bh); pi <- cf_expit(eta); w <- pi * (1 - pi)
  M <- crossprod(X1a, w * X1a) + lam_a * Dpen
  q <- rowSums((X1a %*% solve(M)) * X1a); H <- w * q
  xi <- eta - q / (1 - H) * (y - pi)
  set.seed(SEED + 11 * r); fold <- sample(rep_len(1:10, na)); xk <- numeric(na)
  for (k in 1:10) { tr <- fold != k; bk <- newton(X1a[tr, ], y[tr], lam_a * 9 / 10, Dpen, bh); xk[!tr] <- drop(X1a[!tr, ] %*% bk) }
  list(r = r, y = y, eta = eta, eta_true = drop(X1a %*% ba), xi = xi, xk = xk, q = q, H = H, beta = bh,
       p1 = cav_from_index(xi, y)$p_chi2, pk = cav_from_index(xk, y)$p_chi2)
}

t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c("Xa", "X1a", "ya", "ba", "lam_a", "Dpen", "na", "SEED", "newton", "sim_y", "one_rep", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
reps <- parLapply(cl, seq_len(R), one_rep)
cat(sprintf("part A done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# exact leave-one-out on 240 observations of replicate 1
r1 <- reps[[1]]; set.seed(7); idx <- sample(na, 240)
clusterExport(cl, "r1")
xl <- unlist(parLapply(cl, idx, function(i) { b <- newton(X1a[-i, ], r1$y[-i], lam_a, Dpen, r1$beta); sum(X1a[i, ] * b) }))
stopCluster(cl)
cat(sprintf("part B done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(list(reps = reps, idx = idx, xi_exact = xl, lam_a = lam_a, nfeat = rowSums(Xa)), file.path(OUT, "E35b_ads_diagnosis.rds"))

err <- r1$xi[idx] - xl
cat(sprintf("one-step vs exact LOO (240 obs, rep 1): cor = %.4f, mean abs err = %.3f, max abs err = %.3f, sd(exact xi) = %.3f\n",
            cor(r1$xi[idx], xl), mean(abs(err)), max(abs(err)), sd(xl)))
cat("per-replicate P values, one-step and ten-fold:\n")
print(data.frame(rep = sapply(reps, `[[`, "r"), one_step = round(sapply(reps, `[[`, "p1"), 4), tenfold = round(sapply(reps, `[[`, "pk"), 4)))
cat(sprintf("rejections at 5%%: one-step %d/%d, ten-fold %d/%d\n", sum(sapply(reps, `[[`, "p1") < .05), R, sum(sapply(reps, `[[`, "pk") < .05), R))
