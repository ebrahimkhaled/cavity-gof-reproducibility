# E36b -- how stable is the cross-validated glaucoma penalty, and the corrected cavity P-value at it?
# 10-fold cv.glmnet (ridge) under 50 fold seeds (20261200 + k, L'Ecuyer reset before each), two set-ups:
#   inside: raw predictors, standardize = TRUE (glmnet standardises inside each fold)
#   global: predictors standardised on all n first, standardize = FALSE (the set-up that gave 295.3 in the companion paper)
# For every seed the corrected test (E36 version B: exact refits, standardisation inside each leave-one-out fit, q = 2,
# G = 10, chi2_8) is computed at that seed's lambda_1se and lambda_min.
# Output: results/E36b_glaucoma_cv_seeds.csv (one row per seed x set-up), summary printed to results/E36b_glaucoma_cv_seeds.log.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
suppressPackageStartupMessages({ library(TH.data); library(glmnet) }); data("GlaucomaM", package = "TH.data")
Xraw <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]); yg <- as.numeric(GlaucomaM$Class == "glaucoma"); n <- nrow(Xraw)
ridge_ws <- function(X1, y, lambda, beta = rep(0, ncol(X1)), tol = 1e-10, maxit = 100) {
  Dp <- c(0, rep(1, ncol(X1) - 1))
  for (it in seq_len(maxit)) { eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1, w * X1) + diag(lambda * Dp), crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn); beta <- bn }
  beta
}
std_with <- function(X, mu, s) sweep(sweep(X, 2, mu), 2, s, "/")
loo_inside <- function(Xr, y, lambda) {
  b_full <- ridge_ws(cbind(1, scale(Xr)), y, lambda); xi <- numeric(nrow(Xr))
  for (i in seq_len(nrow(Xr))) { Z1 <- cbind(1, std_with(Xr, colMeans(Xr[-i, ]), apply(Xr[-i, ], 2, sd)))
    xi[i] <- sum(Z1[i, ] * ridge_ws(Z1[-i, ], y[-i], lambda, beta = b_full)) }
  xi
}
f_out <- file.path(OUT, "E36b_glaucoma_cv_seeds.csv")
have <- if (file.exists(f_out)) read.csv(f_out) else NULL
for (k in 1:50) for (su in c("inside", "global")) {
  if (!is.null(have) && any(have$k == k & have$setup == su)) next
  RNGkind("L'Ecuyer-CMRG"); set.seed(20261200 + k)
  cv <- if (su == "inside") cv.glmnet(Xraw, yg, family = "binomial", alpha = 0, standardize = TRUE, nfolds = 10) else
    cv.glmnet(scale(Xraw), yg, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
  l1 <- n * cv$lambda.1se; l0 <- n * cv$lambda.min
  p1 <- cav_from_index(loo_inside(Xraw, yg, l1), yg, 10)$p_chi2; p0 <- cav_from_index(loo_inside(Xraw, yg, l0), yg, 10)$p_chi2
  o <- data.frame(k = k, setup = su, lambda_1se = l1, lambda_min = l0, p_1se = p1, p_min = p0)
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
}
r <- read.csv(f_out)
sink(file.path(OUT, "E36b_glaucoma_cv_seeds.log"), split = TRUE)
for (su in c("inside", "global")) { s <- r[r$setup == su, ]
  cat(sprintf("%-6s (%d seeds): lambda_1se median %.1f [range %.1f, %.1f]; lambda_min median %.1f [%.1f, %.1f]\n", su, nrow(s),
              median(s$lambda_1se), min(s$lambda_1se), max(s$lambda_1se), median(s$lambda_min), min(s$lambda_min), max(s$lambda_min)))
  cat(sprintf("        corrected cavity P at lambda_1se: median %.3f [range %.3f, %.3f], share < .05 = %.2f; at lambda_min: median %.3f [%.3f, %.3f], share < .05 = %.2f\n",
              median(s$p_1se), min(s$p_1se), max(s$p_1se), mean(s$p_1se < .05), median(s$p_min), min(s$p_min), max(s$p_min), mean(s$p_min < .05))) }
sink()
