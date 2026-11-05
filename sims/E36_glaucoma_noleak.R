# E36 -- glaucoma re-analysis WITHOUT leakage (referee item m9 / S13).
# The paper's fifth usage rule: the fit that produces xi_i must not use observation i. E21 standardised the 62
# predictors on all 196 eyes before fitting, so every leave-one-out fit used x_i through the column means and SDs.
# Here the leave-one-out index is computed by EXACT refits (n = 196 ridge fits), each one standardising the predictors
# on the n - 1 eyes it is fitted to and applying that same centring and scaling to the held-out eye:
#   xi_i = (1, (x_i - mu_(-i)) / s_(-i))' beta_hat_(-i).
# The penalty is the same lambda in every refit (the one-step index of the paper approximates exactly this).
# Raw predictors with one scalar penalty were NOT used as the primary analysis: the raw SDs range from 0.003 to 0.78,
# so a scalar penalty on the raw scale is a different model (it is reported as a secondary analysis, version C, which
# is what cavity.gof() computes when it is given the raw matrix).
# Versions, all at G = 10 unless stated:
#   A  one-step index, predictors standardised on all n          (= E21, the paper's numbers)
#   A' exact refits,   predictors standardised on all n          (isolates the one-step approximation)
#   B  exact refits,   predictors standardised inside each fit   (the corrected analysis: PRIMARY)
#   C  one-step index, raw predictors, penalty chosen by CV on the raw scale (cavity.gof() route, no standardisation)
# Penalty grid: 295.3 x (0.5, 1, 2) and the cross-validated lambda_min and lambda_1se (cv.glmnet, standardize = TRUE,
# which standardises inside each fold; 10 folds, seed 20261136).
# Outputs: results/E36_glaucoma_noleak.csv (one row per version x lambda x G), results/E36_glaucoma_deciles.csv,
#          results/E36_glaucoma_noleak.log (sink), figures/fig5_glaucoma_v2.pdf.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
suppressPackageStartupMessages({ library(TH.data); library(glmnet) }); data("GlaucomaM", package = "TH.data")
Xraw <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]); yg <- as.numeric(GlaucomaM$Class == "glaucoma")
n <- nrow(Xraw); p <- ncol(Xraw); SEED <- 20261136

# ridge by Newton with a warm start (same objective as cf_ridge: log-likelihood minus lambda/2 |slopes|^2)
ridge_ws <- function(X1, y, lambda, beta = rep(0, ncol(X1)), tol = 1e-10, maxit = 100) {
  Dp <- c(0, rep(1, ncol(X1) - 1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1, w * X1) + diag(lambda * Dp), crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn)
    beta <- bn
  }
  beta
}
std_with <- function(X, mu, s) sweep(sweep(X, 2, mu), 2, s, "/")
# exact leave-one-out index; inside = TRUE standardises on the n - 1 training rows of each refit
loo_index <- function(Xr, y, lambda, inside = TRUE) {
  n <- nrow(Xr); mu0 <- colMeans(Xr); s0 <- apply(Xr, 2, sd)
  b_full <- ridge_ws(cbind(1, std_with(Xr, mu0, s0)), y, lambda)
  xi <- numeric(n)
  for (i in seq_len(n)) {
    if (inside) { mu <- colMeans(Xr[-i, ]); s <- apply(Xr[-i, ], 2, sd) } else { mu <- mu0; s <- s0 }
    Z1 <- cbind(1, std_with(Xr, mu, s))
    b <- ridge_ws(Z1[-i, ], y[-i], lambda, beta = b_full)
    xi[i] <- sum(Z1[i, ] * b)
  }
  xi
}
# test on a given index: q = 2 (chi2_{G-2} and exact projection), q = 3 with chi2_{G-2}, smooth degree 3
test_index <- function(xi, y, G = 10) {
  c2 <- cav_from_index(xi, y, G, basis_degree = 3)
  c3 <- tryCatch(cav_lognormal(xi, y, G), error = function(e) NULL)
  list(c2 = c2, row = data.frame(S2 = c2$S, p2_chi2 = c2$p_chi2, p2_exact = c2$p_a, p_smooth3 = c2$p_E, b = c2$b, a = c2$alpha,
       S3 = if (is.null(c3)) NA else c3$S, p3_chi2G2 = if (is.null(c3)) NA else pchisq(c3$S, G - 2, lower.tail = FALSE),
       p3_exact = if (is.null(c3)) NA else c3$p_a, tau3 = if (is.null(c3)) NA else c3$tau))
}

sink(file.path(OUT, "E36_glaucoma_noleak.log"), split = TRUE)
cat(sprintf("GlaucomaM: n = %d, p = %d, events = %d; raw predictor SDs %.4f to %.3f\n", n, p, sum(yg), min(apply(Xraw, 2, sd)), max(apply(Xraw, 2, sd))))

# ---- cross-validated penalties (standardisation inside each fold) --------------------------------------------------
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
cvs <- cv.glmnet(Xraw, yg, family = "binomial", alpha = 0, standardize = TRUE, nfolds = 10)
lam_cv <- c(cv_min = n * cvs$lambda.min, cv_1se = n * cvs$lambda.1se)
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
cvr <- cv.glmnet(Xraw, yg, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
lam_raw <- c(raw_cv_min = n * cvr$lambda.min, raw_cv_1se = n * cvr$lambda.1se)
cat(sprintf("CV (standardised inside folds): lambda_min = %.1f, lambda_1se = %.1f (theory scale)\n", lam_cv[1], lam_cv[2]))
cat(sprintf("CV (raw predictors):            lambda_min = %.2f, lambda_1se = %.2f (theory scale)\n", lam_raw[1], lam_raw[2]))
LAMS <- c(half = 295.3 / 2, base = 295.3, double = 2 * 295.3, lam_cv)

rows <- list(); Xs_all <- scale(Xraw)
for (ln in names(LAMS)) {
  lam <- LAMS[[ln]]
  xiA <- cav_gof_fast(Xs_all, yg, lam, G = 10)$xi
  xiA2 <- loo_index(Xraw, yg, lam, inside = FALSE)
  xiB <- loo_index(Xraw, yg, lam, inside = TRUE)
  bfull <- cf_ridge(cbind(1, Xs_all), yg, lam, diag(c(0, rep(1, p)))); tb <- naive_hl(drop(cbind(1, Xs_all) %*% bfull), yg, 10)
  for (G in if (ln == "base") c(5, 10, 20) else 10) {
    for (v in list(list("A_onestep_globalstd", xiA), list("A2_exact_globalstd", xiA2), list("B_exact_insidestd", xiB))) {
      tt <- test_index(v[[2]], yg, G)$row
      rows[[length(rows) + 1]] <- data.frame(version = v[[1]], lambda_name = ln, lambda = lam, G = G, tt, textbook_p = tb)
    }
  }
  cat(sprintf("lambda %-6s %6.1f: max |xi_A - xi_A'| = %.4f, max |xi_A' - xi_B| = %.4f, sd(xi_B) = %.3f\n", ln, lam,
              max(abs(xiA - xiA2)), max(abs(xiA2 - xiB)), sd(xiB)))
  if (ln == "base") { xiB_base <- xiB; xiA_base <- xiA }
}
# version C: raw predictors, cavity.gof() route (one-step index, no standardisation), CV penalties on the raw scale
for (ln in names(lam_raw)) {
  lam <- lam_raw[[ln]]; cc <- cav_gof_fast(Xraw, yg, lam, G = 10)
  braw <- cf_ridge(cbind(1, Xraw), yg, lam, diag(c(0, rep(1, p)))); tb <- naive_hl(drop(cbind(1, Xraw) %*% braw), yg, 10)
  rows[[length(rows) + 1]] <- data.frame(version = "C_onestep_raw", lambda_name = ln, lambda = lam, G = 10, test_index(cc$xi, yg, 10)$row, textbook_p = tb)
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E36_glaucoma_noleak.csv"), row.names = FALSE)
options(width = 220); print(res, digits = 3, row.names = FALSE)

# ---- decile table at lambda = 295.3, versions A and B -----------------------------------------------------------------
dec <- do.call(rbind, lapply(list(list("A_onestep_globalstd", xiA_base), list("B_exact_insidestd", xiB_base)), function(v) {
  c2 <- cav_from_index(v[[2]], yg, 10)
  m <- cf_expit(c2$alpha + c2$b * v[[2]]); g <- pmin(ceiling(rank(v[[2]], ties.method = "first") / (n / 10)), 10)
  idx <- split(seq_len(n), g)
  data.frame(version = v[[1]], g = seq_along(idx), n = sapply(idx, length), obs = sapply(idx, function(I) sum(yg[I])),
             exp = sapply(idx, function(I) sum(m[I])), e = c2$e)
}))
write.csv(dec, file.path(OUT, "E36_glaucoma_deciles.csv"), row.names = FALSE)
cat("\nObserved vs expected by cavity decile, lambda = 295.3, G = 10:\n"); print(dec, digits = 3, row.names = FALSE)
# how many eyes change decile between A and B
gA <- pmin(ceiling(rank(xiA_base, ties.method = "first") / (n / 10)), 10); gB <- pmin(ceiling(rank(xiB_base, ties.method = "first") / (n / 10)), 10)
cat(sprintf("eyes that change decile between A and B: %d of %d; correlation of the two indices %.5f\n", sum(gA != gB), n, cor(xiA_base, xiB_base)))
sink()

# ---- Figure 5, version 2 (same style as make_figures.R; corrected index B) ---------------------------------------------
COL <- c(cav = "#0072B2", ins = "#D55E00", txt = "#CC79A7"); GRID <- "grey88"; INK2 <- "grey35"; FIG <- file.path(ROOT, "figures")
c2 <- cav_from_index(xiB_base, yg, 10); mB <- cf_expit(c2$alpha + c2$b * xiB_base)
Xg1 <- cbind(1, Xs_all); bg <- cf_ridge(Xg1, yg, 295.3, diag(c(0, rep(1, p)))); pih <- cf_expit(drop(Xg1 %*% bg))
idx <- split(seq_along(yg), gB)
ob <- sapply(idx, function(I) mean(yg[I])); ex <- sapply(idx, function(I) mean(mB[I])); nn <- sapply(idx, length); fi <- sapply(idx, function(I) mean(pih[I]))
se <- sqrt(ex * (1 - ex) / nn); worst <- which.max(abs(c2$e))
pdf(file.path(FIG, "fig5_glaucoma_v2.pdf"), width = 4.6, height = 4.0, pointsize = 10)
par(cex = 1, mar = c(4, 4.3, 1, 1), mgp = c(2.4, 0.7, 0), las = 1)
plot(c(0, 1), c(0, 1), type = "n", xlab = "Expected proportion", ylab = "Observed proportion", axes = FALSE)
abline(h = seq(0, 1, .2), v = seq(0, 1, .2), col = GRID); abline(0, 1, lty = 2, col = "grey50")
arrows(ex, pmax(ex - 1.96 * se, 0), ex, pmin(ex + 1.96 * se, 1), angle = 90, code = 3, length = 0.03, col = adjustcolor(COL["cav"], 0.45))
points(fi, ob, pch = 24, bg = "white", col = COL["ins"], cex = 1.05, lwd = 1.4)
points(ex, ob, pch = 21, bg = COL["cav"], col = "white", cex = 1.3)
text(ex[worst], ob[worst], sprintf("decile %d", worst), pos = 2, cex = 0.8, col = "grey20")
axis(1, col = INK2, col.axis = INK2); axis(2, col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
legend("topleft", bty = "n", cex = 0.8, pch = c(21, 24), pt.bg = c(COL["cav"], "white"), col = c("white", COL["ins"]), pt.cex = c(1.3, 1),
       legend = c("against the recalibrated cavity curve", "against the shrunk ridge probabilities"))
invisible(dev.off())
cat("fig5_glaucoma_v2.pdf written\n")
