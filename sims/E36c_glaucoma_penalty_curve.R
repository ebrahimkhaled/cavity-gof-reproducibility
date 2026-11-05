# E36c -- glaucoma: the corrected cavity P-value (E36 version B: exact refits, standardisation inside each leave-one-out
# fit) as a function of the ridge penalty, on a log grid of 25 values from 5 to 1200 (theory scale), G = 10.
# q = 2 (chi2_8), q = 3 (chi2_8), smooth degree 3, the decile-5 residual and the textbook test on the full fit.
# Output: results/E36c_glaucoma_penalty_curve.csv. Deterministic (no random numbers).
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
Xraw <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]); yg <- as.numeric(GlaucomaM$Class == "glaucoma"); p <- ncol(Xraw)
ridge_ws <- function(X1, y, lambda, beta = rep(0, ncol(X1)), tol = 1e-10, maxit = 100) {
  Dp <- c(0, rep(1, ncol(X1) - 1))
  for (it in seq_len(maxit)) { eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1, w * X1) + diag(lambda * Dp), crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn); beta <- bn }
  beta
}
std_with <- function(X, mu, s) sweep(sweep(X, 2, mu), 2, s, "/")
f_out <- file.path(OUT, "E36c_glaucoma_penalty_curve.csv")
have <- if (file.exists(f_out)) read.csv(f_out)$lambda else numeric(0)
for (lam in signif(exp(seq(log(5), log(1200), length = 25)), 4)) {
  if (any(abs(have - lam) < 1e-8)) next
  b_full <- ridge_ws(cbind(1, scale(Xraw)), yg, lam); xi <- numeric(nrow(Xraw))
  for (i in seq_len(nrow(Xraw))) { Z1 <- cbind(1, std_with(Xraw, colMeans(Xraw[-i, ]), apply(Xraw[-i, ], 2, sd)))
    xi[i] <- sum(Z1[i, ] * ridge_ws(Z1[-i, ], yg[-i], lam, beta = b_full)) }
  c2 <- cav_from_index(xi, yg, 10, basis_degree = 3); c3 <- tryCatch(cav_lognormal(xi, yg, 10), error = function(e) NULL)
  o <- data.frame(lambda = lam, p2 = c2$p_chi2, p3 = if (is.null(c3)) NA else pchisq(c3$S, 8, lower.tail = FALSE),
                  p_smooth3 = c2$p_E, e_min = min(c2$e), g_min = which.min(c2$e), e5 = c2$e[5], b = c2$b,
                  textbook = naive_hl(drop(cbind(1, scale(Xraw)) %*% b_full), yg, 10))
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
}
print(read.csv(f_out), digits = 3, row.names = FALSE)
