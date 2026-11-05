# E37 -- glaucoma placebo (E35, paper: 0.054) re-run under the CORRECTED standardisation (E36 version B).
# y is simulated from the fitted ridge model (lambda = 295.3, predictors standardised on all n, as E35's truth) on the
# REAL predictors; the model is refitted and tested. Per replicate, on the same y:
#   old  : one-step index, predictors standardised on all n (E35's test), q = 2, chi2_8 and exact projection
#   new  : exact leave-one-out refits, standardisation inside each fit (the usage rule), q = 2 chi2_8 / exact, q = 3 chi2_8
#   textbook Hosmer-Lemeshow on the full fit
# Seeds: RNGkind("L'Ecuyer-CMRG") reset, then set.seed(20261137 + r), for every replicate. Appended in batches; a rerun
# skips replicates already on disk.  Usage: Rscript E37_glaucoma_placebo_noleak.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000; W <- if (length(args) >= 2) as.integer(args[2]) else 20
SEED <- 20261137
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
Xraw <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]); yg <- as.numeric(GlaucomaM$Class == "glaucoma")
Xs <- scale(Xraw); LAM <- 295.3; p <- ncol(Xraw)
b_true <- cf_ridge(cbind(1, Xs), yg, LAM, diag(c(0, rep(1, p)))); pr_true <- cf_expit(drop(cbind(1, Xs) %*% b_true))

ridge_ws <- function(X1, y, lambda, beta = rep(0, ncol(X1)), tol = 1e-10, maxit = 100) {
  Dp <- c(0, rep(1, ncol(X1) - 1))
  for (it in seq_len(maxit)) { eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1, w * X1) + diag(lambda * Dp), crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn); beta <- bn }
  beta
}
std_with <- function(X, mu, s) sweep(sweep(X, 2, mu), 2, s, "/")
loo_inside <- function(Xr, y, lambda, b_full) {
  xi <- numeric(nrow(Xr))
  for (i in seq_len(nrow(Xr))) { Z1 <- cbind(1, std_with(Xr, colMeans(Xr[-i, ]), apply(Xr[-i, ], 2, sd)))
    xi[i] <- sum(Z1[i, ] * ridge_ws(Z1[-i, ], y[-i], lambda, beta = b_full)) }
  xi
}
one <- function(r) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(SEED + r); y <- rbinom(length(pr_true), 1, pr_true)
  X1 <- cbind(1, Xs); bh <- cf_ridge(X1, y, LAM, diag(c(0, rep(1, p))))
  old <- cav_gof_fast(Xs, y, LAM, G = 10, beta = bh)
  xi <- loo_inside(Xraw, y, LAM, bh)
  n2 <- cav_from_index(xi, y, 10); n3 <- tryCatch(cav_lognormal(xi, y, 10), error = function(e) NULL)
  data.frame(rep = r, old_chi2 = old$p_chi2, old_exact = old$p_a, new_q2_chi2 = n2$p_chi2, new_q2_exact = n2$p_a,
             new_q3_chi2G2 = if (is.null(n3)) NA else pchisq(n3$S, 8, lower.tail = FALSE), tau3 = if (is.null(n3)) NA else n3$tau,
             textbook = naive_hl(drop(X1 %*% bh), y, 10), b = n2$b, ysum = sum(y))
}
f_out <- file.path(OUT, "E37_glaucoma_placebo_noleak_perrep.csv")
have <- if (file.exists(f_out)) read.csv(f_out)$rep else integer(0); todo <- setdiff(seq_len(R), have)
t0 <- Sys.time()
if (length(todo)) {
  cl <- makeCluster(W); clusterExport(cl, c("one", "SEED", "LAM", "p", "Xs", "Xraw", "pr_true", "ridge_ws", "std_with", "loo_inside",
                                            ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
  for (batch in split(todo, ceiling(seq_along(todo) / (5 * W)))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r) tryCatch(one(r), error = function(e) NULL)))
    write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
    cat(sprintf("reps %d-%d (%.1f min)\n", min(batch), max(batch), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  stopCluster(cl)
}
res <- read.csv(f_out); Rn <- nrow(res)
rej <- function(v, a) { m <- mean(v < a, na.rm = TRUE); sprintf("%.3f (%.3f)", m, sqrt(m * (1 - m) / sum(!is.na(v)))) }
tab <- do.call(rbind, lapply(c("old_chi2", "old_exact", "new_q2_chi2", "new_q2_exact", "new_q3_chi2G2", "textbook"), function(v)
  data.frame(test = v, a01 = rej(res[[v]], .01), a05 = rej(res[[v]], .05), a10 = rej(res[[v]], .10))))
write.csv(cbind(R = Rn, tab), file.path(OUT, "E37_glaucoma_placebo_noleak_summary.csv"), row.names = FALSE)
cat(sprintf("\n=== E37 | glaucoma placebo, corrected standardisation | R = %d ===\n", Rn)); print(tab, row.names = FALSE)
cat(sprintf("paired old vs new (q2 chi2) rejections at .05: both %d, old only %d, new only %d\n", sum(res$old_chi2 < .05 & res$new_q2_chi2 < .05),
            sum(res$old_chi2 < .05 & res$new_q2_chi2 >= .05), sum(res$old_chi2 >= .05 & res$new_q2_chi2 < .05)))
