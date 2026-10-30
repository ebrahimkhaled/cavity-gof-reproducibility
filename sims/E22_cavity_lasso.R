# E22 -- the cavity-recalibrated test after LASSO and ELASTIC-NET fits (glmnet), under the null, including p > n.
# The exact cavity law holds for any estimator built from D_(-i); the leave-one-out approximation is on the active set.
# Compared with the textbook Hosmer-Lemeshow test on the in-sample fitted probabilities (what practice does).
# Usage: Rscript E22_cavity_lasso.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261101
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda_g, alpha = 1, ar = 0.7, norm = 2.31, sparse = 0, prev = .5) {
  Sig <- ar^abs(outer(1:p, 1:p, "-"))
  b0 <- if (sparse > 0) c(rep_len(c(.35, -.3, .25, -.2, .15), sparse), rep(0, p - sparse)) else rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0))
  list(name = name, n = n, p = p, lambda_g = lambda_g, alpha = alpha, b0 = b0, Sc = chol(Sig), prev = prev) }
CELLS <- list(
  mk("lasso_dense_k25_l02", 400, 100, 0.02), mk("lasso_dense_k25_l05", 400, 100, 0.05),
  mk("enet5_dense_k25_l05", 400, 100, 0.05, alpha = .5),
  mk("lasso_sparse5_k25", 400, 100, 0.03, sparse = 5), mk("lasso_sparse5_k50", 400, 200, 0.03, sparse = 5),
  mk("lasso_sparse5_p2n", 300, 600, 0.04, sparse = 5), mk("enet5_sparse10_p2n", 300, 600, 0.04, alpha = .5, sparse = 10),
  mk("lasso_sparse5_prev10", 400, 100, 0.02, sparse = 5, prev = .10))
draw <- function(d, r) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta, d$prev)
  list(X = X, y = rbinom(d$n, 1, cf_expit(a + eta)))
}
one <- function(r, d) {
  z <- draw(d, r)
  c3 <- cav_gof_glmnet(z$X, z$y, d$lambda_g, alpha = d$alpha, G = 10, basis_degree = 3)
  data.frame(rep = r, cell = d$name, p_chi2 = c3$p_chi2, p_a = c3$p_a, p_E3 = c3$p_E, naive = c3$naive,
             S = c3$S, tr_a = c3$tr_a, b = c3$b, nnz = c3$nnz, Hbar = c3$Hbar)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(glmnet)})))
f_out <- file.path(OUT, "E22_cavity_lasso_perrep.csv"); if (file.exists(f_out)) file.remove(f_out)
for (d in CELLS) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
  if (is.null(o)) { cat("ALL FAILED", d$name, "\n"); next }
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("done %-22s n=%d  (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out); rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(p_chi2, p_a, p_E3, naive) ~ cell, res, rej, na.action = na.pass)
agg$meanS <- aggregate(S ~ cell, res, mean)$S; agg$tr_a <- aggregate(tr_a ~ cell, res, mean)$tr_a; agg$b <- aggregate(b ~ cell, res, mean)$b
agg$nnz <- aggregate(nnz ~ cell, res, mean)$nnz; agg$Hbar <- aggregate(Hbar ~ cell, res, mean)$Hbar; agg$R <- aggregate(rep ~ cell, res, length)$rep
agg <- agg[match(sapply(CELLS, `[[`, "name"), agg$cell), ]
write.csv(agg, file.path(OUT, "E22_cavity_lasso_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E22 | cavity test after lasso / elastic net (null) | %.1f min ===\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
