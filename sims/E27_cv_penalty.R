# E27 -- the cavity test when the penalty is CHOSEN BY CROSS-VALIDATION on the same data (what practice does),
# under the null. Ridge: cv.glmnet(alpha = 0, 10 folds, standardize = FALSE), lambda_min and lambda_1se, converted to the
# theory scale lambda = n * lambda_glmnet; lasso: lambda_1se. Compared with the textbook test at the same penalty.
# Usage: Rscript E27_cv_penalty.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 300; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20261107
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, ar = 0.7, norm = 2.31, prev = .5, indep = FALSE, sparse = 0) {
  if (indep) { Sig <- diag(p); b0 <- c(.5, -.4, .3, -.3, .2) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-"))
    b0 <- if (sparse > 0) c(rep_len(c(.35, -.3, .25, -.2, .15), sparse), rep(0, p - sparse)) else rep_len(c(.35, -.3, .25, -.2, .15), p)
    b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  list(name = name, n = n, p = p, b0 = b0, Sc = chol(Sig), prev = prev) }
CELLS <- list(mk("A_p5", 500, 5, indep = TRUE), mk("B_p40", 400, 40), mk("B_p100", 400, 100), mk("H5_k25_ar03", 400, 100, ar = .3),
              mk("prev08_k25", 400, 100, prev = .08), mk("sparse5_p2n", 300, 600, sparse = 5))
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta0 <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta0, d$prev); y <- rbinom(d$n, 1, cf_expit(a + eta0))
  out <- data.frame(rep = r, cell = d$name)
  for (s in c("lambda.min", "lambda.1se")) for (k in c("_cav", "_text", "_lam")) out[[paste0("ridge_", s, k)]] <- NA_real_
  if (d$p < d$n) {                                                     # ridge: our own exact fit at the CV-chosen penalty
    cv <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10, standardize = FALSE)
    for (s in c("lambda.min", "lambda.1se")) {
      lam <- d$n * cv[[s]]; c2 <- cav_gof_fast(X, y, lam, G = 10)
      X1 <- cbind(1, X); bh <- cf_ridge(X1, y, lam, diag(c(0, rep(1, d$p))))
      out[[paste0("ridge_", s, "_cav")]] <- c2$p_a; out[[paste0("ridge_", s, "_text")]] <- naive_hl(drop(X1 %*% bh), y, 10)
      out[[paste0("ridge_", s, "_lam")]] <- lam
    }
  }
  cvl <- glmnet::cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10, standardize = FALSE)
  cl1 <- cav_gof_glmnet(X, y, cvl$lambda.1se, alpha = 1, G = 10)
  out$lasso_1se_cav <- cl1$p_a; out$lasso_1se_text <- cl1$naive; out$lasso_1se_nnz <- cl1$nnz
  out
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(glmnet)})))
res <- do.call(rbind, lapply(CELLS, function(d) { o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
  cat(sprintf("done %-12s n=%d (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl)
write.csv(res, file.path(OUT, "E27_cv_penalty_perrep.csv"), row.names = FALSE)
rej <- function(v) if (all(is.na(v))) NA else mean(v < .05, na.rm = TRUE)
pc <- grep("_cav$|_text$", names(res), value = TRUE)
agg <- do.call(rbind, lapply(split(res, res$cell), function(s) data.frame(cell = s$cell[1], t(sapply(pc, function(v) rej(s[[v]]))),
  lam_min = mean(s$ridge_lambda.min_lam), lam_1se = mean(s$ridge_lambda.1se_lam), nnz = mean(s$lasso_1se_nnz), R = nrow(s))))
write.csv(agg, file.path(OUT, "E27_cv_penalty_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E27 | cavity test with the penalty chosen by cross-validation (null) | %.1f min ===\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
