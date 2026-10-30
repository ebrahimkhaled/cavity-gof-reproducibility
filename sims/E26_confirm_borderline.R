# E26 -- confirmation runs for the two borderline cells, on FRESH seeds (never reuse a seed to confirm):
#   (1) E25 bent_n20000: q = 2 size .100 (drift, as predicted) and q = 3 size .073 at R = 300 -- is q = 3 valid?  R = 1000.
#   (2) E22 enet5_sparse10_p2n: decile size .070 at R = 500 -- real or chance?  R = 1000.
# Usage: Rscript E26_confirm_borderline.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20261106
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
p_b <- 50; Sb <- 0.9^abs(outer(1:p_b, 1:p_b, "-")); bb <- rep_len(c(.35, -.3, .25, -.2, .15), p_b); bb <- bb * 2.31 / sqrt(drop(t(bb) %*% Sb %*% bb)); Cb <- chol(Sb)
p_e <- 600; Se <- 0.7^abs(outer(1:p_e, 1:p_e, "-")); be <- c(rep_len(c(.35, -.3, .25, -.2, .15), 10), rep(0, p_e - 10)); be <- be * 2.31 / sqrt(drop(t(be) %*% Se %*% be)); Ce <- chol(Se)
one_bent <- function(r) {
  set.seed(SEED + r); n <- 20000; X <- matrix(rnorm(n * p_b), n) %*% Cb; y <- rbinom(n, 1, cf_expit(drop(X %*% bb)))
  c2 <- cav_gof_fast(X, y, n * 2, G = 10); c3 <- cav_lognormal(c2$xi, y, G = 10)
  data.frame(rep = r, cell = "bent_n20000", q2 = c2$p_a, q3 = c3$p_a, S2 = c2$S, S3 = c3$S, tr3 = c3$tr, tr2 = c2$tr_a)
}
one_enet <- function(r) {
  set.seed(SEED + 1e5 + r); n <- 300; X <- matrix(rnorm(n * p_e), n) %*% Ce; y <- rbinom(n, 1, cf_expit(drop(X %*% be)))
  c3 <- cav_gof_glmnet(X, y, 0.04, alpha = .5, G = 10, basis_degree = 3)
  data.frame(rep = r, cell = "enet5_sparse10_p2n", q2 = c3$p_a, q3 = c3$p_E, S2 = c3$S, S3 = NA, tr3 = NA, tr2 = c3$tr_a)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(glmnet)})))
f_out <- file.path(OUT, "E26_confirm_perrep.csv"); if (file.exists(f_out)) file.remove(f_out)
for (job in list(list(f = one_enet, name = "enet"), list(f = one_bent, name = "bent"))) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, jobfun) tryCatch(jobfun(r), error = function(e) NULL), jobfun = job$f))
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("done %s n=%d (%.1f min)\n", job$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out); rej <- function(v) mean(v < .05, na.rm = TRUE); se <- function(v) sqrt(.05 * .95 / sum(!is.na(v)))
for (cc in unique(res$cell)) { s <- res[res$cell == cc, ]
  cat(sprintf("%-20s R = %d | decile/q2 size %.3f | %s size %.3f | MC SE %.4f | mean S2 %.2f (tr %.2f) | mean S3 %.2f (tr %.2f)\n", cc, nrow(s),
              rej(s$q2), if (cc == "bent_n20000") "q3" else "smooth", rej(s$q3), se(s$q2), mean(s$S2), mean(s$tr2), mean(s$S3), mean(s$tr3))) }
