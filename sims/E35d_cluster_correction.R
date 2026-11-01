# E35d -- validation of the cluster correction (R/cluster_correction.R) on real predictors.
# (1) Internet-Ads placebo, same seeds as E35: size of the cavity test with Omega and with the cluster-corrected
#     reference. (2) Internet-Ads power: responses from the fitted index plus a quadratic departure. (3) glaucoma
#     placebo with the correction. (4) the corrected test on the real Internet-Ads responses.
# Clusters: connected components of leave-one-out cosine > 0.999 (repeated rows; a cut of 0.9 chains most observations
# together when p is small, e.g. p = 5), recomputed from each fit.
# Usage: Rscript E35d_cluster_correction.R [R_size] [R_power] [R_glaucoma] [workers]
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
Rs <- if (length(args) >= 1) as.integer(args[1]) else 120; Rp <- if (length(args) >= 2) as.integer(args[2]) else 60
Rg <- if (length(args) >= 3) as.integer(args[3]) else 500; W <- if (length(args) >= 4) as.integer(args[4]) else 12
SEED <- 20261133; CUT <- 0.999; QUAD <- 0.25          # departure: eta + QUAD * sd(eta) * (z^2 - 1), z = standardised eta
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R")); source(file.path(ROOT, "R", "cluster_correction.R"))

d <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
d <- d[!is.na(d$V4), ]
Xa <- as.matrix(d[, 4:(ncol(d) - 1)]); ya <- as.numeric(d[, ncol(d)] == "ad.")
Xa <- Xa[, colSums(Xa) >= 5]; na <- nrow(Xa); pa <- ncol(Xa)
set.seed(SEED); cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
lam_a <- na * cvf$lambda.1se
ba <- cf_ridge(cbind(1, Xa), ya, lam_a, diag(c(0, rep(1, pa))))
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
Xg <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); yg <- as.numeric(GlaucomaM$Class == "glaucoma")
lam_g <- 295.3; bg <- cf_ridge(cbind(1, Xg), yg, lam_g, diag(c(0, rep(1, ncol(Xg)))))

newton <- function(X1, y, lambda, Dpen, beta, tol = 1e-9, maxit = 50) {
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1, w * X1) + lambda * Dpen, crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn)
    beta <- bn
  }
  beta
}
# one fit, the one-step cavity index, near-twin clusters from the same fit, and both references
test_once <- function(Xm, y, lambda, beta0, seed) {
  X1 <- cbind(1, Xm); Dpen <- diag(c(0, rep(1, ncol(Xm)))); bh <- newton(X1, y, lambda, Dpen, beta0)
  eta <- drop(X1 %*% bh); pi <- cf_expit(eta); w <- pi * (1 - pi)
  q <- rowSums((X1 %*% solve(crossprod(X1, w * X1) + lambda * Dpen)) * X1)
  xi <- eta - q / (1 - w * q) * (y - pi)
  cl <- cav_clusters(Xm, bh, lambda, CUT); set.seed(seed)
  a <- cav_from_index_cluster(xi, y, cl)
  data.frame(omega = a$p_omega, cluster = a$p_cluster, twin = cav_twin_share(cl), tr = a$trace_ratio, textbook = naive_hl(eta, y, 10))
}
job <- function(r, arm) {
  if (arm == "glaucoma") { Xm <- Xg; b <- bg; lam <- lam_g; set.seed(SEED + 2000 + r) } else { Xm <- Xa; b <- ba; lam <- lam_a; set.seed(SEED + 1000 + r) }
  eta <- drop(cbind(1, Xm) %*% b)
  if (arm == "ads_power") { z <- (eta - mean(eta)) / sd(eta); eta <- eta + QUAD * sd(eta) * (z^2 - 1) }
  y <- rbinom(nrow(Xm), 1, cf_expit(eta))
  cbind(data.frame(rep = r, arm = arm), test_once(Xm, y, lam, b, SEED + 7 * r))
}

real <- test_once(Xa, ya, lam_a, ba, SEED)
cat(sprintf("real Internet-Ads data: twin share %.2f, trace ratio %.2f, P(Omega) = %.4f, P(cluster) = %.4f\n", real$twin, real$tr, real$omega, real$cluster))
real_g <- test_once(Xg, yg, lam_g, bg, SEED)
cat(sprintf("real glaucoma data: twin share %.2f, trace ratio %.2f, P(Omega) = %.4f, P(cluster) = %.4f\n", real_g$twin, real_g$tr, real_g$omega, real_g$cluster))

f_out <- file.path(OUT, "E35d_cluster_perrep.csv")
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c("Xa", "ba", "lam_a", "Xg", "bg", "lam_g", "SEED", "CUT", "QUAD", "newton", "test_once", "job", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
for (arm in list(list("glaucoma", Rg), list("ads_size", Rs), list("ads_power", Rp))) {
  have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$arm == arm[[1]]] } else integer(0)
  todo <- setdiff(seq_len(arm[[2]]), have)
  B <- if (arm[[1]] == "glaucoma") 10 * W else W
  for (batch in split(todo, ceiling(seq_along(todo) / B))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r, a) tryCatch(job(r, a), error = function(e) NULL), a = arm[[1]]))
    write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
    cat(sprintf("%-9s reps %d-%d  (%.1f min)\n", arm[[1]], min(batch), max(batch), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
stopCluster(cl)
res <- read.csv(f_out); rj <- function(v) mean(v < .05)
agg <- do.call(rbind, lapply(split(res, res$arm), function(s) data.frame(arm = s$arm[1], R = nrow(s), omega = rj(s$omega),
  cluster = rj(s$cluster), se_cluster = sqrt(rj(s$cluster) * (1 - rj(s$cluster)) / nrow(s)), textbook = rj(s$textbook),
  twin = mean(s$twin), trace_ratio = mean(s$tr))))
write.csv(agg, file.path(OUT, "E35d_cluster_summary.csv"), row.names = FALSE)
cat("\n=== E35d | cluster correction on real predictors ===\n"); print(agg, digits = 3, row.names = FALSE)
