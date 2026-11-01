# E35c -- the sample-split cavity test on the Internet-Ads placebo data sets of E35 (same seeds, so the same y).
# E35 found the full-sample test rejects 57 per cent of correct models there; E35b traced this to repeated
# predictor rows (51 per cent of images share their row), which make the leave-one-out indices depend on each
# other's responses. The split index never depends on a response in the half it is tested on, so this run asks
# whether removing that dependence restores the level. Results are appended one batch at a time.
# Usage: Rscript E35c_ads_split_placebo.R [R] [workers]
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 120; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261133
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

d <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
d <- d[!is.na(d$V4), ]
Xa <- as.matrix(d[, 4:(ncol(d) - 1)]); ya <- as.numeric(d[, ncol(d)] == "ad.")
Xa <- Xa[, colSums(Xa) >= 5]; na <- nrow(Xa); pa <- ncol(Xa)
set.seed(SEED); cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
lam_a <- na * cvf$lambda.1se
ba <- cf_ridge(cbind(1, Xa), ya, lam_a, diag(c(0, rep(1, pa))))
real <- cav_split_gof(Xa, ya, lam_a, seed = SEED)
cat(sprintf("real data: split P = %.4f, cross-fitted P = %.4f\n", real$p_split, real$p_cross))

placebo_split <- function(r) {
  set.seed(SEED + 1000 * 1 + r); y <- rbinom(na, 1, cf_expit(drop(cbind(1, Xa) %*% ba)))   # the E35 "ads" seed
  s <- cav_split_gof(Xa, y, lam_a, seed = SEED + r)
  data.frame(rep = r, split = s$p_split, cross = s$p_cross)
}
f_out <- file.path(OUT, "E35c_ads_split_perrep.csv")
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c("Xa", "ba", "lam_a", "na", "SEED", "placebo_split", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
have <- if (file.exists(f_out)) read.csv(f_out)$rep else integer(0)
todo <- setdiff(seq_len(R), have)
for (batch in split(todo, ceiling(seq_along(todo) / W))) {
  o <- do.call(rbind, parLapply(cl, batch, function(r) tryCatch(placebo_split(r), error = function(e) NULL)))
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("reps %d-%d  (%.1f min)\n", min(batch), max(batch), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out)
cat(sprintf("\n=== E35c | Internet-Ads placebo, %d replications ===\nsplit test size %.3f (SE %.3f); cross-fitted size %.3f (SE %.3f)\n",
            nrow(res), mean(res$split < .05), sqrt(mean(res$split < .05) * (1 - mean(res$split < .05)) / nrow(res)),
            mean(res$cross < .05), sqrt(mean(res$cross < .05) * (1 - mean(res$cross < .05)) / nrow(res))))
