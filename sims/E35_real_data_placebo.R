# E35 -- real data: (1) the UCI Internet Advertisements data (Kushmerick 1999) analysed with a ridge fit at the penalty
# chosen by ten-fold cross-validation (one-standard-error rule), and (2) placebo checks of size on the REAL predictors of
# both data sets: y is simulated from the fitted model on the observed X, the model is refitted at the same penalty, and
# the cavity and textbook tests are applied. The binary Internet-Ads predictors are used as they are (no standardisation,
# so no fit uses observation i through a scaling); the glaucoma predictors are standardised once, as in the companion paper.
# Usage: Rscript E35_real_data_placebo.R [R_placebo] [workers]
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20261133
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

# ---- Internet Advertisements -----------------------------------------------------------------------------------------
d <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
d <- d[!is.na(d$V4), ]                                    # 15 rows with the 'local' indicator missing
Xa <- as.matrix(d[, 4:(ncol(d) - 1)]); ya <- as.numeric(d[, ncol(d)] == "ad.")
Xa <- Xa[, colSums(Xa) >= 5]                              # drop features present in fewer than 5 images
na <- nrow(Xa); pa <- ncol(Xa)
set.seed(SEED); cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
lam_a <- na * cvf$lambda.1se                              # glmnet's penalty on the scale of the paper's objective
ca <- cav_gof_fast(Xa, ya, lam_a, G = 10)
X1a <- cbind(1, Xa); ba <- cf_ridge(X1a, ya, lam_a, diag(c(0, rep(1, pa)))); eta_a <- drop(X1a %*% ba)
la <- cav_lognormal(ca$xi, ya, G = 10)
sink(file.path(OUT, "E35_internet_ads_analysis.txt"))
cat(sprintf("Internet Advertisements: n = %d, p = %d (p/n = %.2f), %d adverts (%.1f%%)\n", na, pa, pa / na, sum(ya), 100 * mean(ya)))
cat(sprintf("ridge penalty (1-se rule): glmnet lambda = %.5f, paper scale lambda = %.2f\n", cvf$lambda.1se, lam_a))
cat(sprintf("cavity test, two-parameter: S = %.2f, P = %.4f (chi2_8), P = %.4f (exact projection); slope b = %.2f, mean leverage %.3f\n",
            ca$S, ca$p_chi2, ca$p_a, ca$b, ca$Hbar))
cat(sprintf("cavity test, three-parameter: S = %.2f, P = %.4f (chi2_8), tau = %.3f\n", la$S, pchisq(la$S, 8, lower.tail = FALSE), la$tau))
cat(sprintf("textbook Hosmer-Lemeshow on the ridge fit: P = %.3g\n", naive_hl(eta_a, ya, 10)))
cat("standardised cavity residuals by decile:", sprintf("%.2f", ca$e), "\n")
sink()

# ---- glaucoma ---------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
Xg <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); yg <- as.numeric(GlaucomaM$Class == "glaucoma")
lam_g <- 295.3; bg <- cf_ridge(cbind(1, Xg), yg, lam_g, diag(c(0, rep(1, ncol(Xg)))))

# ---- placebo: simulate from the fitted model on the real X, refit at the same penalty ---------------------------------
placebo <- function(r, Xm, beta, lambda, tag) {     # Xm, not X: parLapply has its own argument X
  set.seed(SEED + 1000 * match(tag, c("ads", "glaucoma")) + r)
  X1 <- cbind(1, Xm); y <- rbinom(nrow(Xm), 1, cf_expit(drop(X1 %*% beta)))
  bh <- cf_ridge(X1, y, lambda, diag(c(0, rep(1, ncol(Xm)))))   # one fit, shared by both tests
  cv <- cav_gof_fast(Xm, y, lambda, G = 10, beta = bh)
  data.frame(rep = r, data = tag, cavity = cv$p_chi2, cavity_exact = cv$p_a, textbook = naive_hl(drop(X1 %*% bh), y, 10), b = cv$b)
}
f_out <- file.path(OUT, "E35_placebo_perrep.csv")
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, c("placebo", "SEED", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
# an Internet-Ads replicate costs about seven CPU-minutes, so results are appended one batch of W replicates at a
# time and a rerun skips the replicates already on disk
for (job in list(list("glaucoma", Xg, bg, lam_g), list("ads", Xa, ba, lam_a))) {
  have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$data == job[[1]]] } else integer(0)
  todo <- setdiff(seq_len(R), have)
  for (batch in split(todo, ceiling(seq_along(todo) / W))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r, Xm, beta, lambda, tag) tryCatch(placebo(r, Xm, beta, lambda, tag),
                                  error = function(e) NULL), Xm = job[[2]], beta = job[[3]], lambda = job[[4]], tag = job[[1]]))
    write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
    cat(sprintf("placebo %-9s reps %d-%d  (%.1f min)\n", job[[1]], min(batch), max(batch),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
stopCluster(cl)
res <- read.csv(f_out); rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- do.call(rbind, lapply(split(res, res$data), function(s) data.frame(data = s$data[1], R = nrow(s), cavity = rej(s$cavity),
  cavity_exact = rej(s$cavity_exact), textbook = rej(s$textbook), slope = mean(s$b))))
write.csv(agg, file.path(OUT, "E35_placebo_summary.csv"), row.names = FALSE)
cat("\n"); cat(readLines(file.path(OUT, "E35_internet_ads_analysis.txt")), sep = "\n")
cat("\n=== E35 | placebo size on the real predictors ===\n"); print(agg, digits = 3, row.names = FALSE)
