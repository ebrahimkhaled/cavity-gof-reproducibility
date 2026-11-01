# E35f -- on the 12 Internet-Ads placebo data sets of E35b: (1) split each grouped residual into a noise part
# (y - true pi) and a misfit part (true pi - recalibrated m), to see whether the failure is a variance or a mean
# problem; (2) the kernel version of the cluster correction, every pair weighted by its leave-one-out cosine.
# Reads results/E35b_ads_diagnosis.rds. Usage: Rscript E35f_ads_misfit_decomposition.R
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R")); source(file.path(ROOT, "R", "cluster_correction.R"))
o <- readRDS(file.path(OUT, "E35b_ads_diagnosis.rds"))

E <- D <- N <- NULL
for (k in seq_along(o$reps)) {
  r <- o$reps[[k]]; set.seed(k)
  rc <- cav_recal_fit(r$xi, r$y); m <- cf_expit(rc[1] + rc[2] * r$xi); v <- m * (1 - m); pt <- cf_expit(r$eta_true)
  g <- pmin(ceiling(rank(r$xi, ties.method = "random") / (length(r$y) / 10)), 10); sV <- sqrt(tapply(v, g, sum))
  E <- rbind(E, tapply(r$y - m, g, sum) / sV); D <- rbind(D, tapply(pt - m, g, sum) / sV); N <- rbind(N, tapply(r$y - pt, g, sum) / sV)
}

d <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
d <- d[!is.na(d$V4), ]
Xa <- as.matrix(d[, 4:(ncol(d) - 1)]); Xa <- Xa[, colSums(Xa) >= 5]
K <- t(sapply(seq_along(o$reps), function(k) {
  r <- o$reps[[k]]; C <- cav_cosine(Xa, r$beta, o$lam_a); set.seed(k)
  cav_from_index_kernel(r$xi, r$y, C, c(0.5, 0.7, 0.9, 0.99))
}))

sink(file.path(OUT, "E35f_ads_misfit_decomposition.txt"))
cat("=== E35f | misfit part (true pi - m) of the grouped residuals, by decile, 12 placebo data sets ===\nmean:\n")
print(round(colMeans(D), 2)); cat("sd:\n"); print(round(apply(D, 2, sd), 2))
cat("noise part (y - true pi): sd by decile\n"); print(round(apply(N, 2, sd), 2))
cat(sprintf("mean S = %.1f; mean sum misfit^2 = %.1f; mean sum noise^2 = %.1f; 2 mean sum cross = %.1f\n",
            mean(rowSums(E^2)), mean(rowSums(D^2)), mean(rowSums(N^2)), 2 * mean(rowSums(D * N))))
cat("\n=== kernel correction: P values (omega = closed form; k<c0> = pairs weighted by max(0, (cos - c0)/(1 - c0))) ===\n")
print(round(K, 4)); cat("rejections at 5%:", colSums(K < .05), "of", nrow(K), "\n")
sink()
cat(readLines(file.path(OUT, "E35f_ads_misfit_decomposition.txt")), sep = "\n")
