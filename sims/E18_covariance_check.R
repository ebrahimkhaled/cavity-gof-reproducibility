# Matrix-level check of the cavity reference: the empirical covariance of the standardised group residuals
# (E18 per-replicate vectors) should equal Omega_a = I - (rank-2 recalibration projection): eigenvalues 1 (x8) and 0 (x2).
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
res <- read.csv(file.path(ROOT, "results", "E18_cavity_size_perrep.csv"))
ec <- grep("^e[0-9]+$", names(res))
cat("cell               eigenvalues of Cov(e) (sorted)                                  | mean(e) norm2 | max |Cov - I| on the 8-dim complement\n")
for (cc in unique(res$cell)) {
  E <- as.matrix(res[res$cell == cc, ec]); C <- cov(E); ev <- sort(eigen(C, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  cat(sprintf("%-16s %s | %.3f\n", cc, paste(sprintf("%.2f", ev), collapse = " "), sum(colMeans(E)^2)))
}
