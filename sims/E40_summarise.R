# E40_summarise -- cavity (lasso / ridge, q = 2 / q = 3, chi2_8) against GRPtest in the sparse designs of E40.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); RES <- file.path(ROOT, "results")
f <- function(file) { d <- read.csv(file.path(RES, file)); d[!duplicated(d[, c("cell", "rep")]), ] }
rr <- function(v) { v <- v[!is.na(v)]; m <- mean(v < .05); sprintf("%.3f (%.3f)", m, sqrt(m * (1 - m) / length(v))) }
summ <- function(d) do.call(rbind, lapply(split(d, d$cell), function(s) data.frame(cell = s$cell[1], R = nrow(s), R_grp = sum(!is.na(s$grp)),
  lasso_q2 = rr(s$lasso_q2), lasso_q3 = rr(s$lasso_q3), ridge_q2 = rr(s$ridge_q2), ridge_q3 = rr(s$ridge_q3), grp = rr(s$grp),
  grp_median_p = round(median(s$grp, na.rm = TRUE), 3), textbook = rr(s$textbook), nnz = round(mean(s$nnz), 1))))
sink(file.path(RES, "E40_summary.log"), split = TRUE)
cat("=== pilot (null only; seeds 20261150 + 1e5 + r and + 3e5 + r) ===\n"); options(width = 250); print(summ(f("E40_grp_sparse_pilot_perrep.csv")), row.names = FALSE)
m <- f("E40_grp_sparse_perrep.csv"); s <- summ(m); write.csv(s, file.path(RES, "E40_grp_sparse_summary.csv"), row.names = FALSE)
cat("\n=== main (seeds 20261150 + r) ===\n"); print(s, row.names = FALSE)
for (cc in unique(m$cell)) { z <- m[m$cell == cc & !is.na(m$grp) & !is.na(m$lasso_q2), ]; dd <- (z$lasso_q2 < .05) - (z$grp < .05)
  cat(sprintf("%-11s paired cavity(lasso, q=2) - GRP: %+.3f (%.3f); both reject %d, cavity only %d, GRP only %d\n", cc, mean(dd), sd(dd) / sqrt(nrow(z)),
              sum(z$lasso_q2 < .05 & z$grp < .05), sum(z$lasso_q2 < .05 & z$grp >= .05), sum(z$lasso_q2 >= .05 & z$grp < .05))) }
sink()
