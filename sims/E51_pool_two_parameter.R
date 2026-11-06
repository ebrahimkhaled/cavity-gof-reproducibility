HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
# Two-parameter size in the strongly correlated design, pooling E25 (300 replications) and E51 (1000 on new seeds): Section 6.3 and Supplementary Section S6.
a <- read.csv(file.path(ROOT, "results", "E25_q2_vs_q3_perrep.csv")); b <- read.csv(file.path(ROOT, "results", "E51_confirm_perrep.csv"))
for (n in c(2000, 5000, 20000)) { cn <- paste0("bent_n", n); p <- c(a$q2_chi2[a$cell == cn], b$full_T0[b$cell == cn]); s <- qchisq(p, 8, lower.tail = FALSE)
  cat(n, length(p), "size", round(mean(p < .05), 4), "se", round(sqrt(mean(p < .05) * (1 - mean(p < .05)) / length(p)), 4), "mean S", round(mean(s), 3), "se", round(sd(s) / sqrt(length(s)), 3), "\n") }
