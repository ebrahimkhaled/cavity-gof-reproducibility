# Is the excess of the 3-parameter test a weak-identification effect of estimating tau? If so, referring S3 to chi2_{G-2}
# (i.e. not subtracting the tau direction) should be conservative-to-nominal. Pools E25, E26 and E29 per-replicate files.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); R <- paste0(file.path(ROOT, "results"), "/")
a <- read.csv(paste0(R, "E25_q2_vs_q3_perrep.csv")); b <- read.csv(paste0(R, "E26_confirm_perrep.csv")); c <- read.csv(paste0(R, "E29_q3_excess_n20000_perrep.csv"))
rows <- list(
  data.frame(src = "E25", cell = a$cell, S3 = a$S3, tr3 = a$tr3, p3 = a$q3_exact),
  data.frame(src = "E26", cell = b$cell, S3 = b$S3, tr3 = b$tr3, p3 = b$q3),
  data.frame(src = "E29", cell = "bent_n20000", S3 = c$S3, tr3 = c$tr3, p3 = c$p3))
d <- do.call(rbind, rows); d <- d[!is.na(d$S3), ]
d$p_chi8 <- pchisq(d$S3, 8, lower.tail = FALSE)
for (cc in unique(d$cell)) { s <- d[d$cell == cc, ]
  cat(sprintf("%-14s R %4d | exact q3 ref: size %.3f | chi2_8 ref: size %.3f | mean S3 %.2f (SE %.2f) trace %.2f\n",
              cc, nrow(s), mean(s$p3 < .05), mean(s$p_chi8 < .05), mean(s$S3), sd(s$S3) / sqrt(nrow(s)), mean(s$tr3))) }
