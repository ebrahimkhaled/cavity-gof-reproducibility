# E38_summarise -- tables and the pooled P-value QQ plot from results/E38_q3_size_all_perrep.csv.
# (1) fingerprint: the recomputed q = 2 P-values against those stored by E18 / E19 / E22 (same data => same P-values);
# (2) size of q = 2 (chi2_{G-2}) and q = 3 (chi2_{G-2}) at alpha = 0.01, 0.05, 0.10 with Monte Carlo SEs, every design;
# (3) booktabs table results/E38_size_table.tex; (4) figures/figS_qq.pdf (paper figure style).
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
d <- read.csv(file.path(RES, "E38_q3_size_all_perrep.csv")); d <- d[!duplicated(d[, c("cell", "rep")]), ]
# q = 3 statistics that cav_lognormal() could not return (singular information matrix) were recomputed by E42
d$q3_repaired <- FALSE
if (file.exists(file.path(RES, "E42_q3_repair.csv"))) { rp <- read.csv(file.path(RES, "E42_q3_repair.csv")); rp <- rp[rp$source == "E38", ]
  for (k in seq_len(nrow(rp))) { i <- which(d$cell == rp$cell[k] & d$rep == rp$rep[k] & is.na(d$q3_chi2G2))
    d$q3_chi2G2[i] <- rp$q3_chi2G2[k]; d$S3[i] <- rp$S3[k]; d$tau3[i] <- rp$tau3[k]; d$q3_repaired[i] <- TRUE } }
old <- rbind(read.csv(file.path(RES, "E18_cavity_size_perrep.csv"))[, c("rep", "cell", "p_chi2")],
             read.csv(file.path(RES, "E19_cavity_stress_perrep.csv"))[, c("rep", "cell", "p_chi2")],
             read.csv(file.path(RES, "E22_cavity_lasso_perrep.csv"))[, c("rep", "cell", "p_chi2")])
m <- merge(d, old, by = c("cell", "rep"))
fp <- aggregate(cbind(maxdiff = abs(q2_chi2 - p_chi2)) ~ cell, m, max)
fp$n_matched <- aggregate(rep ~ cell, m, length)$rep
fp$rej_flip <- aggregate(cbind(f = (q2_chi2 < .05) != (p_chi2 < .05)) ~ cell, m, sum)$f
sink(file.path(RES, "E38_summary.log"), split = TRUE)
cat("Fingerprint (recomputed q=2 chi2 P vs stored): max abs difference over all matched replicates =", signif(max(fp$maxdiff), 3),
    "; decisions at .05 that differ:", sum(fp$rej_flip), "of", sum(fp$n_matched), "\n")
print(fp[order(-fp$maxdiff), ][1:5, ], row.names = FALSE)
ORDER <- c("A_p5", "B_p40", "B_p100", "B_p160", "H3_k40_lam50", "H5_k25_ar03", "H8_n800_k25", "k10_t3", "prev15_k25", "prev08_k25", "prev08_A",
           "n1600_k25", "n3200_k25", "strong4_k25", "weak1_k25", "lam10_k25", "lam5000_k25", "k50", "binary_k25", "t3_k25", "G5_k25", "G20_k25",
           "prev05_strong3_k25", "A_n5000", "lasso_dense_k25_l02", "lasso_dense_k25_l05", "enet5_dense_k25_l05", "lasso_sparse5_k25",
           "lasso_sparse5_k50", "lasso_sparse5_p2n", "enet5_sparse10_p2n", "lasso_sparse5_prev10")
LAB <- c(A_p5 = "fixed $p$", B_p40 = "$p/n$ 0.10", B_p100 = "$p/n$ 0.25", B_p160 = "$p/n$ 0.40", H3_k40_lam50 = "$p/n$ 0.40, light penalty",
         H5_k25_ar03 = "weak correlation", H8_n800_k25 = "$n = 800$", k10_t3 = "$t_3$ predictors", prev15_k25 = "prevalence 0.15",
         prev08_k25 = "prevalence 0.08", prev08_A = "fixed $p$, prevalence 0.08", n1600_k25 = "$n = 1600$", n3200_k25 = "$n = 3200$",
         strong4_k25 = "strong signal", weak1_k25 = "weak signal", lam10_k25 = "light penalty", lam5000_k25 = "heavy penalty", k50 = "$p/n$ 0.50",
         binary_k25 = "binary predictors", t3_k25 = "$t_3$ predictors", G5_k25 = "$G = 5$", G20_k25 = "$G = 20$", prev05_strong3_k25 = "prevalence 0.05",
         A_n5000 = "fixed $p$, $n = 5000$", lasso_dense_k25_l02 = "lasso, dense", lasso_dense_k25_l05 = "lasso, dense, heavier",
         enet5_dense_k25_l05 = "elastic net, dense", lasso_sparse5_k25 = "lasso, sparse", lasso_sparse5_k50 = "lasso, sparse, $p/n$ 0.5",
         lasso_sparse5_p2n = "lasso, sparse, $p = 2n$", enet5_sparse10_p2n = "elastic net, $p = 2n$", lasso_sparse5_prev10 = "lasso, prevalence 0.10")
PANEL <- setNames(c(rep("A", 11), rep("B", 13), rep("C", 8)), ORDER)
rate <- function(v, a) { v <- v[!is.na(v)]; c(r = mean(v < a), se = sqrt(mean(v < a) * (1 - mean(v < a)) / length(v))) }
tab <- do.call(rbind, lapply(ORDER, function(cc) { s <- d[d$cell == cc, ]
  o <- data.frame(panel = PANEL[[cc]], cell = cc, R = nrow(s), q3_failed = sum(is.na(s$q3_chi2G2)), q3_repaired = sum(s$q3_repaired))
  for (a in c(.01, .05, .10)) for (t in c("q2_chi2", "q3_chi2G2")) { z <- rate(s[[t]], a); o[[sprintf("%s_%03d", t, round(100 * a))]] <- z[["r"]]; o[[sprintf("%s_%03d_se", t, round(100 * a))]] <- z[["se"]] }
  o$q3_exact_05 <- mean(s$q3_exact < .05, na.rm = TRUE); o$meanS2 <- mean(s$S2); o$meanS3 <- mean(s$S3, na.rm = TRUE); o$tr3 <- mean(s$tr3, na.rm = TRUE)
  o$diff05 <- mean((s$q3_chi2G2 < .05) - (s$q2_chi2 < .05), na.rm = TRUE)
  o$diff05_se <- sd((s$q3_chi2G2 < .05) - (s$q2_chi2 < .05), na.rm = TRUE) / sqrt(sum(!is.na(s$q3_chi2G2)))
  o }))
# 500-replicate rates for the three extended designs too (same replicates as Figure 3)
ext <- do.call(rbind, lapply(c("B_p100", "B_p160", "lasso_sparse5_p2n"), function(cc) { s <- d[d$cell == cc & d$rep <= 500, ]
  data.frame(cell = cc, R = nrow(s), q2_05 = mean(s$q2_chi2 < .05), q3_05 = mean(s$q3_chi2G2 < .05, na.rm = TRUE)) }))
write.csv(tab, file.path(RES, "E38_size_summary.csv"), row.names = FALSE)
options(width = 250)
show <- tab[, c("panel", "cell", "R", "q3_failed", "q2_chi2_001", "q3_chi2G2_001", "q2_chi2_005", "q3_chi2G2_005", "q2_chi2_010", "q3_chi2G2_010", "q3_exact_05", "diff05", "diff05_se", "meanS2", "meanS3", "tr3")]
print(show, digits = 3, row.names = FALSE)
cat("\nThe three extended designs at their first 500 replicates (the Figure 3 data):\n"); print(ext, digits = 3, row.names = FALSE)
for (pn in c("A", "B", "C", "all")) { s <- if (pn == "all") tab else tab[tab$panel == pn, ]
  cat(sprintf("panel %-3s: q2 size at .05 range %.3f-%.3f, q3 %.3f-%.3f; at .01 q2 %.3f-%.3f, q3 %.3f-%.3f; at .10 q2 %.3f-%.3f, q3 %.3f-%.3f; mean paired diff q3-q2 at .05 %+.4f\n",
              pn, min(s$q2_chi2_005), max(s$q2_chi2_005), min(s$q3_chi2G2_005), max(s$q3_chi2G2_005), min(s$q2_chi2_001), max(s$q2_chi2_001),
              min(s$q3_chi2G2_001), max(s$q3_chi2G2_001), min(s$q2_chi2_010), max(s$q2_chi2_010), min(s$q3_chi2G2_010), max(s$q3_chi2G2_010), mean(s$diff05))) }
out2 <- function(k, a) { z <- tab[k, ]; sprintf("%.3f", z[[a]]) }
flag <- function(r, se, a, R) if (abs(r - a) > 2 * sqrt(a * (1 - a) / R)) "$^*$" else ""
lines <- c("% generated by sims/E38_summarise.R -- size of the cavity test, q = 2 and q = 3, both with the chi2_{G-2} reference",
           "\\begin{tabular}{llrrrrrrr}", "\\toprule",
           " & & & \\multicolumn{2}{c}{$\\alpha = 0.01$} & \\multicolumn{2}{c}{$\\alpha = 0.05$} & \\multicolumn{2}{c}{$\\alpha = 0.10$} \\\\",
           "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}", "Panel & Design & $R$ & $q=2$ & $q=3$ & $q=2$ & $q=3$ & $q=2$ & $q=3$ \\\\", "\\midrule")
for (k in seq_len(nrow(tab))) { z <- tab[k, ]
  cells <- unlist(lapply(c("001", "005", "010"), function(a) { al <- as.numeric(a) / 100; sapply(c("q2_chi2", "q3_chi2G2"), function(t)
    paste0(sprintf("%.3f", z[[sprintf("%s_%s", t, a)]]), flag(z[[sprintf("%s_%s", t, a)]], 0, al, z$R))) }))
  lines <- c(lines, sprintf("%s & %s & %d & %s \\\\", if (k == 1 || tab$panel[k] != tab$panel[k - 1]) z$panel else "", LAB[[z$cell]], z$R, paste(cells, collapse = " & "))) }
lines <- c(lines, "\\bottomrule", "\\end{tabular}",
           "% $^*$: more than two binomial standard errors sqrt(alpha(1-alpha)/R) from alpha.")
writeLines(lines, file.path(RES, "E38_size_table.tex"))
sink()

# ---- QQ plot of P-values, pooled over the eleven ridge designs of panel A, first 500 replicates each (5500) ----------
COL <- c(cav = "#0072B2", ins = "#D55E00", txt = "#CC79A7"); GRID <- "grey88"; INK2 <- "grey35"
panel_lab <- function(l) mtext(l, side = 3, adj = -0.18, line = 0.2, font = 2, cex = 1.05)
pa <- d[d$cell %in% ORDER[1:11] & d$rep <= 500, ]
p2 <- sort(pa$q2_chi2); p3 <- sort(pa$q3_chi2G2[!is.na(pa$q3_chi2G2)]); N <- length(p2); N3 <- length(p3)
u <- (seq_len(N) - 0.5) / N; u3 <- (seq_len(N3) - 0.5) / N3
k <- seq_len(N); lo <- qbeta(0.025, k, N - k + 1); hi <- qbeta(0.975, k, N - k + 1)
pdf(file.path(FIG, "figS_qq.pdf"), width = 7.2, height = 3.5, pointsize = 10)
layout(matrix(1:2, 1)); par(cex = 1, mar = c(4, 4.3, 1.6, 1), mgp = c(2.4, 0.7, 0), las = 1)
plot(c(0, 1), c(0, 1), type = "n", xlab = "Uniform quantile", ylab = "Ordered P-value", axes = FALSE)
abline(h = seq(0, 1, .2), v = seq(0, 1, .2), col = GRID)
polygon(c(u, rev(u)), c(lo, rev(hi)), col = "grey85", border = NA); abline(0, 1, lty = 2, col = "grey50")
th <- unique(round(seq(1, N, length = 400))); th3 <- unique(round(seq(1, N3, length = 400)))
points(u3[th3], p3[th3], pch = 21, bg = "white", col = COL["cav"], cex = 0.55, lwd = 1.1)
points(u[th], p2[th], pch = 16, col = COL["cav"], cex = 0.45)
axis(1, col = INK2, col.axis = INK2); axis(2, col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
legend("topleft", bty = "n", cex = 0.85, pch = c(16, 21), pt.bg = c(NA, "white"), col = COL[c("cav", "cav")],
       legend = c("two-parameter", "three-parameter"))
panel_lab("A")
xl <- c(5e-4, 0.25)
plot(xl, xl, type = "n", log = "xy", xlab = "Uniform quantile", ylab = "Ordered P-value", axes = FALSE)
abline(v = c(.001, .01, .05, .1), h = c(.001, .01, .05, .1), col = GRID)
sel <- u <= 0.25 & u >= 5e-4; polygon(c(u[sel], rev(u[sel])), c(lo[sel], rev(hi[sel])), col = "grey85", border = NA); abline(0, 1, lty = 2, col = "grey50")
s3 <- u3 <= 0.25 & u3 >= 5e-4
points(u3[s3], pmax(p3[s3], 1e-4), pch = 21, bg = "white", col = COL["cav"], cex = 0.55, lwd = 1.1)
points(u[sel], pmax(p2[sel], 1e-4), pch = 16, col = COL["cav"], cex = 0.45)
tk <- c(.001, .01, .05, .1, .25); axis(1, at = tk, labels = c("0.001", "0.01", "0.05", "0.1", "0.25"), col = INK2, col.axis = INK2)
axis(2, at = tk, labels = c("0.001", "0.01", "0.05", "0.1", "0.25"), col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
panel_lab("B")
invisible(dev.off())
cat(sprintf("figS_qq.pdf: pooled %d replicates (q=2), %d (q=3) from the 11 designs of Figure 3 panel A\n", N, N3))
