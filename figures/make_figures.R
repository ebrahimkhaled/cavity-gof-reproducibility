# Figures for the cavity paper. Base R, vector PDF, no titles inside images, par(cex = 1) after layout().
# Palette validated with the dataviz validator (light surface): cavity #0072B2, in-sample #D55E00, textbook #CC79A7;
# identity is ALSO carried by shape (circle / triangle / square), never colour alone. Every number read from results/.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "figures")
COL <- c(cav = "#0072B2", ins = "#D55E00", txt = "#CC79A7"); GRID <- "grey88"; INK2 <- "grey35"
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
panel_lab <- function(l) mtext(l, side = 3, adj = -0.18, line = 0.2, font = 2, cex = 1.05)

# ---- Figure 1: the mechanism ---------------------------------------------------------------------------------------
# A: one data set (design B, p/n = 0.25): self-influence eta_hat - xi against xi, by outcome, with decile boundaries of eta_hat
set.seed(20261112); n <- 400; p <- 100; S <- 0.7^abs(outer(1:p, 1:p, "-")); X <- matrix(rnorm(n * p), n) %*% chol(S)
b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(drop(t(b0) %*% S %*% b0)); y <- rbinom(n, 1, cf_expit(drop(X %*% b0)))
X1 <- cbind(1, X); D <- diag(c(0, rep(1, p))); bh <- cf_ridge(X1, y, 416, D); eta <- drop(X1 %*% bh); pi <- cf_expit(eta); w <- pi * (1 - pi)
q <- rowSums((X1 %*% solve(crossprod(X1, w * X1) + 416 * D)) * X1); xi <- eta - q / (1 - w * q) * (y - pi)
# B: E23 exact leave-one-out indices (10,000 observations): observed proportion vs Proposition-2 curve vs the shrunk curve pi(xi)
e23 <- read.csv(file.path(RES, "E23_verify_prop1_perobs_seed20261103.csv")); e23 <- e23[e23$cell == "B_p100", ]
br <- quantile(e23$xi_exact, 0:20 / 20); g <- cut(e23$xi_exact, br, include.lowest = TRUE)
bx <- tapply(e23$xi_exact, g, mean); by <- tapply(e23$y, g, mean); bn <- as.vector(table(g))
xs <- seq(min(e23$xi_exact), max(e23$xi_exact), length = 200)
lam <- sapply(xs, function(x) { k <- which.min(abs(e23$xi_exact - x)); NA })
ord <- order(e23$xi_exact); lam_s <- stats::lowess(e23$xi_exact[ord], e23$Lam_exact[ord], f = 0.05)
pdf(file.path(FIG, "fig1_mechanism.pdf"), width = 7.2, height = 3.3, pointsize = 10)
layout(matrix(1:2, 1)); par(cex = 1, mar = c(4, 5.2, 1.6, 1), mgp = c(2.4, 0.7, 0), las = 1)
plot(xi, eta - xi, type = "n", xlab = "Leave-one-out index", ylab = "", axes = FALSE); title(ylab = "Fitted minus leave-one-out index", line = 3.6)
abline(v = quantile(eta, 1:9 / 10), col = GRID, lwd = 1); abline(h = 0, col = "grey70")
points(xi[y == 0], (eta - xi)[y == 0], pch = 17, cex = 0.55, col = adjustcolor(COL["ins"], 0.8))
points(xi[y == 1], (eta - xi)[y == 1], pch = 16, cex = 0.55, col = adjustcolor(COL["cav"], 0.8))
axis(1, col = INK2, col.axis = INK2); axis(2, col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
legend("topright", bty = "n", cex = 0.85, pch = c(16, 17), col = COL[c("cav", "ins")], legend = c("y = 1", "y = 0"))
panel_lab("A")
plot(range(xs), c(0, 1), type = "n", xlab = "Leave-one-out index", ylab = "Probability of an event", axes = FALSE)
grid(nx = NA, ny = NULL, col = GRID, lty = 1)
lines(xs, cf_expit(xs), col = COL["ins"], lwd = 2, lty = 2)
lines(lam_s, col = COL["cav"], lwd = 2)
points(bx, by, pch = 21, bg = "white", col = "black", cex = 0.9)
axis(1, col = INK2, col.axis = INK2); axis(2, col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
# short labels so the legend stays right of x = 0, where both curves are above 0.5
legend("bottomright", bty = "n", cex = 0.85, seg.len = 1.6, lty = c(1, 2, NA), pch = c(NA, NA, 21), lwd = c(2, 2, NA), col = c(COL["cav"], COL["ins"], "black"),
       legend = c("Proposition 2", "ridge fit", "observed"))
panel_lab("B")
invisible(dev.off())

# ---- Figure 2: self-selection, predicted vs simulated (E24) ---------------------------------------------------------
s <- read.csv(file.path(RES, "E24_selection_formula_pergroup.csv"))
cells <- c("A_p5", "B_p40", "B_p100", "B_p160", "H3_k40_lam50", "H5_k25_ar03", "H8_n800_k25", "prev15_k25", "prev08_k25")
obs <- sapply(cells, function(cc) sum(tapply(s$e_ins[s$cell == cc], s$g[s$cell == cc], mean)^2))
pre <- sapply(cells, function(cc) sum(tapply(s$pred[s$cell == cc], s$g[s$cell == cc], mean)^2))
pdf(file.path(FIG, "fig2_selection.pdf"), width = 7.2, height = 3.3, pointsize = 10)
layout(matrix(1:2, 1), widths = c(1, 1.25)); par(cex = 1, mar = c(4, 4.3, 1.6, 1), mgp = c(2.4, 0.7, 0), las = 1)
plot(obs, pre, log = "xy", xlim = c(0.5, 80), ylim = c(0.5, 80), type = "n", xlab = "Simulated non-centrality", ylab = "Predicted non-centrality", axes = FALSE)
abline(h = c(1, 3, 10, 30), v = c(1, 3, 10, 30), col = GRID); abline(0, 1, lty = 2, col = "grey50")
points(obs, pre, pch = 21, bg = COL["cav"], col = "white", cex = 1.3)
text(obs, pre, labels = seq_along(obs), pos = 4, cex = 0.7, col = INK2)
axis(1, at = c(1, 3, 10, 30), col = INK2, col.axis = INK2); axis(2, at = c(1, 3, 10, 30), col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
panel_lab("A")
gg <- 1:10
o1 <- tapply(s$e_ins[s$cell == "H5_k25_ar03"], s$g[s$cell == "H5_k25_ar03"], mean); p1 <- tapply(s$pred[s$cell == "H5_k25_ar03"], s$g[s$cell == "H5_k25_ar03"], mean)
o2 <- tapply(s$e_ins[s$cell == "prev08_k25"], s$g[s$cell == "prev08_k25"], mean); p2 <- tapply(s$pred[s$cell == "prev08_k25"], s$g[s$cell == "prev08_k25"], mean)
c1 <- tapply(s$e_cav[s$cell == "H5_k25_ar03"], s$g[s$cell == "H5_k25_ar03"], mean)
plot(gg, o1, ylim = range(c(o1, p1, o2, p2)), type = "n", xlab = "Group (decile)", ylab = "Mean standardised residual", axes = FALSE)
abline(h = pretty(range(c(o1, p1, o2, p2))), col = GRID); abline(h = 0, col = "grey60")
lines(gg, p1, col = COL["cav"], lty = 2, lwd = 1.6); points(gg, o1, pch = 16, col = COL["cav"])
lines(gg, p2, col = COL["ins"], lty = 2, lwd = 1.6); points(gg, o2, pch = 17, col = COL["ins"])
points(gg, c1, pch = 4, col = "black", lwd = 1.4)
axis(1, at = gg, gap.axis = 0, col = INK2, col.axis = INK2); axis(2, col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
legend("topleft", bty = "n", cex = 0.8, pch = c(16, 17, NA, 4), lty = c(NA, NA, 2, NA), col = c(COL["cav"], COL["ins"], "grey30", "black"),
       legend = c("weak correlation", "prevalence 0.08", "predicted", "grouped on the cavity index"))
panel_lab("B")
invisible(dev.off())

# ---- Figure 3: size in every design (E32 = E18 + textbook, E19 stress, E22 lasso) -----------------------------------
sz <- read.csv(file.path(RES, "E32_textbook_on_E18_summary.csv"))
labA <- c(A_p5 = "fixed p", B_p40 = "p/n 0.10", B_p100 = "p/n 0.25", B_p160 = "p/n 0.40", H3_k40_lam50 = "p/n 0.40, light penalty",
          H5_k25_ar03 = "weak correlation", H8_n800_k25 = "n = 800", k10_t3 = "t3 predictors", prev15_k25 = "prevalence 0.15",
          prev08_k25 = "prevalence 0.08", prev08_A = "fixed p, prevalence 0.08")
sz <- sz[match(names(labA), sz$cell), ]
st <- read.csv(file.path(RES, "E19_cavity_stress_summary.csv"))
labB <- c(n1600_k25 = "n = 1600", n3200_k25 = "n = 3200", strong4_k25 = "strong signal", weak1_k25 = "weak signal", lam10_k25 = "light penalty",
          lam5000_k25 = "heavy penalty", k50 = "p/n 0.50", binary_k25 = "binary predictors", t3_k25 = "t3 predictors", G5_k25 = "G = 5",
          G20_k25 = "G = 20", prev05_strong3_k25 = "prevalence 0.05", A_n5000 = "fixed p, n = 5000")
st <- st[match(names(labB), st$cell), ]
la <- read.csv(file.path(RES, "E22_cavity_lasso_summary.csv"))
labC <- c(lasso_dense_k25_l02 = "lasso, dense", lasso_dense_k25_l05 = "lasso, dense, heavier", enet5_dense_k25_l05 = "elastic net, dense",
          lasso_sparse5_k25 = "lasso, sparse", lasso_sparse5_k50 = "lasso, sparse, p/n 0.5", lasso_sparse5_p2n = "lasso, sparse, p = 2n",
          enet5_sparse10_p2n = "elastic net, p = 2n", lasso_sparse5_prev10 = "lasso, prevalence 0.10")
la <- la[match(names(labC), la$cell), ]
xl <- c(0.015, 1); tk <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
dotpanel <- function(lab, series, title_letter) {
  k <- length(lab); yy <- rev(seq_len(k))
  plot(1, type = "n", log = "x", xlim = xl, ylim = c(0.5, k + 0.5), xlab = "Rejection rate under a correct model", ylab = "", axes = FALSE)
  rect(0.031, 0, 0.069, k + 1, col = "grey92", border = NA); abline(v = 0.05, col = "grey55", lty = 2); abline(h = yy, col = GRID)
  for (sv in series) points(pmax(sv$v, 0.016), yy, pch = sv$pch, col = sv$col, bg = sv$col, cex = 1.1)
  axis(1, at = tk, labels = c("0.02", "0.05", "0.1", "0.2", "0.5", "1"), col = INK2, col.axis = INK2)
  axis(2, at = yy, labels = unname(lab), tick = FALSE, cex.axis = 0.82, col.axis = "grey15")
  panel_lab(title_letter)
}
pdf(file.path(FIG, "fig3_size.pdf"), width = 7.2, height = 8.6, pointsize = 10)
layout(matrix(1:3, 3), heights = c(11.6, 13.6, 8.6)); par(cex = 1, mgp = c(2.3, 0.6, 0), las = 1)
par(mar = c(3.6, 11, 2.6, 1))
dotpanel(labA, list(list(v = sz$textbook, pch = 22, col = COL["txt"]), list(v = sz$insample, pch = 24, col = COL["ins"]), list(v = sz$cavity, pch = 21, col = COL["cav"])), "A")
# legend above the panel: inside the plot it hid the textbook squares at a rejection rate of 1
legend("top", inset = c(0, -0.16), xpd = NA, horiz = TRUE, bty = "n", cex = 0.85, pch = c(21, 24, 22), pt.bg = COL, col = COL,
       legend = c("cavity test", "in-sample closed form", "textbook test"))
par(mar = c(3.6, 11, 1.6, 1))
dotpanel(labB, list(list(v = st$p_a, pch = 21, col = COL["cav"])), "B")
dotpanel(labC, list(list(v = la$naive, pch = 22, col = COL["txt"]), list(v = la$p_a, pch = 21, col = COL["cav"])), "C")
invisible(dev.off())

# ---- Figure 4: power against larger departures (E30, size 1.5) ------------------------------------------------------
pw <- read.csv(file.path(RES, "E30_power_strong_summary.csv")); pw <- pw[pw$rel == 1.5, ]
cellsP <- c(A_p5 = "fixed p", B_p40 = "p/n 0.10", B_p100 = "p/n 0.25"); deps <- c(cubic = "cubic", quadratic = "quadratic", interaction = "interaction")
pdf(file.path(FIG, "fig4_power.pdf"), width = 7.2, height = 3.2, pointsize = 10)
# legend in its own strip under the panels: inside panel 1 it covered the quadratic points
layout(matrix(c(1, 2, 3, 4, 4, 4), 2, byrow = TRUE), heights = c(1, 0.13)); par(cex = 1, mar = c(2.6, 3.6, 1.8, 0.6), mgp = c(2.3, 0.6, 0), las = 1)
for (k in seq_along(cellsP)) {
  d <- pw[pw$cell == names(cellsP)[k], ]; d <- d[match(names(deps), d$dep), ]
  plot(c(0.5, 3.5), c(0, 1), type = "n", xlab = "", ylab = if (k == 1) "Power" else "", axes = FALSE)
  abline(h = seq(0, 1, 0.25), col = GRID); abline(h = 0.05, col = "grey55", lty = 2)
  xo <- c(-0.21, -0.07, 0.07, 0.21)
  points(1:3 + xo[1], d$cav_hl, pch = 21, bg = COL["cav"], col = COL["cav"], cex = 1.15)
  points(1:3 + xo[2], d$cav_e3, pch = 21, bg = "white", col = COL["cav"], cex = 1.15, lwd = 1.5)
  points(1:3 + xo[3], d$calm_hl, pch = 24, bg = COL["ins"], col = COL["ins"], cex = 1.05)
  points(1:3 + xo[4], d$calm_edge, pch = 24, bg = "white", col = COL["ins"], cex = 1.05, lwd = 1.5)
  axis(1, at = 1:3, labels = FALSE, tick = FALSE); mtext(unname(deps), side = 1, at = 1:3, line = 0.8, cex = 0.72); axis(2, at = seq(0, 1, 0.25), col = INK2, col.axis = INK2, cex.axis = 0.85)
  mtext(cellsP[k], side = 3, line = 0.3, cex = 0.9, col = "grey15")
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", horiz = TRUE, bty = "n", cex = 0.85, pch = c(21, 21, 24, 24), pt.bg = c(COL["cav"], "white", COL["ins"], "white"), col = COL[c("cav", "cav", "ins", "ins")],
       pt.lwd = 1.5, legend = c("cavity, decile", "cavity, smooth", "in-sample, decile", "in-sample, smooth"))
invisible(dev.off())

# ---- Figure 5: glaucoma ---------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
Xg <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); yg <- as.numeric(GlaucomaM$Class == "glaucoma")
cv <- cav_gof_fast(Xg, yg, 295.3, G = 10)
Xg1 <- cbind(1, Xg); bg <- cf_ridge(Xg1, yg, 295.3, diag(c(0, rep(1, ncol(Xg))))); pih <- cf_expit(drop(Xg1 %*% bg))
idx <- split(seq_along(yg), cv$groups)
ob <- sapply(idx, function(I) mean(yg[I])); ex <- sapply(idx, function(I) mean(cv$m[I])); nn <- sapply(idx, length); fi <- sapply(idx, function(I) mean(pih[I]))
se <- sqrt(ex * (1 - ex) / nn)
pdf(file.path(FIG, "fig5_glaucoma.pdf"), width = 4.6, height = 4.0, pointsize = 10)
par(cex = 1, mar = c(4, 4.3, 1, 1), mgp = c(2.4, 0.7, 0), las = 1)
plot(c(0, 1), c(0, 1), type = "n", xlab = "Expected proportion", ylab = "Observed proportion", axes = FALSE)
abline(h = seq(0, 1, .2), v = seq(0, 1, .2), col = GRID); abline(0, 1, lty = 2, col = "grey50")
arrows(ex, pmax(ex - 1.96 * se, 0), ex, pmin(ex + 1.96 * se, 1), angle = 90, code = 3, length = 0.03, col = adjustcolor(COL["cav"], 0.45))
points(fi, ob, pch = 24, bg = "white", col = COL["ins"], cex = 1.05, lwd = 1.4)
points(ex, ob, pch = 21, bg = COL["cav"], col = "white", cex = 1.3)
text(ex[5], ob[5], "decile 5", pos = 2, cex = 0.8, col = "grey20")
axis(1, col = INK2, col.axis = INK2); axis(2, col = INK2, col.axis = INK2); box(bty = "l", col = INK2)
legend("topleft", bty = "n", cex = 0.8, pch = c(21, 24), pt.bg = c(COL["cav"], "white"), col = c("white", COL["ins"]), pt.cex = c(1.3, 1),
       legend = c("against the recalibrated cavity curve", "against the shrunk ridge probabilities"))
invisible(dev.off())
cat("figures written\n")
