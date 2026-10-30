# Exact check (Davies) that chi2_{G-3} + w1 chi2_1 + w2 chi2_1 is conservative against the chi2_{G-2} critical value
# over the whole region w1, w2 in [0, 1], w1 + w2 <= 1, for G = 5, 10, 20 and alpha = 0.01, 0.05, 0.10.
suppressPackageStartupMessages(library(CompQuadForm))
grid <- seq(0, 1, by = 0.01)
for (G in c(5, 10, 20)) for (a in c(0.01, 0.05, 0.10)) {
  crit <- qchisq(1 - a, G - 2); worst <- 0; arg <- c(NA, NA)
  for (w1 in grid) for (w2 in grid[grid <= 1 - w1 + 1e-12]) {
    lam <- c(rep(1, G - 3), w1, w2); lam <- lam[lam > 1e-12]
    pr <- davies(crit, lam, acc = 1e-9, lim = 1e6)$Qq
    if (pr > worst) { worst <- pr; arg <- c(w1, w2) }
  }
  cat(sprintf("G %2d alpha %.2f : max size over w1 + w2 <= 1 is %.5f at (w1, w2) = (%.2f, %.2f)\n", G, a, worst, arg[1], arg[2]))
}
