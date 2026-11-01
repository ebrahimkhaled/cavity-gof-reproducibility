# Conservativeness of the chi-squared(G-2) critical value when the limit of S is
# chi-squared(G-3) + w0 X0 + w1 X1 + w2 X2 (three non-unit weights, X_j independent chi-squared(1)), w0 + w1 + w2 <= 1.
# Exact tail probabilities by Davies' algorithm on a grid of step 0.05.
suppressPackageStartupMessages(library(CompQuadForm))
step <- 0.05; worst <- NULL
for (G in c(5, 10, 20)) for (alpha in c(0.01, 0.05, 0.10)) {
  crit <- qchisq(1 - alpha, G - 2); mx <- 0; arg <- NA
  for (w0 in seq(0, 1, step)) for (w1 in seq(0, 1 - w0, step)) for (w2 in seq(0, 1 - w0 - w1, step)) {
    lam <- c(rep(1, G - 3), w0, w1, w2); lam <- lam[lam > 0]
    pr <- davies(crit, lam, acc = 1e-9, lim = 1e6)$Qq
    if (pr > mx) { mx <- pr; arg <- c(w0, w1, w2) }
  }
  cat(sprintf("G=%2d alpha=%.2f  max tail prob %.5f at (w0,w1,w2) = (%.2f, %.2f, %.2f)\n", G, alpha, mx, arg[1], arg[2], arg[3]))
}
