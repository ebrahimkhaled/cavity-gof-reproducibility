# Eigenvalues of Omega = I - U F^{-1} U' for the logistic-normal (q = 3) recalibration, at the true parameters.
# Question: is there an exact zero eigenvalue (as for the logistic curve), and do the non-unit weights sum to <= 1?
# Population quantities by quadrature over xi ~ N(m, t^2), decile cells; no Monte Carlo.
herm <- function(k) {           # Gauss-Hermite (probabilists') nodes and weights via Golub-Welsch
  J <- matrix(0, k, k); for (i in 1:(k - 1)) J[i, i + 1] <- J[i + 1, i] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2)
}
q <- herm(40)
expit <- function(u) 1 / (1 + exp(-u))
lam_parts <- function(x, a, b, tau) {
  L <- 0; da <- 0; dt <- 0
  for (k in seq_along(q$x)) { p <- expit(a + b * x + tau * q$x[k]); d <- q$w[k] * p * (1 - p)
    L <- L + q$w[k] * p; da <- da + d; dt <- dt + d * q$x[k] }
  list(L = L, D = cbind(da, x * da, dt))        # d/da, d/db, d/dtau
}
omega_eigen <- function(a, b, tau, m = 0, t = 1, G = 10, N = 200000) {
  x <- qnorm((seq_len(N) - 0.5) / N, m, t)       # a fine quantile grid stands in for the integral over xi
  lp <- lam_parts(x, a, b, tau); v <- lp$L * (1 - lp$L)
  g <- pmin(ceiling(seq_len(N) / (N / G)), G)
  Vg <- as.vector(tapply(v, g, sum))
  U <- apply(lp$D, 2, function(col) as.vector(tapply(col, g, sum))) / sqrt(Vg)
  Fi <- crossprod(lp$D, lp$D / v)
  ev <- sort(eigen(diag(G) - U %*% solve(Fi, t(U)), symmetric = TRUE, only.values = TRUE)$values)
  ev
}
cases <- expand.grid(b = c(0.7, 1.5, 3), tau = c(0.25, 0.75, 1.5), a = c(0, -1.5))
for (r in seq_len(nrow(cases))) {
  ev <- omega_eigen(cases$a[r], cases$b[r], cases$tau[r])
  nonunit <- ev[ev < 0.999]
  cat(sprintf("a=%5.1f b=%4.1f tau=%4.2f | smallest three: %.4f %.4f %.4f | unit count %d | sum of non-unit %.3f\n",
              cases$a[r], cases$b[r], cases$tau[r], ev[1], ev[2], ev[3], sum(ev >= 0.999), sum(nonunit)))
}
