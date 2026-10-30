# How far is the logistic-normal curve Lambda(xi) = E pi(b xi + tau Z) from the best two-parameter logistic curve,
# over the distribution of the cavity index? Gives the bias of the q = 2 recalibration, and the n at which it
# starts to matter for a G = 10 test (noncentrality = n * sum_g p_g * (gap_g)^2 / v_g  >= 1).
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
source(file.path(ROOT, "R", "closedform_engine.R"))
gh <- statmod_gh(60); expit <- function(z) 1 / (1 + exp(-z))
for (cfg in list(c(3.00, 1.62, 0.57), c(5.93, 1.66, 0.29), c(2, 1, 0.8), c(2, 2.5, 0.8), c(1.5, 4, 1))) {
  b <- cfg[1]; tau <- cfg[2]; sdxi <- cfg[3]
  x <- qnorm((1:4000 - 0.5) / 4000) * sdxi                     # cavity index distribution (Gaussian, mean 0)
  L <- rowSums(sapply(seq_along(gh$x), function(k) gh$w[k] * expit(b * x + tau * gh$x[k])))
  fit <- suppressWarnings(glm(L ~ x, family = quasibinomial())); Lg <- fitted(fit)   # best logistic in KL sense
  g <- ceiling((1:4000) / 400); gap <- tapply(L - Lg, g, mean); v <- tapply(L * (1 - L), g, mean)
  nc_per_n <- sum(0.1 * gap^2 / v)
  cat(sprintf("b %.2f tau %.2f sd(xi) %.2f | max |Lambda - logistic| = %.4f | decile gaps max %.4f | n for noncentrality 1: %8.0f\n",
              b, tau, sdxi, max(abs(L - Lg)), max(abs(gap)), 1 / nc_per_n))
}
