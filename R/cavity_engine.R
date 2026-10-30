# =====================================================================
# cavity_engine.R -- a grouped goodness-of-fit test built ENTIRELY on the cavity (leave-one-out) index.
#
# Every earlier attempt at removing the self-selection term mixed a cavity quantity with an in-sample one:
#   E3e grouped on the leave-one-out index but kept the full-fit residual y - pi_hat  (anti-selection);
#   E11 differenced group sums with and without the one-step displacement            (over-corrects);
#   E3d predicted the raw residual mean by mean-field and subtracted mu_hat           (1% accuracy needed).
# Here the grouping, the fitted probability and the variance all come from the same cavity index
#   xi_i = x_i' beta_hat_(-i)   (one Newton step from the ridge fit, exact derivative, no refits),
# recalibrated along itself by a two-parameter logistic fit  m_i = expit(alpha + b xi_i).
# Under the structural null the cavity index is independent of y_i given the true index, so
# E(y_i | xi_i) is a smooth function of xi_i that the logistic recalibration captures; there is no
# self-selection to correct, because y_i neither moves its own group nor its own fitted value.
#
# Covariance of the standardised group residuals, by the delta method in y with the EXACT Jacobian
# of the ridge fit (d beta_hat / d y = M^{-1} X'):
#   d xi_k / d y_j = x_k' M^{-1} x_j / (1 - H_k)   (j != k),   0   (j = k: the cavity property)
#   J = I - D_v Z F_r^{-1} Z' (I - b D_v L) - b D_v L,   L = off-diagonal cavity Jacobian
#   Cov(y - m) ~= J diag(v) J'.
# =====================================================================
suppressPackageStartupMessages(library(CompQuadForm))
if (!exists("cf_expit")) cf_expit <- function(z) 1 / (1 + exp(-z))
# Gauss-Hermite nodes and weights for a STANDARD NORMAL (weights sum to one); self-contained copy
if (!exists("statmod_gh")) statmod_gh <- function(n) {
  i <- seq_len(n - 1); J <- matrix(0, n, n); J[cbind(i, i + 1)] <- sqrt(i); J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2)
}
# recalibration fit with a degeneracy guard: a constant or perfectly separating cavity index (e.g. a lasso fit
# with an empty active set) makes the recalibration meaningless; stop rather than return a silent number.
cav_recal_fit <- function(xi, y) {
  if (sd(xi) < 1e-8) stop("cavity index is constant: the fit is degenerate (Assumption 2(b) fails)")
  f <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial()))
  if (!f$converged || any(!is.finite(f$coefficients)) || abs(f$coefficients[2]) > 50)
    stop("recalibration did not converge or separates: the fit is degenerate (Assumption 2(b) fails)")
  f$coefficients
}
if (!exists("cf_ridge")) cf_ridge <- function(X1, y, lambda, Dpen, tol = 1e-10, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    bn <- drop(solve(crossprod(X1, w * X1) + lambda * Dpen, crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  beta
}

cav_pvalue <- function(S, Sigma) {
  lam <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values; lam <- lam[lam > 1e-9]
  if (!length(lam)) return(NA_real_)
  d <- CompQuadForm::davies(S, lambda = lam, lim = 1e5, acc = 1e-7); p <- d$Qq
  if (d$ifault != 0 || is.na(p) || p < 0 || p > 1) p <- CompQuadForm::imhof(S, lambda = lam)$Qq
  min(max(p, 0), 1)
}

# the cavity pieces from one ridge fit
cav_pieces <- function(X, y, lambda, G = 10, beta = NULL) {
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, ncol(X))))
  if (is.null(beta)) beta <- cf_ridge(X1, y, lambda, Dpen)
  eta <- drop(X1 %*% beta); pi <- pmin(pmax(cf_expit(eta), 1e-9), 1 - 1e-9); w <- pi * (1 - pi)
  M <- crossprod(X1, w * X1) + lambda * Dpen
  B <- X1 %*% solve(M, t(X1))                      # n x n, B_kj = x_k' M^{-1} x_j
  q <- diag(B); H <- w * q                         # leverages, in [0, 1)
  psi <- y - pi
  xi <- eta - q / (1 - H) * psi                     # cavity index (one Newton step, exact derivative)
  rc <- cav_recal_fit(xi, y)
  m <- pmin(pmax(cf_expit(rc[1] + rc[2] * xi), 1e-9), 1 - 1e-9); v <- m * (1 - m)
  groups <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G)
  list(X1 = X1, B = B, H = H, q = q, xi = xi, eta = eta, pi = pi, w = w, alpha = rc[1], b = rc[2],
       m = m, v = v, groups = groups, beta = beta, n = n, G = G)
}

# Fast version for large n: no n x n matrix. Reference (a) in closed form:
#   J D_v J' = D_v - D_v Z F_r^{-1} Z' D_v   for J = I - D_v Z F_r^{-1} Z',  so
#   Omega_a = I_G - U_r F_r^{-1} U_r',   U_r = A D_v Z   (G x 2): the classical Hosmer-Lemeshow projection
#   for the two-parameter recalibration model.
cav_gof_fast <- function(X, y, lambda, G = 10, beta = NULL, basis_degree = 0) {
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, ncol(X))))
  if (is.null(beta)) beta <- cf_ridge(X1, y, lambda, Dpen)
  eta <- drop(X1 %*% beta); pi <- pmin(pmax(cf_expit(eta), 1e-9), 1 - 1e-9); w <- pi * (1 - pi)
  M <- crossprod(X1, w * X1) + lambda * Dpen
  q <- rowSums((X1 %*% solve(M)) * X1); H <- w * q
  xi <- eta - q / (1 - H) * (y - pi)
  rc <- cav_recal_fit(xi, y)
  m <- pmin(pmax(cf_expit(rc[1] + rc[2] * xi), 1e-9), 1 - 1e-9); v <- m * (1 - m)
  groups <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), groups); G <- length(idx)
  Vg <- vapply(idx, function(I) sum(v[I]), 0.0)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0.0) / sqrt(Vg)
  Z <- cbind(1, xi); Fr <- crossprod(Z, v * Z)
  Ur <- t(vapply(idx, function(I) colSums(v[I] * Z[I, , drop = FALSE]), numeric(2))) / sqrt(Vg)
  Om <- diag(G) - Ur %*% solve(Fr, t(Ur))
  S <- sum(e^2)
  out <- list(S = S, e = e, p_chi2 = pchisq(S, G - 2, lower.tail = FALSE), p_a = cav_pvalue(S, Om),
              tr_a = sum(diag(Om)), b = unname(rc[2]), alpha = unname(rc[1]), Hbar = mean(H), xi = xi, m = m, groups = groups)
  if (basis_degree > 0) {                                     # smooth-basis version: polynomials in the group-mean cavity index
    xbar <- vapply(idx, function(I) mean(xi[I]), 0.0)
    Pz <- as.matrix(stats::poly(xbar, basis_degree))
    zr <- crossprod(Pz, e); S_E <- as.numeric(t(zr) %*% solve(crossprod(Pz)) %*% zr)
    e2 <- eigen(Om, symmetric = TRUE); Sh <- e2$vectors %*% (sqrt(pmax(e2$values, 0)) * t(e2$vectors))
    Pp <- Pz %*% solve(crossprod(Pz), t(Pz))
    out$S_E <- S_E; out$p_E <- cav_pvalue(S_E, Sh %*% Pp %*% Sh)
  }
  out
}

# ---- the same test after a LASSO / ELASTIC-NET fit (glmnet). The exact cavity law holds for any estimator built
# from D_(-i); only the leave-one-out approximation changes: on the active set E (Rad & Maleki 2020, non-smooth case)
#   xi_i = eta_i - q_i / (1 - w_i q_i) (y_i - pi_i),   q_i = x_iE' M_E^{-1} x_iE,
#   M_E = X_E' W X_E + n lambda (1 - alpha) D_E     (glmnet objective: mean deviance / 2 + lambda P_alpha(beta)).
# X is standardised once by the caller and glmnet is run with standardize = FALSE, so the penalty acts on the X given.
cav_gof_glmnet <- function(X, y, lambda_g, alpha = 1, G = 10, basis_degree = 0) {
  X <- as.matrix(X); n <- nrow(X)
  fit <- glmnet::glmnet(X, y, family = "binomial", alpha = alpha, lambda = lambda_g, standardize = FALSE, thresh = 1e-12, maxit = 1e6)
  b <- as.numeric(coef(fit, s = lambda_g)); E <- which(b[-1] != 0)
  X1 <- cbind(1, X); eta <- drop(X1 %*% b); pi <- pmin(pmax(cf_expit(eta), 1e-9), 1 - 1e-9); w <- pi * (1 - pi)
  XE <- cbind(1, X[, E, drop = FALSE]); DE <- diag(c(0, rep(1, length(E))), length(E) + 1)
  ME <- crossprod(XE, w * XE) + n * lambda_g * (1 - alpha) * DE
  q <- rowSums((XE %*% solve(ME)) * XE); H <- w * q
  xi <- eta - q / (1 - H) * (y - pi)
  out <- cav_from_index(xi, y, G, basis_degree)
  out$nnz <- length(E); out$Hbar <- mean(H); out$naive <- naive_hl(eta, y, G)
  out
}

# shared back end: recalibrate y on a given index, group on it, classical projection reference
cav_from_index <- function(xi, y, G = 10, basis_degree = 0) {
  n <- length(y)
  rc <- cav_recal_fit(xi, y)
  m <- pmin(pmax(cf_expit(rc[1] + rc[2] * xi), 1e-9), 1 - 1e-9); v <- m * (1 - m)
  groups <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), groups); G <- length(idx)
  Vg <- vapply(idx, function(I) sum(v[I]), 0.0)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0.0) / sqrt(Vg)
  Z <- cbind(1, xi); Fr <- crossprod(Z, v * Z)
  Ur <- t(vapply(idx, function(I) colSums(v[I] * Z[I, , drop = FALSE]), numeric(2))) / sqrt(Vg)
  Om <- diag(G) - Ur %*% solve(Fr, t(Ur)); S <- sum(e^2)
  out <- list(S = S, e = e, p_chi2 = pchisq(S, G - 2, lower.tail = FALSE), p_a = cav_pvalue(S, Om), tr_a = sum(diag(Om)),
              b = unname(rc[2]), alpha = unname(rc[1]))
  if (basis_degree > 0) {
    xbar <- vapply(idx, function(I) mean(xi[I]), 0.0); Pz <- as.matrix(stats::poly(xbar, basis_degree))
    zr <- crossprod(Pz, e); out$S_E <- as.numeric(t(zr) %*% solve(crossprod(Pz)) %*% zr)
    e2 <- eigen(Om, symmetric = TRUE); Sh <- e2$vectors %*% (sqrt(pmax(e2$values, 0)) * t(e2$vectors))
    out$p_E <- cav_pvalue(out$S_E, Sh %*% (Pz %*% solve(crossprod(Pz), t(Pz))) %*% Sh)
  }
  out
}

# EXACT recalibration family of Proposition 1:  Lambda(xi; a, b, tau) = E pi(a + b xi + tau Z), Z ~ N(0,1),
# fitted by maximum likelihood (Gauss-Hermite, 24 nodes). Reference = classical projection for a q = 3 model:
#   Omega = I - U F^{-1} U',  U = A dm/dtheta,  F = sum_i (dm_i/dtheta)(dm_i/dtheta)' / v_i.
cav_lognormal <- function(xi, y, G = 10, nodes = 24, start = NULL) {
  n <- length(y); gh <- statmod_gh(nodes)
  lam <- function(th) { a <- th[1]; b <- th[2]; tau <- exp(th[3])
    P <- 0; for (k in seq_len(nodes)) P <- P + gh$w[k] * cf_expit(a + b * xi + tau * gh$x[k]); pmin(pmax(P, 1e-9), 1 - 1e-9) }
  nll <- function(th) { m <- lam(th); -sum(y * log(m) + (1 - y) * log(1 - m)) }
  if (is.null(start)) { rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients); start <- c(rc, log(0.5)) }
  op <- optim(start, nll, method = "L-BFGS-B", lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)))
  th <- op$par; m <- lam(th); v <- m * (1 - m)
  a <- th[1]; b <- th[2]; tau <- exp(th[3])
  D <- matrix(0, n, 3)                                   # dm/d(a, b, log tau)
  for (k in seq_len(nodes)) { pk <- cf_expit(a + b * xi + tau * gh$x[k]); dk <- gh$w[k] * pk * (1 - pk)
    D[, 1] <- D[, 1] + dk; D[, 2] <- D[, 2] + dk * xi; D[, 3] <- D[, 3] + dk * gh$x[k] * tau }
  groups <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G); idx <- split(seq_len(n), groups); G <- length(idx)
  Vg <- vapply(idx, function(I) sum(v[I]), 0.0)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0.0) / sqrt(Vg)
  if (tau <= 1.001e-3) D <- D[, 1:2, drop = FALSE]    # tau on its lower bound: the noise column is ~0 and F is singular;
                                                     # drop it (boundary case; the fit is then the logistic one)
  Fi <- crossprod(D, D / v)
  U <- t(vapply(idx, function(I) colSums(D[I, , drop = FALSE]), numeric(ncol(D)))) / sqrt(Vg)
  Om <- diag(G) - U %*% solve(Fi, t(U)); S <- sum(e^2)
  list(S = S, e = e, p_a = cav_pvalue(S, Om), p_chi2 = pchisq(S, G - 3, lower.tail = FALSE), tr = sum(diag(Om)),
       a = a, b = b, tau = tau)
}

# ---- SAMPLE-SPLIT version (fully proved): fit on half A, index on half B is x_B' beta_A EXACTLY (no approximation);
# conditionally on A the pairs (xi_i, y_i), i in B, are i.i.d. with the exact logistic-normal curve of Proposition 2.
# lambda is rescaled to the half sample so that lambda / n is unchanged. Cross-fit: swap the halves and combine the
# two p-values by the Cauchy combination (valid under arbitrary dependence, asymptotically).
cav_split_index <- function(X, y, idxA, fit_fun) {
  bA <- fit_fun(X[idxA, , drop = FALSE], y[idxA])
  drop(cbind(1, X[-idxA, , drop = FALSE]) %*% bA)
}
cav_split_gof <- function(X, y, lambda = NULL, G = 10, lambda_g = NULL, alpha = 1, seed = 1, q = 2) {
  X <- as.matrix(X); n <- nrow(X)
  fit_fun <- if (is.null(lambda_g)) {
    function(Xa, ya) cf_ridge(cbind(1, Xa), ya, lambda * nrow(Xa) / n, diag(c(0, rep(1, ncol(Xa)))))
  } else {
    function(Xa, ya) as.numeric(coef(glmnet::glmnet(Xa, ya, family = "binomial", alpha = alpha, lambda = lambda_g,
                                                    standardize = FALSE, thresh = 1e-12, maxit = 1e6), s = lambda_g))
  }
  set.seed(seed); idxA <- sort(sample.int(n, floor(n / 2))); idxB <- setdiff(seq_len(n), idxA)
  one_half <- function(train, test) {
    xi <- cav_split_index(X, y, train, fit_fun)
    if (q == 2) { r <- cav_from_index(xi, y[test], G); r$p_chi2 } else { r <- cav_lognormal(xi, y[test], G); pchisq(r$S, G - 2, lower.tail = FALSE) }
  }
  p1 <- one_half(idxA, idxB); p2 <- one_half(idxB, idxA)
  cauchy <- 0.5 - atan(mean(tan((0.5 - c(p1, p2)) * pi))) / pi
  list(p_split = p1, p_cross = cauchy)
}

# the textbook Hosmer-Lemeshow test on the in-sample fitted probabilities, for comparison (chi2_{G-2})
naive_hl <- function(eta, y, G = 10) {
  n <- length(y); pi <- cf_expit(eta); groups <- pmin(ceiling(rank(pi, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), groups)
  S <- sum(vapply(idx, function(I) (sum(y[I]) - sum(pi[I]))^2 / (length(I) * mean(pi[I]) * (1 - mean(pi[I]))), 0.0))
  pchisq(S, length(idx) - 2, lower.tail = FALSE)
}

cav_gof <- function(X, y, lambda, G = 10, keep = FALSE) {
  pc <- cav_pieces(X, y, lambda, G)
  n <- pc$n; idx <- split(seq_len(n), pc$groups); G <- length(idx)
  Vg <- vapply(idx, function(I) sum(pc$v[I]), 0.0)
  e <- vapply(idx, function(I) sum(y[I] - pc$m[I]), 0.0) / sqrt(Vg)
  S <- sum(e^2)
  C <- matrix(0, G, n); for (g in seq_len(G)) C[g, idx[[g]]] <- 1
  A <- C / sqrt(Vg)                                  # G x n aggregation
  Z <- cbind(1, pc$xi); Fr <- crossprod(Z, pc$v * Z)
  # (a) classical: recalibration parameters estimated, cavity index treated as fixed
  P <- Z %*% solve(Fr, t(Z * pc$v))                  # n x n:  Z F_r^{-1} Z' D_v
  Ja <- diag(n) - pc$v * Z %*% solve(Fr, t(Z))       # I - D_v Z F^-1 Z'
  Om_a <- A %*% (Ja %*% (pc$v * t(Ja))) %*% t(A)
  # (b) plus the dependence of every cavity index on the OTHER responses (exact ridge Jacobian)
  L <- pc$B / (1 - pc$H); diag(L) <- 0               # d xi_k / d y_j
  bL <- pc$b * (pc$v * L)                            # b D_v L
  Jb <- diag(n) - pc$v * Z %*% solve(Fr, t(Z) %*% (diag(n) - bL)) - bL
  Om_b <- A %*% (Jb %*% (pc$v * t(Jb))) %*% t(A)
  out <- list(S = S, e = e,
              p_I = cav_pvalue(S, diag(G)),
              p_chi2 = pchisq(S, G - 2, lower.tail = FALSE),
              p_a = cav_pvalue(S, Om_a), p_b = cav_pvalue(S, Om_b),
              tr_a = sum(diag(Om_a)), tr_b = sum(diag(Om_b)),
              b = pc$b, alpha = pc$alpha, Hbar = mean(pc$H))
  if (keep) { out$Om_a <- Om_a; out$Om_b <- Om_b; out$pc <- pc }
  out
}
