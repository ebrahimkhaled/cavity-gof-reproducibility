# Cluster correction for the cavity test (E35d). With strongly clustered fixed predictors (repeated or nearly
# repeated rows), observations in one cluster share the error of their index, so their residuals are positively
# correlated and the grouped residuals are overdispersed relative to Omega. The correction adds the within-cluster
# cross-products of the projected residuals to Omega:
#   Omega_c = Omega + sum_c sum_{i != j in c} k_i k_j' r_i r_j,   k_i = a_i - U_r F_r^{-1} z_i,
# where a_i is the group indicator of i scaled by V_g^{-1/2}, z_i = (1, xi_i) and r_i = y_i - m_i. Singleton
# clusters add nothing, so without near-twins the test is unchanged.

# near-twin clusters: connected components of the graph with an edge where the leave-one-out cosine
# B_kj / sqrt(B_kk B_jj) exceeds `cut` (B = X1 M^{-1} X1' from the fit)
cav_clusters <- function(X, beta, lambda, cut = 0.9) {
  X1 <- cbind(1, X); pi <- cf_expit(drop(X1 %*% beta)); w <- pi * (1 - pi)
  L <- X1 %*% chol2inv(chol(crossprod(X1, w * X1) + lambda * diag(c(0, rep(1, ncol(X))))))
  q <- rowSums(L * X1); n <- nrow(X); cl <- seq_len(n)
  find <- function(i) { while (cl[i] != i) { cl[i] <<- cl[cl[i]]; i <- cl[i] }; i }
  for (k in seq_len(n - 1)) {                              # row by row, so no n x n matrix is stored
    bk <- drop(X1[(k + 1):n, , drop = FALSE] %*% L[k, ]) / sqrt(q[k] * q[(k + 1):n])
    for (j in k + which(bk > cut)) { a <- find(k); b <- find(j); if (a != b) cl[b] <- a }
  }
  vapply(seq_len(n), find, 1L)
}

# share of observations with a near-twin: the diagnostic reported before the test
cav_twin_share <- function(cl) mean(tabulate(cl)[cl] > 1)

cav_from_index_cluster <- function(xi, y, cl, G = 10) {
  n <- length(y); rc <- cav_recal_fit(xi, y)
  m <- pmin(pmax(cf_expit(rc[1] + rc[2] * xi), 1e-9), 1 - 1e-9); v <- m * (1 - m)
  groups <- pmin(ceiling(rank(xi, ties.method = "random") / (n / G)), G)
  idx <- split(seq_len(n), groups); G <- length(idx)
  Vg <- vapply(idx, function(I) sum(v[I]), 0.0)
  r <- y - m; e <- vapply(idx, function(I) sum(r[I]), 0.0) / sqrt(Vg)
  Z <- cbind(1, xi); Fr <- crossprod(Z, v * Z)
  Ur <- t(vapply(idx, function(I) colSums(v[I] * Z[I, , drop = FALSE]), numeric(2))) / sqrt(Vg)
  Om <- diag(G) - Ur %*% solve(Fr, t(Ur))
  K <- -Z %*% solve(Fr, t(Ur)); K[cbind(seq_len(n), groups)] <- K[cbind(seq_len(n), groups)] + 1 / sqrt(Vg[groups])   # n x G, rows k_i
  KR <- K * r
  big <- which(tabulate(cl)[cl] > 1)
  Ex <- matrix(0, G, G)
  if (length(big)) {
    s <- rowsum(KR[big, , drop = FALSE], cl[big])           # per-cluster sums of k_i r_i
    Ex <- crossprod(s) - crossprod(KR[big, , drop = FALSE]) # off-diagonal (i != j) part only
  }
  Oc <- Om + Ex; ev <- eigen(Oc, symmetric = TRUE); Oc <- ev$vectors %*% (pmax(ev$values, 0) * t(ev$vectors))
  S <- sum(e^2)
  list(S = S, e = e, p_chi2 = pchisq(S, G - 2, lower.tail = FALSE), p_omega = cav_pvalue(S, Om), p_cluster = cav_pvalue(S, Oc),
       trace_ratio = sum(diag(Oc)) / sum(diag(Om)))
}
