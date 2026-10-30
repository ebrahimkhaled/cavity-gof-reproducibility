# E24 -- the SELECTION term of in-sample grouping, predicted from one fit and checked by simulation.
# Same calibrated cavity residuals eps_i = y_i - m_i, grouped two ways:
#   (cavity)    on deciles of the cavity index xi        -> mean zero (Proposition 1);
#   (in-sample) on deciles of the in-sample index eta_hat = xi + V psi  -> mean m_g (self-selection).
# Prediction (first order in the self-influence V_i psi_i):
#   E sum_{i in g} eps_i ~= sum_i V_i (1 - H_i) v_i [K_h(xi_i - a_{g-1}) - K_h(xi_i - a_g)],   v = m(1-m),
# a_g the in-sample quantile boundaries, K_h a Gaussian kernel (the boundary density, estimated from the same fit).
# Usage: Rscript E24_selection_formula.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 300; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261104
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, prev = .5, indep = FALSE) {
  if (indep) { Sig <- diag(p); b0 <- c(.5, -.4, .3, -.3, .2) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), prev = prev) }
CELLS <- list(mk("A_p5", 500, 5, 100, indep = TRUE), mk("B_p40", 400, 40, 416), mk("B_p100", 400, 100, 416), mk("B_p160", 400, 160, 416),
              mk("H3_k40_lam50", 400, 160, 50), mk("H5_k25_ar03", 400, 100, 416, ar = 0.3), mk("H8_n800_k25", 800, 200, 832),
              mk("prev15_k25", 400, 100, 416, prev = .15), mk("prev08_k25", 400, 100, 416, prev = .08))
G <- 10
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta0 <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta0, d$prev); y <- rbinom(d$n, 1, cf_expit(a + eta0))
  n <- d$n; X1 <- cbind(1, X); Dp <- diag(c(0, rep(1, d$p))); bh <- cf_ridge(X1, y, d$lambda, Dp)
  eta <- drop(X1 %*% bh); pi <- cf_expit(eta); w <- pi * (1 - pi); M <- crossprod(X1, w * X1) + d$lambda * Dp
  q <- rowSums((X1 %*% solve(M)) * X1); H <- w * q; V <- q / (1 - H); xi <- eta - V * (y - pi)
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients); m <- cf_expit(rc[1] + rc[2] * xi); v <- m * (1 - m)
  eps <- y - m
  gc <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G)
  gi <- pmin(ceiling(rank(eta, ties.method = "first") / (n / G)), G)
  sc <- function(g) { idx <- split(seq_len(n), g); vapply(idx, function(I) sum(eps[I]) / sqrt(sum(v[I])), 0.0) }
  e_cav <- sc(gc); e_ins <- sc(gi)
  # prediction: boundaries = in-sample quantiles mapped to the xi scale by xi ~ eta - V psi (use the quantiles of xi itself)
  bnd <- quantile(xi, (1:(G - 1)) / G); h <- 1.06 * sd(xi) * n^(-1 / 5)
  Kd <- function(b) dnorm((xi - b) / h) / h
  num <- numeric(G); for (k in seq_len(G - 1)) { s <- sum(V * (1 - H) * v * Kd(bnd[k])); num[k + 1] <- num[k + 1] + s; num[k] <- num[k] - s }
  Vg <- vapply(split(seq_len(n), gi), function(I) sum(v[I]), 0.0)
  pred <- num / sqrt(Vg)
  data.frame(rep = r, cell = d$name, g = seq_len(G), e_cav = e_cav, e_ins = e_ins, pred = pred, Vbar = mean(V), sdxi = sd(xi))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))))
stopCluster(cl)
write.csv(res, file.path(OUT, "E24_selection_formula_pergroup.csv"), row.names = FALSE)
cat(sprintf("=== E24 | self-selection of in-sample grouping: predicted vs observed | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
out <- do.call(rbind, lapply(unique(res$cell), function(cc) {
  s <- res[res$cell == cc, ]; ob <- tapply(s$e_ins, s$g, mean); pr <- tapply(s$pred, s$g, mean); cv <- tapply(s$e_cav, s$g, mean)
  data.frame(cell = cc, norm2_observed = sum(ob^2), norm2_predicted = sum(pr^2), ratio = sum(pr^2) / sum(ob^2),
             corr = cor(ob, pr), norm2_cavity = sum(cv^2), V2_over_sd2 = mean(s$Vbar^2 / s$sdxi^2)) }))
print(out, digits = 3, row.names = FALSE)
cat("\nPer-group means, in-sample grouping (observed / predicted), first two cells with large terms:\n")
for (cc in c("H5_k25_ar03", "prev08_k25", "B_p100")) { s <- res[res$cell == cc, ]
  cat(cc, ": ", paste(sprintf("%+.2f/%+.2f", tapply(s$e_ins, s$g, mean), tapply(s$pred, s$g, mean)), collapse = "  "), "\n") }
