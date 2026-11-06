# E50 -- tau from a SEPARATE estimate, not from the curve. See README.
#
# PRE-REGISTRATION (written 2026-11-06 before any E50 result).
# Why: the audit found the d = 3 fit has an interior maximiser only 26% to 54% of the time and sd(tau-hat^2) about 15 at n = 400. A binary y gives
#   tau through the shape of one curve only. But with xi standardised, Prop 2 gives b = gamma rho and tau = gamma sqrt(1 - rho^2), so
#       tau = kappa b,  kappa = sqrt((1 - rho^2) / rho^2),  rho = corr(true index, fitted index) given the fit.
#   This is measurement error: rho^2 is the reliability of the index. With kappa from outside the curve, the fit has two parameters, (a, b).
# Estimates of rho^2 (clamped to [.02, .98]; kappa = 0 is the reference "T0"):
#   R2: rho^2 = corr(y, L)^2 / adjusted R^2. Under Gaussian rows the best linear predictor of y is along the true index, and L is linear in a fresh x,
#       so corr(y, L) = rho x sqrt(R^2_pop). Needs p < n - 2.
#   XF: split-half reliability. Three folds, each fitted alone (penalty scaled by fold size); r = mean correlation of two fold indices on the third fold;
#       Spearman-Brown to full size, rho^2 = 3 r / (1 + 2 r). Any p.
#   OR: oracle, rho^2 from the true beta and Sigma, full-fit coefficients. Diagnostic only, never a candidate.
# Same cells, seeds, index standardisation, G = 10 and chi2_{G-2} as E47c. d3 = the current test, recomputed here.
# PASS for a candidate (XF, or R2 where p < n) on I1 needs ALL of:
#   (a) size in [0.0305, 0.0695] at 5% in A_p5, B_p40, B_p100, B_p160, H5_k25_ar03, lasso_sparse5_p2n (R2: cells with p < n);
#   (b) quadratic power at B_p40_quad1 and B_p100_quad1 not below d3 by more than 2 paired SE;
#   (c) size at most 0.0695 at binary_k25 and t3_k25;
#   (d) stability: sd of the implied tau at B_p100 below one third of its oracle mean, at n = 400.
# Rare or sparse binary cells and k75, ridge_p2n, prev08 are reported; failure there narrows the scope and does not fail the trial.
# If (a), (b) or (c) fails: stop, no tuning. The default logistic index is judged by the same rule and reported beside I1.
suppressPackageStartupMessages({ library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) as.integer(args[1]) else 14; PHASE <- if (length(args) >= 2) args[2] else "run"
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
OUTF <- file.path(ROOT, "results", "E50_tau_perrep.csv")
src <- readLines(file.path(ROOT, "sims", "E47c_standardised.R"))        # cells, gen and loo_linear exactly as E47c
eval(parse(text = src[grep("^SEED0", src):(grep("^# default statistic", src) - 1)]))
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
eval(parse(text = src[grep("^loo_linear", src):grep("^std <-", src)]))

cells_e <- function(xi, y, m, G = 10) { n <- length(y); v <- m * (1 - m); g <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G); idx <- split(seq_len(n), g)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0) / sqrt(vapply(idx, function(I) sum(v[I]), 0)); pchisq(sum(e^2), length(idx) - 2, lower.tail = FALSE) }
d3 <- function(xi, y, nodes = 24) { gh <- statmod_gh(nodes)                 # E43's q3, also returning b and tau
  lam <- function(th) { P <- 0; for (k in seq_len(nodes)) P <- P + gh$w[k] * cf_expit(th[1] + th[2] * xi + exp(th[3]) * gh$x[k]); pmin(pmax(P, 1e-9), 1 - 1e-9) }
  nll <- function(th) { m <- lam(th); -sum(y * log(m) + (1 - y) * log(1 - m)) }
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients)
  op <- optim(c(rc, log(0.5)), nll, method = "L-BFGS-B", lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)))
  c(p = cells_e(xi, y, lam(op$par)), tau = exp(op$par[3])) }
tied <- function(xi, y, kappa, nodes = 24) { if (is.na(kappa)) return(c(p = NA, tau = NA)); gh <- statmod_gh(nodes)   # two parameters, tau = kappa |b|
  lam <- function(th) { P <- 0; for (k in seq_len(nodes)) P <- P + gh$w[k] * cf_expit(th[1] + th[2] * xi + abs(th[2]) * kappa * gh$x[k]); pmin(pmax(P, 1e-9), 1 - 1e-9) }
  nll <- function(th) { m <- lam(th); -sum(y * log(m) + (1 - y) * log(1 - m)) }
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients)
  op <- nlminb(rc, nll, lower = c(-20, -20), upper = c(20, 20))   # CORRECTION after the first run: L-BFGS-B jumped to b = 20 in about 1% of fits, where the clamp in lam() kills the gradient. Only the optimiser changed.
  c(p = cells_e(xi, y, lam(op$par)), tau = kappa * abs(op$par[2])) }
kap <- function(r2) { if (is.na(r2)) return(NA_real_); r2 <- min(max(r2, .02), .98); sqrt((1 - r2) / r2) }
adjR2 <- function(X, y) { n <- length(y); p <- ncol(X); if (p >= n - 2) return(NA_real_); rs <- qr.resid(qr(cbind(1, X)), y)
  1 - (sum(rs^2) / (n - p - 1)) / var(y) }
rho_R2 <- function(L, y, a) { if (is.na(a)) return(NA_real_); c0 <- cor(y, L); if (a <= .001 || c0 <= 0) return(.02); c0^2 / a }
rho_XF <- function(X, y, seed, fit) {                                       # fit(Xa, ya, share) returns coefficients without the intercept
  set.seed(seed); n <- length(y); f <- sample(rep_len(1:3, n)); B <- lapply(1:3, function(k) fit(X[f == k, , drop = FALSE], y[f == k], mean(f == k)))
  if (any(vapply(B, function(b) anyNA(b), TRUE))) return(NA_real_)
  r <- mean(c(cor(X[f == 3, ] %*% B[[1]], X[f == 3, ] %*% B[[2]]), cor(X[f == 2, ] %*% B[[1]], X[f == 2, ] %*% B[[3]]), cor(X[f == 1, ] %*% B[[2]], X[f == 1, ] %*% B[[3]])))
  r <- max(r, 0); 3 * r / (1 + 2 * r) }
lin_b <- function(X, y, pen) { X1 <- cbind(1, X); solve(crossprod(X1) + diag(c(0, rep(pen, ncol(X)))), crossprod(X1, y))[-1] }
rho_OR <- function(b, d) { if (!d$law %in% c("gaussian", "t3", "binary")) return(NA_real_); u <- d$Sc %*% d$b0; w <- d$Sc %*% b; sum(u * w)^2 / (sum(u^2) * sum(w^2)) }

score <- function(tag, L, y, rXF, rR2, rOR) { xi <- std(L); s <- function(f) tryCatch(f, error = function(e) c(p = NA, tau = NA))
  a <- s(d3(xi, y)); t0 <- s(tied(xi, y, 0)); o <- s(tied(xi, y, kap(rOR))); x <- s(tied(xi, y, kap(rXF))); q <- s(tied(xi, y, kap(rR2)))
  v <- c(a[1], t0[1], o[1], x[1], q[1], a[2], o[2], x[2], q[2], rOR, rXF, rR2)
  names(v) <- paste0(tag, "_", c("d3", "T0", "OR", "XF", "R2", "tau_d3", "tau_OR", "tau_XF", "tau_R2", "rho_OR", "rho_XF", "rho_R2")); v }
one <- function(r, d) {
  z <- gen(r, d); X <- z$X; y <- z$y; n <- d$n; sx <- d$seed + r + 5000000; a <- adjR2(X, y); pen <- 4 * d$lamI
  out <- score("I1", loo_linear(X, y, pen), y, rho_XF(X, y, sx, function(Xa, ya, s) lin_b(Xa, ya, pen * s)), rho_R2(loo_linear(X, y, pen), y, a), rho_OR(lin_b(X, y, pen), d))
  if (d$p < n - 2) { L2 <- loo_linear(X, y, 0)
    out <- c(out, score("I2", L2, y, rho_XF(X, y, sx, function(Xa, ya, s) if (ncol(Xa) < nrow(Xa) - 2) lin_b(Xa, ya, 0) else NA), rho_R2(L2, y, a), rho_OR(lin_b(X, y, 0), d))) }
  if (!is.null(d$lambda)) { pc <- cav_pieces(X, y, d$lambda)
    fl <- function(Xa, ya, s) cf_ridge(cbind(1, Xa), ya, d$lambda * s, diag(c(0, rep(1, ncol(Xa)))))[-1]
    out <- c(out, score("full", pc$xi, y, rho_XF(X, y, sx, fl), rho_R2(pc$xi, y, a), rho_OR(pc$beta[-1], d))) }
  cols <- as.vector(outer(c("I1", "I2", "full"), c("d3", "T0", "OR", "XF", "R2", "tau_d3", "tau_OR", "tau_XF", "tau_R2", "rho_OR", "rho_XF", "rho_R2"), paste, sep = "_"))
  res <- setNames(rep(NA_real_, length(cols)), cols); res[names(out)] <- out
  data.frame(rep = r, cell = d$name, t(res))
}
if (PHASE == "smoke") { old <- read.csv(file.path(ROOT, "results", "E47c_std_perrep.csv")); d <- CELLS[[3]]
  for (r in 1:3) { o <- one(r, d); print(round(c(o$I1_d3, old$I1[old$cell == d$name & old$rep == r], o$full_d3, old$full[old$cell == d$name & old$rep == r], o$I1_XF, o$I1_R2, o$I1_OR, o$I1_rho_OR, o$I1_rho_XF, o$I1_rho_R2, o$I1_tau_OR, o$I1_tau_XF, o$I1_tau_R2, o$I1_tau_d3), 3)) }
  for (k in c(3, 7, 12)) print(system.time(one(4, CELLS[[k]]))); quit("no") }
done <- if (file.exists(OUTF)) { x <- read.csv(OUTF); paste(x$cell, x$rep) } else character()
cl <- makeCluster(W); clusterExport(cl, c("std", "gen", "one", "score", "d3", "tied", "kap", "cells_e", "adjR2", "rho_R2", "rho_XF", "rho_OR", "lin_b", "loo_linear", "a_for", "cf_expit", "statmod_gh", "cav_pieces", "cav_recal_fit", "cf_ridge"))
for (d in CELLS) for (st in seq(1, d$R, by = 200)) {
  rs <- st:min(st + 199, d$R); rs <- rs[!paste(d$name, rs) %in% done]; if (!length(rs)) next
  out <- do.call(rbind, parLapply(cl, rs, one, d = d))     # every batch written to disk: resume by re-running
  write.table(out, OUTF, sep = ",", row.names = FALSE, col.names = !file.exists(OUTF), append = file.exists(OUTF))
  r5 <- function(v) round(mean(v < .05, na.rm = TRUE), 3)
  cat(d$name, max(rs), "I1 d3", r5(out$I1_d3), "T0", r5(out$I1_T0), "OR", r5(out$I1_OR), "XF", r5(out$I1_XF), "R2", r5(out$I1_R2), "\n")
}
stopCluster(cl)
