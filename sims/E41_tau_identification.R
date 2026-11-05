# E41 -- why the three-parameter fit reports tau-hat ~ 0.49 at n = 5000 and 20000 in design A (fixed p, Sigma = I), where
# the true tau is close to zero (referee m6 / S11; supplement l.222-224 quotes E28's MEDIAN tau-hat).
# FINDING THAT MOTIVATED THIS SCRIPT (profile check, 3 replicates): cav_lognormal() (and the same code in
# ebrahim.gof::cavity.gof) starts L-BFGS-B at tau = 0.5 and, because the likelihood is almost flat in tau, stops next to
# the start: tau-hat = 0.483 at n = 20000 in every replicate, and 0.050 when started at 0.05. So the 'estimate' is the
# starting value. This script measures how much that matters for the TEST, by refitting with a converged optimiser:
#   L-BFGS-B on (a, b, log tau) with the SAME box as cav_lognormal (|a|, |b| <= 20, 1e-3 <= tau <= 20), analytic gradient,
#   factr = 10 and pgtol = 1e-12 (instead of the defaults 1e7 and 0), three starts (tau = 0.05, 0.5, 2), best value kept.
# Data: design A as E28 (Mersenne-Twister reset, set.seed(20261108 + r + n)), n = 500, 5000, 20000, first 300 replicates;
#       design B p/n = 0.25 and 0.40 as E18 (set.seed(20261030 + r)), 500 replicates;
#       the bent design of E25 at n = 5000 (set.seed(20261105 + r)), 300 replicates.
# Output: results/E41_tau_identification_perrep.csv (appended per batch; resumable), summary in E41_tau_identification.log.
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); W <- if (length(args) >= 1) as.integer(args[1]) else 8; RA <- if (length(args) >= 2) as.integer(args[2]) else 300
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
GH <- statmod_gh(24)
lnorm_curve <- function(th, xi) { m <- 0; d <- matrix(0, length(xi), 3); tau <- exp(th[3])
  for (k in 1:24) { pk <- cf_expit(th[1] + th[2] * xi + tau * GH$x[k]); m <- m + GH$w[k] * pk; dk <- GH$w[k] * pk * (1 - pk)
    d[, 1] <- d[, 1] + dk; d[, 2] <- d[, 2] + dk * xi; d[, 3] <- d[, 3] + dk * GH$x[k] * tau }
  list(m = pmin(pmax(m, 1e-12), 1 - 1e-12), d = d) }
# converged logistic-normal fit; returns S, chi2_{G-2} P and the fitted parameters (grouping as cav_lognormal)
lnorm_converged <- function(xi, y, G = 10) {
  nll <- function(th) { m <- lnorm_curve(th, xi)$m; -sum(y * log(m) + (1 - y) * log(1 - m)) }
  gr <- function(th) { c <- lnorm_curve(th, xi); -colSums((y - c$m) / (c$m * (1 - c$m)) * c$d) }
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients)
  fits <- lapply(c(0.05, 0.5, 2), function(t0) { b0 <- rc[2] * sqrt(1 + pi * t0^2 / 8)
    st <- pmin(pmax(c(rc[1] * sqrt(1 + pi * t0^2 / 8), b0, log(t0)), c(-19.9, -19.9, log(1e-3))), c(19.9, 19.9, log(20)))
    optim(st, nll, gr, method = "L-BFGS-B", lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)),
          control = list(factr = 10, pgtol = 1e-12, maxit = 5000)) })
  op <- fits[[which.min(sapply(fits, `[[`, "value"))]]; th <- op$par
  m <- lnorm_curve(th, xi)$m; v <- m * (1 - m)
  g <- pmin(ceiling(rank(xi, ties.method = "first") / (length(y) / G)), G); idx <- split(seq_along(y), g)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0) / sqrt(vapply(idx, function(I) sum(v[I]), 0)); S <- sum(e^2)
  list(S = S, p = pchisq(S, G - 2, lower.tail = FALSE), tau = exp(th[3]), b = th[2], nll = op$value,
       spread = diff(range(sapply(fits, `[[`, "value"))))
}
DES <- list(
  list(name = "A_n500", n = 500, kind = "A", R = RA), list(name = "A_n5000", n = 5000, kind = "A", R = RA), list(name = "A_n20000", n = 20000, kind = "A", R = RA),
  list(name = "B_p100", n = 400, p = 100, kind = "B", R = 500), list(name = "B_p160", n = 400, p = 160, kind = "B", R = 500),
  list(name = "bent_n5000", n = 5000, p = 50, kind = "bent", R = 300))
gen <- function(d, r) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  if (d$kind == "A") { b0 <- c(.5, -.4, .3, -.3, .2); set.seed(20261108 + r + d$n); X <- matrix(rnorm(d$n * 5), d$n)
    return(list(X = X, y = rbinom(d$n, 1, cf_expit(drop(X %*% b0))), lambda = d$n / 5, b0 = b0, Sig = diag(5))) }
  ar <- if (d$kind == "B") 0.7 else 0.9; Sig <- ar^abs(outer(1:d$p, 1:d$p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), d$p)
  b0 <- b0 * 2.31 / sqrt(drop(t(b0) %*% Sig %*% b0))
  set.seed((if (d$kind == "B") 20261030 else 20261105) + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% chol(Sig)
  list(X = X, y = rbinom(d$n, 1, cf_expit(drop(X %*% b0))), lambda = if (d$kind == "B") 416 else 2 * d$n, b0 = b0, Sig = Sig)
}
one <- function(r, d) {
  z <- gen(d, r); c2 <- cav_gof_fast(z$X, z$y, z$lambda, G = 10); c3 <- cav_lognormal(c2$xi, z$y, G = 10)
  cc <- lnorm_converged(c2$xi, z$y, 10)
  nll_def <- { m <- lnorm_curve(c(c3$a, c3$b, log(c3$tau)), c2$xi)$m; -sum(z$y * log(m) + (1 - z$y) * log(1 - m)) }
  m2 <- pmin(pmax(cf_expit(c2$alpha + c2$b * c2$xi), 1e-12), 1 - 1e-12); nll2 <- -sum(z$y * log(m2) + (1 - z$y) * log(1 - m2))
  bh <- cf_ridge(cbind(1, z$X), z$y, z$lambda, diag(c(0, rep(1, ncol(z$X)))))[-1]
  s2 <- drop(t(z$b0) %*% z$Sig %*% z$b0); t2 <- drop(t(bh) %*% z$Sig %*% bh); cb <- drop(t(z$b0) %*% z$Sig %*% bh)
  data.frame(rep = r, design = d$name, q2_p = c2$p_chi2, S2 = c2$S,
             def_tau = c3$tau, def_b = c3$b, def_S = c3$S, def_p = pchisq(c3$S, 8, lower.tail = FALSE), def_nll_excess = nll_def - cc$nll,
             conv_tau = cc$tau, conv_b = cc$b, conv_S = cc$S, conv_p = cc$p, start_spread = cc$spread,
             llgain_conv = nll2 - cc$nll, tau_prop = sqrt(max(s2 - cb^2 / t2, 0)))
}
f_out <- file.path(OUT, "E41_tau_identification_perrep.csv")
cl <- makeCluster(W); clusterExport(cl, c("one", "gen", "GH", "lnorm_curve", "lnorm_converged", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
for (d in DES) {
  have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$design == d$name] } else integer(0)
  todo <- setdiff(seq_len(d$R), have)
  for (batch in split(todo, ceiling(seq_along(todo) / (4 * W)))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
    if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  }
  cat(sprintf("done %s\n", d$name))
}
stopCluster(cl)
res <- read.csv(f_out)
sink(file.path(OUT, "E41_tau_identification.log"), split = TRUE)
for (dn in unique(res$design)) { s <- res[res$design == dn, ]; Rn <- nrow(s)
  q <- function(v) paste(sprintf("%.3f", quantile(v, c(.1, .25, .5, .75, .9))), collapse = " ")
  rr <- function(v) { m <- mean(v < .05); sprintf("%.3f (%.3f)", m, sqrt(m * (1 - m) / length(v))) }
  cat(sprintf("%-11s R %d | true (Prop 2) tau median %.3f\n", dn, Rn, median(s$tau_prop)))
  cat(sprintf("   default fit (start 0.5): tau quantiles 10/25/50/75/90 %s; nll above the converged optimum: median %.4f, max %.4f\n", q(s$def_tau), median(s$def_nll_excess), max(s$def_nll_excess)))
  cat(sprintf("   converged fit:           tau quantiles 10/25/50/75/90 %s; share tau < 0.01 %.2f; log-lik gain over q=2 median %.3f (90%% %.3f)\n", q(s$conv_tau), mean(s$conv_tau < .01), median(s$llgain_conv), quantile(s$llgain_conv, .9)))
  cat(sprintf("   size at .05 (chi2_8): q=2 %s | q=3 default %s | q=3 converged %s | max |S_def - S_conv| %.3f, median %.4f\n",
              rr(s$q2_p), rr(s$def_p), rr(s$conv_p), max(abs(s$def_S - s$conv_S)), median(abs(s$def_S - s$conv_S))))
}
sink()
