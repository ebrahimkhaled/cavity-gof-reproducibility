# E44 -- the q = 3 default (logistic-normal recalibration + chi2_{G-2}) next to q = 2 on
#   phase "binary": the six sparse binary / skewed designs of E34 (500 replicates, seeds 20261132 + r) and the rare strong
#                   binary design of E34b (1000 more replicates, seeds 21261199 + r), so 1500 for that design;
#   phase "ads"   : the Internet-Advertisements full-sample placebo of E35 (seeds 20261133 + 1000 + r, 100 replicates).
# Data regenerated exactly: Mersenne-Twister / Inversion / Rejection reset before every set.seed (the originals ran with that
# default). q = 3 uses cav_lognormal()'s fit and computes S only (the E42 route), so a singular information matrix cannot stop
# it. Fingerprint: q = 2 P-values against the stored E34 / E34b / E35 ones (E44_summarise.R).
# Appended per batch, resumable. Usage: Rscript E44_q3_binary_ads.R [workers] [phase: binary | ads | smoke] [R_ads]
suppressPackageStartupMessages({ library(parallel); library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) as.integer(args[1]) else 12; PHASE <- if (length(args) >= 2) args[2] else "binary"
R_ADS <- if (length(args) >= 3) as.integer(args[3]) else 100
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- if (PHASE == "smoke") tempdir() else file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
q3_S <- function(xi, y, G = 10, nodes = 24) {        # cav_lognormal()'s fit, statistic and chi2_{G-2} P only
  n <- length(y); gh <- statmod_gh(nodes)
  lam <- function(th) { P <- 0; for (k in seq_len(nodes)) P <- P + gh$w[k] * cf_expit(th[1] + th[2] * xi + exp(th[3]) * gh$x[k]); pmin(pmax(P, 1e-9), 1 - 1e-9) }
  nll <- function(th) { m <- lam(th); -sum(y * log(m) + (1 - y) * log(1 - m)) }
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients)
  op <- optim(c(rc, log(0.5)), nll, method = "L-BFGS-B", lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)))
  m <- lam(op$par); v <- m * (1 - m); g <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G); idx <- split(seq_len(n), g)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0) / sqrt(vapply(idx, function(I) sum(v[I]), 0))
  c(p = pchisq(sum(e^2), length(idx) - 2, lower.tail = FALSE), tau = exp(op$par[3]))
}
# ---- E34 / E34b designs, copied --------------------------------------------------------------------------------------
.a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda, law, beta, prev, pr = 0.3) list(name = name, n = n, p = p, lambda = lambda, law = law, beta = beta, prev = prev, pr = pr)
CELLS <- list(
  mk("bin5_p100", 400, 100, 100, "binary", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("bin_rare_strong_p100", 400, 100, 100, "binary_rare", c(2.5, 0.6, -0.6, 0.6, -0.6), prev = 0.2),
  mk("bin5_p200", 400, 200, 100, "binary", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("bin5_p800", 400, 800, 100, "binary", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("skew5_p100", 400, 100, 100, "skew", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("mixed5_p100", 400, 100, 100, "mixed", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3))
draw <- function(d, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(seed); n <- d$n; p <- d$p
  X <- switch(d$law, binary = matrix(rbinom(n * p, 1, d$pr), n),
    binary_rare = cbind(rbinom(n, 1, 0.1), matrix(rbinom(n * (p - 1), 1, d$pr), n)),
    skew = matrix(rexp(n * p) - 1, n), mixed = cbind(matrix(rbinom(n * 5, 1, d$pr), n), matrix(rnorm(n * (p - 5)), n)))
  eta <- drop(X[, seq_along(d$beta), drop = FALSE] %*% d$beta)
  list(X = X, y = rbinom(n, 1, cf_expit(.a_for(eta, d$prev) + eta)))
}
one_bin <- function(r, d, src) {          # src E34: seed 20261132 + r; E34b (rare design only): seed 21261199 + r
  z <- draw(d, if (src == "E34") 20261132 + r else 21261199 + r); cv <- cav_gof_fast(z$X, z$y, d$lambda, G = 10); q3 <- q3_S(cv$xi, z$y, 10)
  data.frame(src = src, rep = r, cell = d$name, q2_chi2 = cv$p_chi2, q2_exact = cv$p_a, q3_chi2G2 = q3[["p"]], tau3 = q3[["tau"]])
}
# ---- Internet Ads (E35) ------------------------------------------------------------------------------------------------
if (PHASE == "ads") {
  SEED_A <- 20261133
  dd <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
  dd <- dd[!is.na(dd$V4), ]; Xa <- as.matrix(dd[, 4:(ncol(dd) - 1)]); ya <- as.numeric(dd[, ncol(dd)] == "ad.")
  Xa <- Xa[, colSums(Xa) >= 5]; na <- nrow(Xa); pa <- ncol(Xa)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED_A)
  cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10); lam_a <- na * cvf$lambda.1se
  ba <- cf_ridge(cbind(1, Xa), ya, lam_a, diag(c(0, rep(1, pa))))
  one_ads <- function(r) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED_A + 1000 + r)
    X1 <- cbind(1, Xa); y <- rbinom(na, 1, cf_expit(drop(X1 %*% ba)))
    bh <- cf_ridge(X1, y, lam_a, diag(c(0, rep(1, pa))))           # as E35 (cold start, tol 1e-10)
    cv <- cav_gof_fast(Xa, y, lam_a, G = 10, beta = bh); q3 <- q3_S(cv$xi, y, 10)
    data.frame(rep = r, q2_chi2 = cv$p_chi2, q2_exact = cv$p_a, q3_chi2G2 = q3[["p"]], tau3 = q3[["tau"]])
  }
}
cl <- makeCluster(W)
clusterExport(cl, c("q3_S", ".a_for", "draw", "one_bin", ls(pattern = "^(cf_|cav_|naive_|statmod_)"),
                    if (PHASE == "ads") c("one_ads", "Xa", "ba", "lam_a", "na", "pa", "SEED_A")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
t0 <- Sys.time()
if (PHASE %in% c("binary", "smoke")) {
  f_out <- file.path(OUT, "E44_q3_binary_perrep.csv")
  jobs <- c(lapply(CELLS, function(d) list(d = d, src = "E34", R = if (PHASE == "smoke") 2 else 500)),
            list(list(d = CELLS[[2]], src = "E34b", R = if (PHASE == "smoke") 2 else 1000)))
  for (jb in jobs) {
    have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$cell == jb$d$name & o$src == jb$src] } else integer(0)
    todo <- setdiff(seq_len(jb$R), have)
    for (batch in split(todo, ceiling(seq_along(todo) / (4 * W)))) {
      o <- do.call(rbind, parLapply(cl, batch, function(r, d, src) tryCatch(one_bin(r, d, src), error = function(e) NULL), d = jb$d, src = jb$src))
      if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out)) }
    cat(sprintf("done %-22s %s (%.1f min)\n", jb$d$name, jb$src, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  if (PHASE == "smoke") print(read.csv(f_out))
}
if (PHASE == "ads") {
  f_out <- file.path(OUT, "E44_ads_placebo_q3_perrep.csv")
  have <- if (file.exists(f_out)) read.csv(f_out)$rep else integer(0); todo <- setdiff(seq_len(R_ADS), have)
  for (batch in split(todo, ceiling(seq_along(todo) / W))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r) tryCatch(one_ads(r), error = function(e) NULL)))
    if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
    cat(sprintf("ads reps %d-%d (%.1f min)\n", min(batch), max(batch), as.numeric(difftime(Sys.time(), t0, units = "mins")))) }
}
stopCluster(cl); cat(sprintf("E44 %s finished (%.1f min)\n", PHASE, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
