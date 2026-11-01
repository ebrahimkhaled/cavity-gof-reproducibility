# E33 -- the cavity test against two further comparators, on identical data:
#   (1) a K-fold cross-fitted recalibrated Hosmer-Lemeshow test (K = 5, 10): each fold's index comes from the fit on
#       the other folds (exact out-of-fold, no one-step approximation); then the same recalibration, grouping and
#       chi-squared(G-2) reference as the cavity test, on the whole sample;
#   (2) the residual-prediction test of Jankova, Shah, Buhlmann and Samworth (2020), GRPtests::GRPtest, package defaults.
# Size in designs A and B (index sd 2.31, as E18) and at p = 2n; power against the quadratic and the interaction of size 1.5
# on E30's generator (companion-paper scaling, ||beta0|| = 2.31). Every cell is appended to disk as it finishes; a rerun
# skips cells already written.
# Usage: Rscript E33_comparators.R [R] [workers] [R_grp]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 20
R_GRP <- if (length(args) >= 3) as.integer(args[3]) else 200
# GRPtest grows random forests (about two minutes per call on one thread at p = 100), so it runs on these cells only
GRP_CELLS <- c("B_p100", "B_p800", "B_p100_quad", "B_p100_int")
SEED <- 20261131
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

mk <- function(name, n, p, lambda, scaling = c("sd", "norm"), indep = FALSE, dep = "null", rel = 0) {
  scaling <- match.arg(scaling)
  if (indep) { Sig <- diag(p); b0 <- c(.5, -.4, .3, -.3, .2) } else {
    Sig <- 0.7^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p)
    b0 <- if (scaling == "sd") b0 * 2.31 / sqrt(drop(t(b0) %*% Sig %*% b0)) else b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), dep = dep, rel = rel,
       sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (indep) 1 else sqrt(1 + .49))
}
CELLS <- list(
  mk("A_p5", 500, 5, 100, indep = TRUE), mk("B_p40", 400, 40, 416), mk("B_p100", 400, 100, 416),
  mk("B_p160", 400, 160, 416), mk("B_p800", 400, 800, 416),
  mk("B_p40_quad", 400, 40, 416, "norm", dep = "quadratic", rel = 1.5), mk("B_p100_quad", 400, 100, 416, "norm", dep = "quadratic", rel = 1.5),
  mk("B_p40_int", 400, 40, 416, "norm", dep = "interaction", rel = 1.5), mk("B_p100_int", 400, 100, 416, "norm", dep = "interaction", rel = 1.5))

kfold_index <- function(X, y, lambda, K, seed) {
  n <- length(y); set.seed(seed); fold <- sample(rep_len(seq_len(K), n)); xi <- numeric(n)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, ncol(X))))
  for (k in seq_len(K)) { tr <- fold != k
    b <- cf_ridge(X1[tr, , drop = FALSE], y[tr], lambda * sum(tr) / n, Dpen)   # same lambda / n as the full fit
    xi[!tr] <- drop(X1[!tr, , drop = FALSE] %*% b) }
  xi
}
one <- function(r, d, with_grp) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0)
  if (d$dep != "null") { Zs <- eta / d$sd_eta
    eta <- eta + switch(d$dep, quadratic = d$sd_eta * d$rel * (Zs^2 - 1) / sqrt(2),
                        interaction = (d$rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2]) }
  y <- rbinom(d$n, 1, cf_expit(eta))
  cv <- cav_gof_fast(X, y, d$lambda, G = 10)
  k5 <- cav_from_index(kfold_index(X, y, d$lambda, 5, SEED + 7 * r), y)$p_chi2
  k10 <- cav_from_index(kfold_index(X, y, d$lambda, 10, SEED + 11 * r), y)$p_chi2
  X1 <- cbind(1, X); bh <- cf_ridge(X1, y, d$lambda, diag(c(0, rep(1, d$p))))
  grp <- if (with_grp) tryCatch(as.numeric(suppressWarnings(suppressMessages(GRPtests::GRPtest(X, y, fam = "binomial")))),
                                error = function(e) NA_real_) else NA_real_
  data.frame(rep = r, cell = d$name, cavity = cv$p_chi2, kfold5 = k5, kfold10 = k10, grp = grp,
             textbook = naive_hl(drop(X1 %*% bh), y, 10))
}

f_out <- file.path(OUT, "E33_comparators_perrep.csv")
done <- if (file.exists(f_out)) unique(read.csv(f_out)$cell) else character(0)
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls())
# ranger (inside GRPtest) uses every core by default; with parallel workers that oversubscribes and crashes
invisible(clusterEvalQ(cl, { suppressPackageStartupMessages({ library(CompQuadForm); library(GRPtests) }); options(ranger.num.threads = 1) }))
for (d in CELLS) {
  if (d$name %in% done) next
  rg <- if (d$name %in% GRP_CELLS) R_GRP else 0
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d, rg) tryCatch(one(r, d, r <= rg), error = function(e) NULL),
                                d = d, rg = rg))
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("done %-12s n=%d  (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out); rej <- function(v) if (all(is.na(v))) NA else mean(v < .05, na.rm = TRUE)
agg <- do.call(rbind, lapply(split(res, res$cell), function(s) data.frame(cell = s$cell[1], R = nrow(s), R_grp = sum(!is.na(s$grp)),
  cavity = rej(s$cavity), kfold5 = rej(s$kfold5), kfold10 = rej(s$kfold10), grp = rej(s$grp), textbook = rej(s$textbook))))
agg <- agg[match(vapply(CELLS, `[[`, "", "name"), agg$cell), ]
write.csv(agg, file.path(OUT, "E33_comparators_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E33 | cavity vs K-fold recalibration vs GRPtest | %.1f min ===\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
