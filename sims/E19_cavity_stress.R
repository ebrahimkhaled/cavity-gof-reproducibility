# E19 -- STRESS TEST of the cavity-recalibrated test (cav_gof_fast) under the null, in designs chosen to break it:
#   larger n at fixed p/n (the logistic-normal recalibration curve is not exactly logistic: a drift would show here),
#   strong / weak signal, light / heavy penalty, p/n = 0.5, binary and t3 designs, G = 5 / 20,
#   prevalence 0.05 with a strong signal, and a large-n fixed-p cell.
# Decile statistic (chi2_{G-2} and the exact projection reference) and the smooth-basis statistic (degree 2 and 3).
# Usage: Rscript E19_cavity_stress.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261031
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, law = "gaussian", prev = 0.5, G = 10, indep = FALSE, R = NULL) {
  if (indep) { Sig <- diag(p); b0 <- rep_len(c(.5, -.4, .3, -.3, .2), p) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), law = law, prev = prev, G = G, R = R) }
CELLS <- list(
  mk("n1600_k25", 1600, 400, 1664, R = 300), mk("n3200_k25", 3200, 800, 3328, R = 200),
  mk("strong4_k25", 400, 100, 416, norm = 4), mk("weak1_k25", 400, 100, 416, norm = 1),
  mk("lam10_k25", 400, 100, 10), mk("lam5000_k25", 400, 100, 5000), mk("k50", 400, 200, 416),
  mk("binary_k25", 400, 100, 416, law = "binary"), mk("t3_k25", 400, 100, 416, law = "t3"),
  mk("G5_k25", 400, 100, 416, G = 5), mk("G20_k25", 400, 100, 416, G = 20),
  mk("prev05_strong3_k25", 400, 100, 416, norm = 3, prev = .05), mk("A_n5000", 5000, 5, 1000, indep = TRUE, R = 300))
draw <- function(d, r) {
  set.seed(SEED + r)
  Zm <- switch(d$law, gaussian = matrix(rnorm(d$n * d$p), d$n), t3 = matrix(rt(d$n * d$p, 3) / sqrt(3), d$n),
               binary = matrix((rbinom(d$n * d$p, 1, .3) - .3) / sqrt(.21), d$n))
  X <- Zm %*% d$Sc; eta <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else a_for(eta, d$prev)
  list(X = X, y = rbinom(d$n, 1, cf_expit(a + eta)))
}
one <- function(r, d) {
  z <- draw(d, r)
  c3 <- cav_gof_fast(z$X, z$y, d$lambda, G = d$G, basis_degree = 3)
  c2 <- cav_gof_fast(z$X, z$y, d$lambda, G = d$G, basis_degree = 2)
  data.frame(rep = r, cell = d$name, p_chi2 = c3$p_chi2, p_a = c3$p_a, p_E3 = c3$p_E, p_E2 = c2$p_E,
             S = c3$S, tr_a = c3$tr_a, b = c3$b, Hbar = c3$Hbar, t(setNames(c3$e, paste0("e", seq_along(c3$e)))))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
f_out <- file.path(OUT, "E19_cavity_stress_perrep.csv"); if (file.exists(f_out)) file.remove(f_out)
for (d in CELLS) {
  Rd <- if (is.null(d$R)) R else min(R, d$R)
  o <- do.call(rbind, parLapply(cl, seq_len(Rd), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
  if (is.null(o)) { cat("ALL FAILED", d$name, "\n"); next }
  o <- o[, c("rep", "cell", "p_chi2", "p_a", "p_E3", "p_E2", "S", "tr_a", "b", "Hbar")]   # G differs across cells: keep fixed columns on disk
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("done %-20s n=%d  (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(p_chi2, p_a, p_E3, p_E2) ~ cell, res, rej, na.action = na.pass)
agg$meanS <- aggregate(S ~ cell, res, mean)$S; agg$tr_a <- aggregate(tr_a ~ cell, res, mean)$tr_a
agg$b <- aggregate(b ~ cell, res, mean)$b; agg$Hbar <- aggregate(Hbar ~ cell, res, mean)$Hbar; agg$R <- aggregate(rep ~ cell, res, length)$rep
agg <- agg[match(sapply(CELLS, `[[`, "name"), agg$cell), ]
write.csv(agg, file.path(OUT, "E19_cavity_stress_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E19 | cavity test stress (null) | %.1f min ===\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
