# E18 -- size of the CAVITY-RECALIBRATED grouped test (R/cavity_engine.R) against CALM's decile statistic,
# on the same datasets, in the registered cells AND in every cell where CALM-decile failed:
#   E10 held-out failures (kappa .40 lambda 50; AR 0.3; n = 800; t3 at kappa .10) and E13 low prevalence.
# Per replicate: p-values under four references, the statistic, the trace of each reference and the
# standardised group residuals (to measure the null mean directly).
# Usage: Rscript E18_cavity_size.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 300; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261030
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

.a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-12, 12))$root
mk <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, law = "gaussian", prev = 0.5, G = 10, indep = FALSE) {
  if (indep) { Sig <- diag(p); b0 <- c(.5, -.4, .3, -.3, .2) } else {
    Sig <- ar^abs(outer(1:p, 1:p, "-")); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)) }
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), law = law, prev = prev, G = G) }
CELLS <- list(
  mk("A_p5", 500, 5, 100, indep = TRUE), mk("B_p40", 400, 40, 416), mk("B_p100", 400, 100, 416), mk("B_p160", 400, 160, 416),
  mk("H3_k40_lam50", 400, 160, 50), mk("H5_k25_ar03", 400, 100, 416, ar = 0.3), mk("H8_n800_k25", 800, 200, 832),
  mk("k10_t3", 400, 40, 416, law = "t3"), mk("prev15_k25", 400, 100, 416, prev = .15), mk("prev08_k25", 400, 100, 416, prev = .08),
  mk("prev08_A", 500, 5, 100, indep = TRUE, prev = .08))
draw <- function(d, r) {
  set.seed(SEED + r)
  Zm <- switch(d$law, gaussian = matrix(rnorm(d$n * d$p), d$n), t3 = matrix(rt(d$n * d$p, 3) / sqrt(3), d$n))
  X <- Zm %*% d$Sc; eta <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else .a_for(eta, d$prev)
  list(X = X, y = rbinom(d$n, 1, cf_expit(a + eta)))
}
one <- function(r, d) {
  z <- draw(d, r)
  cv <- cav_gof(z$X, z$y, d$lambda, G = d$G)
  cf <- tryCatch(calm.gof(z$X, z$y, d$lambda, G = d$G, basis = "decile")$SC.HL$p.value, error = function(e) NA_real_)
  data.frame(rep = r, cell = d$name, p_I = cv$p_I, p_chi2 = cv$p_chi2, p_a = cv$p_a, p_b = cv$p_b, calm_hl = cf,
             S = cv$S, tr_a = cv$tr_a, tr_b = cv$tr_b, b = cv$b, Hbar = cv$Hbar, t(setNames(cv$e, paste0("e", seq_along(cv$e)))))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, c(ls(), ".a_for")); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
f_out <- file.path(OUT, "E18_cavity_size_perrep.csv"); if (file.exists(f_out)) file.remove(f_out)
for (d in CELLS) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("done %-14s n=%d  (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(p_I, p_chi2, p_a, p_b, calm_hl) ~ cell, res, rej, na.action = na.pass)
ecols <- grep("^e[0-9]+$", names(res))
agg$meanS <- aggregate(S ~ cell, res, mean)$S; agg$tr_a <- aggregate(tr_a ~ cell, res, mean)$tr_a; agg$tr_b <- aggregate(tr_b ~ cell, res, mean)$tr_b
agg$Ee2 <- sapply(agg$cell, function(cc) sum(colMeans(res[res$cell == cc, ecols], na.rm = TRUE)^2))
agg$b <- aggregate(b ~ cell, res, mean)$b; agg$R <- aggregate(rep ~ cell, res, length)$rep
agg <- agg[match(sapply(CELLS, `[[`, "name"), agg$cell), ]
write.csv(agg, file.path(OUT, "E18_cavity_size_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E18 | cavity-recalibrated test vs CALM decile | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
