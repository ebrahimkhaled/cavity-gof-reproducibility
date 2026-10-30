# E20 -- POWER of the cavity-recalibrated test on E4b's EXACT datasets (SEED 20260908, same generator),
# paired with CALM (decile + adaptive EDGE) and prepivoting from results/E4b2_final_perrep.csv.
# Cells A (n500 p5 lam100), B kappa .10/.25/.40 (lam 416). Departures: null; index Hermite cubic rel .4/.7;
# off-index x1*x2 rel .4/.7.  Also E4b's harder departures at rel 1.0/1.5 are NOT in E4b2, so only the five above.
# Usage: Rscript E20_cavity_power.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20260908                       # E4b stream -> identical datasets
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (tag == "A") 1 else sqrt(1 + .49)) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416), mk("B", 400, 160, 416))
DEPS <- rbind(data.frame(dep = "null", rel = 0), data.frame(dep = "index", rel = c(.4, .7)), data.frame(dep = "offindex", rel = c(.4, .7)))
one <- function(r, d, dep, rel) {
  # data generation copied from E4b line for line
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(dep, null = 0, index = d$sd_eta * rel * (Zs^3 - 3 * Zs) / sqrt(6), offindex = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(d$n, 1, cf_expit(eta))
  c3 <- cav_gof_fast(X, y, d$lambda, G = 10, basis_degree = 3)
  c2 <- cav_gof_fast(X, y, d$lambda, G = 10, basis_degree = 2)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, cav_hl = c3$p_a, cav_chi2 = c3$p_chi2, cav_e3 = c3$p_E, cav_e2 = c2$p_E)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d, dep, rel) tryCatch(one(r, d, dep, rel), error = function(e) NULL), d = d, dep = DEPS$dep[k], rel = DEPS$rel[k]))
  cat(sprintf("done %-7s %-8s %.1f n=%d (%.1f min)\n", d$name, DEPS$dep[k], DEPS$rel[k], nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))))
stopCluster(cl); rownames(res) <- NULL
cmp <- read.csv(file.path(OUT, "E4b2_final_perrep.csv"))[, c("rep", "cell", "dep", "rel", "calm_hl", "calm_edge", "pp_hl", "pp_edge")]
m <- merge(res, cmp, by = c("rep", "cell", "dep", "rel")); stopifnot(nrow(m) == nrow(res))
write.csv(m, file.path(OUT, "E20_cavity_power_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(cav_hl, cav_e3, cav_e2, calm_hl, calm_edge, pp_hl, pp_edge) ~ cell + dep + rel, m, rej)
pd <- function(a, b) aggregate(cbind(d = (m[[a]] < .05) - (m[[b]] < .05)) ~ cell + dep + rel, m, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
d1 <- pd("cav_hl", "pp_hl"); d2 <- pd("cav_e3", "pp_edge")
agg$cavhl_minus_pp <- sprintf("%+.3f (%.3f)", d1$d[, "m"], d1$d[, "se"]); agg$cave3_minus_ppedge <- sprintf("%+.3f (%.3f)", d2$d[, "m"], d2$d[, "se"]); agg$R <- R
write.csv(agg, file.path(OUT, "E20_cavity_power_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E20 | cavity test power vs CALM vs prepivot | E4b data | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg[order(agg$dep, agg$cell, agg$rel), ], digits = 3, row.names = FALSE)
