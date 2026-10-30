# E30 -- power against STRONGER departures (E20's were too weak to separate the tests outside fixed p).
# Cavity (decile, smooth degree 2 and 3) vs CALM (decile, adaptive EDGE) on identical data. Departures (E4b dial REL):
# Hermite cubic and quadratic along the true index, interaction x1*x2 orthogonal to it; REL 1.0 and 1.5 (kept < 0.8165 for
# monotonicity is NOT required here since we only need a correctly generated binary outcome).
# Usage: Rscript E30_power_strong.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 12
SEED <- 20261110
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (tag == "A") 1 else sqrt(1 + .49)) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416))
DEPS <- rbind(data.frame(dep = "cubic", rel = c(1, 1.5)), data.frame(dep = "quadratic", rel = c(1, 1.5)), data.frame(dep = "interaction", rel = c(1, 1.5)))
one <- function(r, d, dep, rel) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(dep, cubic = d$sd_eta * rel * (Zs^3 - 3 * Zs) / sqrt(6), quadratic = d$sd_eta * rel * (Zs^2 - 1) / sqrt(2),
                      interaction = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(d$n, 1, cf_expit(eta))
  c3 <- cav_gof_fast(X, y, d$lambda, G = 10, basis_degree = 3); c2 <- cav_gof_fast(X, y, d$lambda, G = 10, basis_degree = 2)
  cf <- tryCatch(calm.gof(X, y, d$lambda), error = function(e) NULL)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, cav_hl = c3$p_a, cav_e3 = c3$p_E, cav_e2 = c2$p_E,
             calm_hl = if (is.null(cf)) NA else cf$SC.HL$p.value, calm_edge = if (is.null(cf)) NA else cf$SC.EDGE.adaptive$p.value)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d, dep, rel) tryCatch(one(r, d, dep, rel), error = function(e) NULL), d = d, dep = DEPS$dep[k], rel = DEPS$rel[k]))
  cat(sprintf("done %-7s %-12s %.1f n=%d (%.1f min)\n", d$name, DEPS$dep[k], DEPS$rel[k], nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))))
stopCluster(cl); write.csv(res, file.path(OUT, "E30_power_strong_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(cav_hl, cav_e3, cav_e2, calm_hl, calm_edge) ~ cell + dep + rel, res, rej, na.action = na.pass)
d <- aggregate(cbind(dd = (res$cav_hl < .05) - (res$calm_hl < .05)) ~ cell + dep + rel, res, function(v) c(m = mean(v, na.rm = TRUE), se = sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))))
agg$cavhl_minus_calmhl <- sprintf("%+.3f (%.3f)", d$dd[, "m"], d$dd[, "se"]); agg$R <- R
write.csv(agg, file.path(OUT, "E30_power_strong_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E30 | power against stronger departures | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg[order(agg$cell, agg$dep, agg$rel), ], digits = 3, row.names = FALSE)
