# E39 -- power of the q = 3 cavity test (chi2_{G-2}) next to q = 2, and the parametric bootstrap (prepivoting, the
# companion paper's comparator, ebrahim.gof::shrink.gof, B = 399 refits) on E30's IDENTICAL data sets
# (referee items M2 / S2 and M7(a) / S4(a)). E30 = Figure 4 and Supplement tab:powerstrong.
# Data regenerated with E30's generator: RNGkind reset to Mersenne-Twister / Inversion / Rejection (the default E30 ran
# with), set.seed(20261110 + r). The q = 2 exact-reference P-value is checked against E30's stored cav_hl in E39_summarise.
# Phase "q3": all 18 (cell, departure, size) configurations, 500 replicates -> results/E39_power_q3_perrep.csv
# Phase "boot": size 1.5 only (9 configurations), 500 replicates, bootstrap seed 20261110 + 1e6 + r (as E4b)
#               -> results/E39_power_boot_perrep.csv
# Every batch is appended; a rerun skips replicates already written.
# Usage: Rscript E39_power_q3_boot.R [workers] [phase: q3 | boot | all | smoke] [R]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) as.integer(args[1]) else 10; PHASE <- if (length(args) >= 2) args[2] else "all"
R <- if (length(args) >= 3) as.integer(args[3]) else 500; SMOKE <- PHASE == "smoke"; if (SMOKE) R <- 2
SEED <- 20261110; B_BOOT <- 399
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- if (SMOKE) tempdir() else file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (tag == "A") 1 else sqrt(1 + .49)) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416))
DEPS <- rbind(data.frame(dep = "cubic", rel = c(1, 1.5)), data.frame(dep = "quadratic", rel = c(1, 1.5)), data.frame(dep = "interaction", rel = c(1, 1.5)))
gen <- function(r, d, dep, rel) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(SEED + r)
  X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(dep, cubic = d$sd_eta * rel * (Zs^3 - 3 * Zs) / sqrt(6), quadratic = d$sd_eta * rel * (Zs^2 - 1) / sqrt(2),
                      interaction = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  list(X = X, y = rbinom(d$n, 1, cf_expit(eta)))
}
one_q3 <- function(r, d, dep, rel) {
  z <- gen(r, d, dep, rel); c2 <- cav_gof_fast(z$X, z$y, d$lambda, G = 10)
  c3 <- tryCatch(cav_lognormal(c2$xi, z$y, G = 10), error = function(e) NULL)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, q2_exact = c2$p_a, q2_chi2 = c2$p_chi2,
             q3_chi2G2 = if (is.null(c3)) NA else pchisq(c3$S, 8, lower.tail = FALSE), q3_exact = if (is.null(c3)) NA else c3$p_a,
             tau3 = if (is.null(c3)) NA else c3$tau)
}
one_boot <- function(r, d, dep, rel) {
  z <- gen(r, d, dep, rel)
  pp <- ebrahim.gof::shrink.gof(z$X, z$y, lambda = d$lambda, G = 10, basis = c("edge", "decile"), B = B_BOOT, seed = SEED + 1e6 + r)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, pp_hl = pp$SC.HL$p.value, pp_edge = pp$SC.EDGE$p.value)
}
cl <- makeCluster(W); clusterExport(cl, c("gen", "one_q3", "one_boot", "SEED", "B_BOOT", ls(pattern = "^(cf_|cav_|naive_|statmod_)")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({ library(CompQuadForm); library(ebrahim.gof) })))
t0 <- Sys.time()
run <- function(fun, f_out, rels) for (d in CELLS) for (k in seq_len(nrow(DEPS))) {
  dep <- DEPS$dep[k]; rel <- DEPS$rel[k]; if (!(rel %in% rels)) next
  have <- if (file.exists(f_out)) { o <- read.csv(f_out); o$rep[o$cell == d$name & o$dep == dep & o$rel == rel] } else integer(0)
  todo <- setdiff(seq_len(R), have)
  for (batch in split(todo, ceiling(seq_along(todo) / (2 * W)))) {
    o <- do.call(rbind, parLapply(cl, batch, function(r, d, dep, rel, ff) tryCatch(ff(r, d, dep, rel), error = function(e) NULL),
                                  d = d, dep = dep, rel = rel, ff = fun))
    if (!is.null(o)) write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  }
  cat(sprintf("done %-7s %-12s %.1f (%.1f min)\n", d$name, dep, rel, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (PHASE %in% c("q3", "all", "smoke")) run(one_q3, file.path(OUT, "E39_power_q3_perrep.csv"), c(1, 1.5))
if (PHASE %in% c("boot", "all", "smoke")) run(one_boot, file.path(OUT, "E39_power_boot_perrep.csv"), 1.5)
if (SMOKE) { print(read.csv(file.path(OUT, "E39_power_q3_perrep.csv"))); print(read.csv(file.path(OUT, "E39_power_boot_perrep.csv"))) }
stopCluster(cl); cat(sprintf("E39 finished (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
