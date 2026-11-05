# E42 -- fill in the q = 3 statistic where cav_lognormal() stopped with an error in E38 / E39.
# The error is "system is computationally singular" in the information matrix F of the EXACT projection reference
# (it happens when the fitted tau sits at its upper bound 20 or the curve is degenerate, mostly under strong departures).
# The chi2_{G-2} reference does not need F, so the statistic S3 and its chi2_{G-2} P-value can always be computed: this
# script refits exactly as cav_lognormal() does (same start, optimiser and bounds) and stops before F.
# It re-uses the design definitions of E38 and E39 by evaluating those scripts' top sections (everything before the
# cluster is started), regenerates the data of every replicate with a missing q = 3 value, and writes
# results/E42_q3_repair.csv (source, cell, dep, rel, rep, S3, q3_chi2G2, tau3). The summarisers fill the gaps from it.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); RES <- file.path(ROOT, "results")
top_env <- function(file, stop_at) {
  # source() inside these sections loads the engines into the global environment, where cav_gof_fast etc. are found
  ex <- parse(file); txt <- vapply(ex, function(e) paste(deparse(e), collapse = " "), "")
  k <- which(grepl(stop_at, txt, fixed = TRUE))[1]; env <- new.env()
  for (i in seq_len(k - 1)) eval(ex[[i]], env); env
}
cav_lognormal_S <- function(xi, y, G = 10, nodes = 24) {     # the fit of cav_lognormal(), statistic only
  n <- length(y); gh <- statmod_gh(nodes)
  lam <- function(th) { P <- 0; for (k in seq_len(nodes)) P <- P + gh$w[k] * cf_expit(th[1] + th[2] * xi + exp(th[3]) * gh$x[k]); pmin(pmax(P, 1e-9), 1 - 1e-9) }
  nll <- function(th) { m <- lam(th); -sum(y * log(m) + (1 - y) * log(1 - m)) }
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients)
  op <- optim(c(rc, log(0.5)), nll, method = "L-BFGS-B", lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)))
  m <- lam(op$par); v <- m * (1 - m)
  g <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G); idx <- split(seq_len(n), g)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0) / sqrt(vapply(idx, function(I) sum(v[I]), 0))
  list(S = sum(e^2), tau = exp(op$par[3]))
}
out <- list()
# ---- E38 (size) ---------------------------------------------------------------------------------------------------
e38 <- top_env(file.path(ROOT, "sims", "E38_q3_size_all.R"), "f_out <-")
d38 <- read.csv(file.path(RES, "E38_q3_size_all_perrep.csv")); miss <- d38[is.na(d38$q3_chi2G2), ]
for (k in seq_len(nrow(miss))) { d <- Filter(function(z) z$name == miss$cell[k], e38$CELLS)[[1]]; r <- miss$rep[k]
  z <- e38$draw(d, r)
  xi <- if (d$src == "E22") e38$glmnet_xi(z$X, z$y, d$lambda_g, d$alpha) else cav_gof_fast(z$X, z$y, d$lambda, G = d$G)$xi
  f <- cav_lognormal_S(xi, z$y, d$G)
  out[[length(out) + 1]] <- data.frame(source = "E38", cell = d$name, dep = "null", rel = 0, rep = r, S3 = f$S,
                                       q3_chi2G2 = pchisq(f$S, d$G - 2, lower.tail = FALSE), tau3 = f$tau) }
# ---- E39 (power) --------------------------------------------------------------------------------------------------
e39 <- top_env(file.path(ROOT, "sims", "E39_power_q3_boot.R"), "cl <- makeCluster")
d39 <- read.csv(file.path(RES, "E39_power_q3_perrep.csv")); miss <- d39[is.na(d39$q3_chi2G2), ]
for (k in seq_len(nrow(miss))) { d <- Filter(function(z) z$name == miss$cell[k], e39$CELLS)[[1]]; r <- miss$rep[k]
  z <- e39$gen(r, d, miss$dep[k], miss$rel[k]); xi <- cav_gof_fast(z$X, z$y, d$lambda, G = 10)$xi
  f <- cav_lognormal_S(xi, z$y, 10)
  out[[length(out) + 1]] <- data.frame(source = "E39", cell = d$name, dep = miss$dep[k], rel = miss$rel[k], rep = r, S3 = f$S,
                                       q3_chi2G2 = pchisq(f$S, 8, lower.tail = FALSE), tau3 = f$tau) }
res <- do.call(rbind, out); write.csv(res, file.path(RES, "E42_q3_repair.csv"), row.names = FALSE)
cat(sprintf("E42: repaired %d replicates (E38 %d, E39 %d); tau at the upper bound 20 in %d; P < .05 in %d\n", nrow(res),
            sum(res$source == "E38"), sum(res$source == "E39"), sum(res$tau3 > 19.9), sum(res$q3_chi2G2 < .05)))
print(aggregate(cbind(n = 1, rej = q3_chi2G2 < .05, tau_at_20 = tau3 > 19.9) ~ source + cell + dep + rel, res, sum), row.names = FALSE)
