# E34 -- size of the cavity test where the Gaussian calibration law is least exact: a few strong binary or skewed
# predictors (sparse signal), the typical clinical prediction model. Ridge fits at fixed penalty; G = 10; the cavity
# test with the two-parameter recalibration and the chi-squared(G-2) reference (the paper's rule for n < 5000), its
# exact-projection version, and the textbook Hosmer-Lemeshow test. Each cell is appended to disk; a rerun skips cells done.
# Usage: Rscript E34_sparse_binary_stress.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20261132
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))

.a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
# law: "binary" 0/1 predictors with success probability `pr`; "skew" centred exponential; "mixed" binary signal columns
# plus Gaussian noise columns. beta: the non-zero coefficients, placed on the first columns.
mk <- function(name, n, p, lambda, law, beta, prev, pr = 0.3) list(name = name, n = n, p = p, lambda = lambda, law = law,
                                                                   beta = beta, prev = prev, pr = pr)
CELLS <- list(
  mk("bin5_p100", 400, 100, 100, "binary", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("bin_rare_strong_p100", 400, 100, 100, "binary_rare", c(2.5, 0.6, -0.6, 0.6, -0.6), prev = 0.2),
  mk("bin5_p200", 400, 200, 100, "binary", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("bin5_p800", 400, 800, 100, "binary", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("skew5_p100", 400, 100, 100, "skew", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3),
  mk("mixed5_p100", 400, 100, 100, "mixed", c(1.2, -1.0, 0.9, -0.8, 0.7), prev = 0.3))

draw <- function(d, r) {
  set.seed(SEED + r); n <- d$n; p <- d$p
  X <- switch(d$law,
    binary = matrix(rbinom(n * p, 1, d$pr), n),
    binary_rare = cbind(rbinom(n, 1, 0.1), matrix(rbinom(n * (p - 1), 1, d$pr), n)),   # one rare, strong predictor
    skew = matrix(rexp(n * p) - 1, n),
    mixed = cbind(matrix(rbinom(n * 5, 1, d$pr), n), matrix(rnorm(n * (p - 5)), n)))
  eta <- drop(X[, seq_along(d$beta), drop = FALSE] %*% d$beta)
  list(X = X, y = rbinom(n, 1, cf_expit(.a_for(eta, d$prev) + eta)))
}
one <- function(r, d) {
  z <- draw(d, r); cv <- cav_gof_fast(z$X, z$y, d$lambda, G = 10)
  X1 <- cbind(1, z$X); bh <- cf_ridge(X1, z$y, d$lambda, diag(c(0, rep(1, d$p))))
  data.frame(rep = r, cell = d$name, cavity = cv$p_chi2, cavity_exact = cv$p_a, textbook = naive_hl(drop(X1 %*% bh), z$y, 10),
             b = cv$b, Hbar = cv$Hbar)
}

f_out <- file.path(OUT, "E34_sparse_binary_perrep.csv")
done <- if (file.exists(f_out)) unique(read.csv(f_out)$cell) else character(0)
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, c(ls(), ".a_for"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
for (d in CELLS) {
  if (d$name %in% done) next
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one(r, d), error = function(e) NULL), d = d))
  write.table(o, f_out, sep = ",", row.names = FALSE, col.names = !file.exists(f_out), append = file.exists(f_out))
  cat(sprintf("done %-22s n=%d  (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
res <- read.csv(f_out); rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- do.call(rbind, lapply(split(res, res$cell), function(s) data.frame(cell = s$cell[1], R = nrow(s),
  cavity = rej(s$cavity), cavity_exact = rej(s$cavity_exact), textbook = rej(s$textbook), slope = mean(s$b), Hbar = mean(s$Hbar))))
agg <- agg[match(vapply(CELLS, `[[`, "", "name"), agg$cell), ]
write.csv(agg, file.path(OUT, "E34_sparse_binary_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E34 | sparse binary / skewed predictors (null) | %.1f min ===\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg, digits = 3, row.names = FALSE)
