# E34b -- repeat of E34's one borderline cell (one rare, strong binary predictor; size 0.072 from 500 replications)
# with 1000 replications on seeds disjoint from E34, to separate a real excess from chance.
# Usage: Rscript E34b_confirm_rare_binary.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 21261199   # offset 10^6 from E34 (seeds 20261132 + r), so no data set is shared with it
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
.a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
one <- function(r) {
  set.seed(SEED + r); n <- 400; p <- 100; beta <- c(2.5, 0.6, -0.6, 0.6, -0.6)
  X <- cbind(rbinom(n, 1, 0.1), matrix(rbinom(n * (p - 1), 1, 0.3), n))
  eta <- drop(X[, 1:5] %*% beta); y <- rbinom(n, 1, cf_expit(.a_for(eta, 0.2) + eta))
  cv <- cav_gof_fast(X, y, 100, G = 10)
  data.frame(rep = r, cavity = cv$p_chi2, cavity_exact = cv$p_a)
}
cl <- makeCluster(W); clusterExport(cl, c(ls(), ".a_for")); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r), error = function(e) NULL))); stopCluster(cl)
write.csv(res, file.path(OUT, "E34b_confirm_rare_binary_perrep.csv"), row.names = FALSE)
m <- mean(res$cavity < .05); cat(sprintf("E34b | rare strong binary predictor | R = %d | size %.3f (SE %.3f); exact projection %.3f\n",
                                         nrow(res), m, sqrt(m * (1 - m) / nrow(res)), mean(res$cavity_exact < .05)))
