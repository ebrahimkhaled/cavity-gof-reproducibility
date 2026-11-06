# E49 -- the linear (one-step) index on the real glaucoma predictors. Trial 3 of the study (see README), real-data check.
# Pre-set reading, before running: compare with the default's leak-free P = .021 and its placebo .041 (E36 B, E37). Standardisation inside each
# fit, exact leave-one-out, penalty 4 x lambda on the linear fit, lambda in 295.3 x (.5, 1, 2). The index is standardised (E47b).
# Earlier warning: a heavy penalty re-ranked 38% of patients and the misfit vanished, so the agreement of the ranking is reported as well.
suppressPackageStartupMessages({ library(parallel); library(TH.data) }); data("GlaucomaM", package = "TH.data")
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
Xraw <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]); yg <- as.numeric(GlaucomaM$Class == "glaucoma"); p <- ncol(Xraw); n <- nrow(Xraw)
q3p <- function(xi, y, G = 10, nodes = 24) { xi <- (xi - mean(xi)) / sd(xi); n <- length(y); gh <- statmod_gh(nodes)
  lam <- function(th) { P <- 0; for (k in seq_len(nodes)) P <- P + gh$w[k] * cf_expit(th[1] + th[2] * xi + exp(th[3]) * gh$x[k]); pmin(pmax(P, 1e-9), 1 - 1e-9) }
  nll <- function(th) { m <- lam(th); -sum(y * log(m) + (1 - y) * log(1 - m)) }
  rc <- suppressWarnings(glm.fit(cbind(1, xi), y, family = binomial())$coefficients)
  op <- optim(c(rc, log(0.5)), nll, method = "L-BFGS-B", lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)))
  m <- lam(op$par); v <- m * (1 - m); g <- pmin(ceiling(rank(xi, ties.method = "first") / (n / G)), G); idx <- split(seq_len(n), g)
  e <- vapply(idx, function(I) sum(y[I] - m[I]), 0) / sqrt(vapply(idx, function(I) sum(v[I]), 0)); pchisq(sum(e^2), length(idx) - 2, lower.tail = FALSE) }
# per-eye standardised designs, computed once: Z[[i]] is the n - 1 training rows and z the held-out row
prep <- lapply(seq_len(n), function(i) { mu <- colMeans(Xraw[-i, ]); s <- apply(Xraw[-i, ], 2, sd); Z <- cbind(1, sweep(sweep(Xraw, 2, mu), 2, s, "/"))
  list(Zi = Z[-i, ], zi = Z[i, ]) })
Aof <- function(lam) lapply(prep, function(q) { D <- diag(c(0, rep(4 * lam, p))); drop(q$zi %*% solve(crossprod(q$Zi) + D, t(q$Zi))) })   # row: index weights on y_(-i)
idx_of <- function(A, y) vapply(seq_len(n), function(i) sum(A[[i]] * y[-i]), 0)
LAMS <- 295.3 * c(.5, 1, 2); AS <- lapply(LAMS, Aof)
real <- sapply(AS, function(A) q3p(idx_of(A, yg), yg)); names(real) <- LAMS
cat("real glaucoma, linear index, P at lambda", LAMS, ":", round(real, 4), "\n")
# default index on the same data, for rank agreement (one-step, leak-free version is in E36; ranking compared with the exact refit index)
Xs <- scale(Xraw); bfit <- cav_pieces(Xs, yg, 295.3); L <- idx_of(AS[[2]], yg)
d0 <- pmin(ceiling(rank(bfit$xi, ties.method = "first") / (n / 10)), 10); d1 <- pmin(ceiling(rank(L, ties.method = "first") / (n / 10)), 10)
cat("corr with logistic cavity index", round(cor(bfit$xi, L), 4), " Spearman", round(cor(bfit$xi, L, method = "spearman"), 4), " same decile", round(mean(d0 == d1), 3), "\n")
for (G in c(5)) cat("G=5, lambda 295.3:", round(q3p(L, yg, G), 4), "\n")
b_true <- cf_ridge(cbind(1, Xs), yg, 295.3, diag(c(0, rep(1, p)))); pr <- cf_expit(drop(cbind(1, Xs) %*% b_true))
one <- function(r) { RNGkind("L'Ecuyer-CMRG"); set.seed(20261137 + r); y <- rbinom(n, 1, pr)     # E37's y
  c(r, tryCatch(q3p(idx_of(AS[[2]], y), y), error = function(e) NA_real_)) }
cl <- makeCluster(6); clusterExport(cl, c("one", "q3p", "idx_of", "AS", "pr", "n", "cf_expit", "statmod_gh"))
o <- do.call(rbind, parLapply(cl, 1:1000, one)); stopCluster(cl)
write.csv(data.frame(rep = o[, 1], p = o[, 2]), file.path(OUT, "E49_placebo_perrep.csv"), row.names = FALSE)
e37 <- read.csv(file.path(OUT, "E37_glaucoma_placebo_noleak_perrep.csv")); m <- merge(data.frame(rep = o[, 1], lin = o[, 2]), e37, by = "rep")
cat("placebo linear index", mean(o[, 2] < .05, na.rm = TRUE), " at 1%", mean(o[, 2] < .01, na.rm = TRUE), " default d=3 same y", mean(m$new_q3_chi2G2 < .05, na.rm = TRUE), " NA", sum(is.na(o[, 2])), "\n")
