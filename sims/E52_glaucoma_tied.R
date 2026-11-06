# E52 -- tied recalibration on the real glaucoma data. See README, real-data check. Descriptive: no pass rule, and it decides nothing.
# Same indices as E49 (standardisation inside each fit, penalty 4 x lambda on the linear fit). rho^2 = corr(y, L)^2 / adjusted R^2 as in E50.
suppressPackageStartupMessages({ library(TH.data) }); data("GlaucomaM", package = "TH.data")
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
s50 <- readLines(file.path(ROOT, "sims", "E50_tau_separate.R")); eval(parse(text = s50[grep("^cells_e", s50):(grep("^lin_b", s50) - 1)]))
std <- function(v) (v - mean(v)) / sd(v)
Xraw <- as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]); y <- as.numeric(GlaucomaM$Class == "glaucoma"); p <- ncol(Xraw); n <- nrow(Xraw)
prep <- lapply(seq_len(n), function(i) { mu <- colMeans(Xraw[-i, ]); s <- apply(Xraw[-i, ], 2, sd); Z <- cbind(1, sweep(sweep(Xraw, 2, mu), 2, s, "/")); list(Zi = Z[-i, ], zi = Z[i, ]) })
idx <- function(pen) vapply(seq_len(n), function(i) { q <- prep[[i]]; D <- diag(c(0, rep(pen, p))); sum(drop(q$zi %*% solve(crossprod(q$Zi) + D, t(q$Zi))) * y[-i]) }, 0)
a <- adjR2(scale(Xraw), y); cat("n", n, "p", p, "adjusted R2", round(a, 3), "\n")
Xs <- scale(Xraw); full <- cav_pieces(Xs, y, 295.3)$xi
L <- c(lapply(295.3 * c(.5, 1, 2), function(l) idx(4 * l)), list(idx(0)), list(full)); names(L) <- c("I1_lam147.65", "I1_lam295.3", "I1_lam590.6", "I2_OLS", "full_295.3")
for (nm in names(L)) { xi <- std(L[[nm]]); r <- rho_R2(L[[nm]], y, a); k <- kap(r); o <- d3(xi, y); t0 <- tied(xi, y, 0); q <- tied(xi, y, k)
  cat(sprintf("%-14s rho2 %.3f kappa %.3f | P d3 %.4f (tau %.3f) | T0 %.4f | tied %.4f (tau %.3f)\n", nm, r, k, o[1], o[2], t0[1], q[1], q[2])) }
