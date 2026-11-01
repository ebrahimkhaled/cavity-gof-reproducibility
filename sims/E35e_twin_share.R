# E35e -- the twin-share diagnostic of Supplementary Section S9.4 and the simulation-calibrated P value of the
# Internet-Ads data.  Observation j is a twin of k when x_k' M^{-1} x_j / sqrt(q_k q_j) > 0.999 (cut-off 0.9 is
# also shown, to make the point that it flags low-dimensional designs where the test holds).
# Usage: Rscript E35e_twin_share.R
suppressPackageStartupMessages({ library(glmnet); library(TH.data) })
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R")); source(file.path(ROOT, "R", "cluster_correction.R"))

share <- function(X, lambda, y) {
  b <- cf_ridge(cbind(1, X), y, lambda, diag(c(0, rep(1, ncol(X)))))
  c(cut_0.9 = cav_twin_share(cav_clusters(X, b, lambda, 0.9)), cut_0.999 = cav_twin_share(cav_clusters(X, b, lambda, 0.999)))
}

# one data set of each simulated design family, n = 400
set.seed(11); n <- 400
designs <- list(
  "A, p = 5"            = list(matrix(rnorm(n * 5), n), 1),
  "B, p = 40"           = list(matrix(rnorm(n * 40), n), 100),
  "B, p = 100"          = list(matrix(rnorm(n * 100), n), 100),
  "B, p = 800"          = list(matrix(rnorm(n * 800), n), 100),
  "t3, p = 100"         = list(matrix(rt(n * 100, 3), n), 100),
  "binary, p = 100"     = list(matrix(rbinom(n * 100, 1, 0.3), n), 100),
  "binary, p = 800"     = list(matrix(rbinom(n * 800, 1, 0.3), n), 100),
  "rare binary, p = 100" = list(cbind(rbinom(n, 1, 0.1), matrix(rbinom(n * 99, 1, 0.3), n)), 100),
  "exponential, p = 100" = list(matrix(rexp(n * 100) - 1, n), 100))
tab <- do.call(rbind, lapply(names(designs), function(nm) {
  X <- designs[[nm]][[1]]; p <- ncol(X)
  y <- rbinom(n, 1, cf_expit(drop(cbind(1, X) %*% c(0, rnorm(p) * 2.31 / sqrt(p)))))
  data.frame(design = nm, t(share(X, designs[[nm]][[2]], y)))
}))

data("GlaucomaM", package = "TH.data")
Xg <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); yg <- as.numeric(GlaucomaM$Class == "glaucoma")
tab <- rbind(tab, data.frame(design = "glaucoma (real)", t(share(Xg, 295.3, yg))))

d <- read.csv(file.path(ROOT, "data_raw", "internet_ads", "ad.data"), header = FALSE, na.strings = "?", strip.white = TRUE)
d <- d[!is.na(d$V4), ]
Xa <- as.matrix(d[, 4:(ncol(d) - 1)]); ya <- as.numeric(d[, ncol(d)] == "ad.")
Xa <- Xa[, colSums(Xa) >= 5]
set.seed(20261133); cvf <- cv.glmnet(Xa, ya, family = "binomial", alpha = 0, standardize = FALSE, nfolds = 10)
tab <- rbind(tab, data.frame(design = "Internet Ads (real)", t(share(Xa, nrow(Xa) * cvf$lambda.1se, ya))))
rows <- apply(Xa, 1, paste, collapse = ""); blocks <- table(rows)

# simulation-calibrated P: the real cavity P value (0.0090, E35) against the 100 placebo P values of E35
pl <- read.csv(file.path(OUT, "E35_placebo_perrep.csv")); pl <- pl$cavity[pl$data == "ads"]
k <- sum(pl <= 0.0090)

sink(file.path(OUT, "E35e_twin_share.txt"))
cat("=== E35e | share of observations with a twin (leave-one-out cosine above the cut-off) ===\n")
print(tab, digits = 3, row.names = FALSE)
cat(sprintf("\nInternet Ads: %d distinct rows among %d images; %.0f%% of images share their row; largest block %d\n",
            length(blocks), nrow(Xa), 100 * mean(blocks[rows] > 1), max(blocks)))
cat(sprintf("simulation-calibrated P for the real Internet-Ads data: %d of %d placebo statistics at least as large, P = %.2f\n",
            k, length(pl), (k + 1) / (length(pl) + 1)))
sink()
cat(readLines(file.path(OUT, "E35e_twin_share.txt")), sep = "\n")
