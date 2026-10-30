# E32 -- the TEXTBOOK Hosmer-Lemeshow test (in-sample fitted probabilities, chi2_{G-2}) on E18's exact datasets, so the
# size figure can show all three tests on identical data. Generator copied from E18 (SEED 20261030, same draw()).
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 16
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
src <- readLines(file.path(HERE, "E18_cavity_size.R"))
a <- grep("^SEED <- ", src); b <- grep("^one <- function", src) - 1
eval(parse(text = src[a]))                                   # SEED
eval(parse(text = src[grep("^.a_for <- function", src):b]))  # .a_for, mk, CELLS, draw (verbatim)
one_tb <- function(r, d) {
  z <- draw(d, r); X1 <- cbind(1, z$X); bh <- cf_ridge(X1, z$y, d$lambda, diag(c(0, rep(1, d$p))))
  data.frame(rep = r, cell = d$name, textbook = naive_hl(drop(X1 %*% bh), z$y, d$G))
}
cl <- makeCluster(W); clusterExport(cl, c(ls(), ".a_for")); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), function(r, d) tryCatch(one_tb(r, d), error = function(e) NULL), d = d))))
stopCluster(cl)
e18 <- read.csv(file.path(OUT, "E18_cavity_size_perrep.csv"))[, c("rep", "cell", "p_chi2", "calm_hl")]
m <- merge(res, e18, by = c("rep", "cell")); stopifnot(nrow(m) == nrow(res))
write.csv(m, file.path(OUT, "E32_textbook_on_E18_perrep.csv"), row.names = FALSE)
agg <- aggregate(cbind(cavity = p_chi2 < .05, insample = calm_hl < .05, textbook = textbook < .05) ~ cell, m, mean)
print(agg, digits = 3, row.names = FALSE)
write.csv(agg, file.path(OUT, "E32_textbook_on_E18_summary.csv"), row.names = FALSE)
