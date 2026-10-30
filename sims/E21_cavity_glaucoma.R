# E21 -- the cavity-recalibrated test on the glaucoma data (TH.data::GlaucomaM, n = 196, p = 62, kappa = 0.32),
# at the four penalties of the companion paper (lambda_min, lambda_1se / 2, lambda_1se, 2 lambda_1se; lambda_1se = 295.3 on the
# theory scale), G = 10 and 5/20. Compared with CALM (decile, adaptive EDGE) on the same fit.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
source(file.path(ROOT, "R", "closedform_engine.R")); source(file.path(ROOT, "R", "cavity_engine.R"))
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
Xg <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); yg <- as.numeric(GlaucomaM$Class == "glaucoma")
n <- nrow(Xg)
for (lam in c(295.3 / 2, 295.3, 2 * 295.3)) for (G in c(5, 10, 20)) {
  c3 <- cav_gof_fast(Xg, yg, lam, G = G, basis_degree = 3); c2 <- cav_gof_fast(Xg, yg, lam, G = G, basis_degree = 2)
  cf <- calm.gof(Xg, yg, lam, G = G)
  cat(sprintf("lambda %6.1f G %2d | cavity: HL p = %.4f (chi2 %.4f)  E3 p = %.4f  E2 p = %.4f  slope b = %.2f  Hbar = %.3f | CALM: HL p = %.4f  EDGE p = %.4f\n",
              lam, G, c3$p_a, c3$p_chi2, c3$p_E, c2$p_E, c3$b, c3$Hbar, cf$SC.HL$p.value, cf$SC.EDGE.adaptive$p.value))
}
c10 <- cav_gof_fast(Xg, yg, 295.3, G = 10)
cat("\nObserved vs expected by cavity decile at lambda_1se (G = 10):\n")
idx <- split(seq_len(n), c10$groups)
print(round(data.frame(g = seq_along(idx), n = sapply(idx, length), obs = sapply(idx, function(I) sum(yg[I])),
                       exp = sapply(idx, function(I) sum(c10$m[I])), e = c10$e), 2), row.names = FALSE)
