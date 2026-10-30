HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
d <- read.csv(file.path(ROOT, "results", "E23_verify_prop1_perobs_seed20261103.csv"))
for (cc in unique(d$cell)) {
  s <- d[d$cell == cc, ]; n <- nrow(s)
  f <- function(v) sprintf("%+.4f (SE %.4f)", mean(v), sd(v) / sqrt(n))
  cat(sprintf("[%s] y - pi_true: %s | pi_true - Lambda(xi_exact): %s | y - Lambda: %s | mean pi_true %.4f, mean Lambda %.4f\n",
              cc, f(s$y - s$pi_true), f(s$pi_true - s$Lam_exact), f(s$y - s$Lam_exact), mean(s$pi_true), mean(s$Lam_exact)))
  # regression of the true probability on the exact cavity index: Proposition 1 says E(pi_true | xi) = Lambda(xi)
  br <- quantile(s$xi_exact, 0:10 / 10); g <- cut(s$xi_exact, br, include.lowest = TRUE)
  print(round(data.frame(decile = 1:10, pi_true = tapply(s$pi_true, g, mean), Lambda = tapply(s$Lam_exact, g, mean),
                         diff = tapply(s$pi_true - s$Lam_exact, g, mean), se = tapply(s$pi_true - s$Lam_exact, g, function(v) sd(v) / sqrt(length(v)))), 4), row.names = FALSE)
}
