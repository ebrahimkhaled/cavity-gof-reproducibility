# Post-hoc look at the drift cells of E51, after the verdict rule was read. Descriptive only: it does not change the pre-set verdict.
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
x <- read.csv(file.path(ROOT, "results", "E51_confirm_perrep.csv"))
S <- function(p) qchisq(p, 8, lower.tail = FALSE)
for (cn in c("bent_n2000", "bent_n5000", "bent_n20000", "k25_n3200")) for (ix in c("I1", "I2", "full")) { d <- x[x$cell == cn, ]; g <- function(s) d[[paste0(ix, "_", s)]]
  if (all(is.na(g("d3")))) next
  f <- function(s) { v <- S(g(s)); sprintf("%s size %.3f (se %.3f) 1%% %.3f meanS %.2f (se %.2f)", s, mean(g(s) < .05), sqrt(mean(g(s) < .05) * (1 - mean(g(s) < .05)) / nrow(d)), mean(g(s) < .01), mean(v), sd(v) / sqrt(nrow(d))) }
  dd <- S(g("T0")) - S(g("R2")); cat(cn, ix, "\n  ", f("d3"), "\n  ", f("T0"), "\n  ", f("R2"), "\n   paired S(T0) - S(R2):", round(mean(dd), 3), "se", round(sd(dd) / sqrt(length(dd)), 3), "\n") }
