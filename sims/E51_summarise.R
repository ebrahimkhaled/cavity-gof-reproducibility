HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
x <- read.csv(file.path(ROOT, "results", "E51_confirm_perrep.csv"))
r <- function(p, a = .05) mean(p < a, na.rm = TRUE)
pd <- function(a, b) { d <- (a < .05) - (b < .05); d <- d[!is.na(d)]; if (!length(d)) return(c(NA, NA)); c(mean(d), sd(d) / sqrt(length(d))) }
ms <- function(p) mean(qchisq(p, 8, lower.tail = FALSE), na.rm = TRUE)
for (ix in c("I1", "I2", "full")) { cat("\n==== index", ix, "====\n"); g <- function(d, s) d[[paste0(ix, "_", s)]]
  tab <- do.call(rbind, lapply(unique(x$cell), function(cn) { d <- x[x$cell == cn, ]; if (all(is.na(g(d, "d3")))) return(NULL); a <- pd(g(d, "R2"), g(d, "d3"))
    data.frame(cell = cn, R = sum(!is.na(g(d, "d3"))), nNA_R2 = sum(is.na(g(d, "R2"))), d3 = r(g(d, "d3")), T0 = r(g(d, "T0")), OR = r(g(d, "OR")), R2 = r(g(d, "R2")),
               R2md3 = a[1], se = a[2], S_d3 = ms(g(d, "d3")), S_T0 = ms(g(d, "T0")), S_R2 = ms(g(d, "R2")),
               rhoOR = mean(g(d, "rho_OR"), na.rm = TRUE), rhoR2 = mean(g(d, "rho_R2"), na.rm = TRUE),
               tauOR = median(g(d, "tau_OR"), na.rm = TRUE), tauR2 = median(g(d, "tau_R2"), na.rm = TRUE), sdR2 = sd(g(d, "tau_R2"), na.rm = TRUE),
               medD3 = median(g(d, "tau_d3"), na.rm = TRUE), sdD3 = sd(g(d, "tau_d3"), na.rm = TRUE)) }))
  print(format(tab, digits = 3), row.names = FALSE)
  write.csv(tab, file.path(ROOT, "results", paste0("E51_summary_", ix, ".csv")), row.names = FALSE) }
