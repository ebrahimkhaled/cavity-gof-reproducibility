HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
x <- read.csv(file.path(ROOT, "results", "E50_tau_perrep.csv"))
r <- function(p, a = .05) mean(p < a, na.rm = TRUE)
pd <- function(a, b) { d <- (a < .05) - (b < .05); d <- d[!is.na(d)]; if (!length(d)) return(c(NA, NA)); c(mean(d), sd(d) / sqrt(length(d))) }
for (ix in c("I1", "I2", "full")) { cat("\n==== index", ix, "====\n"); g <- function(d, s) d[[paste0(ix, "_", s)]]
  tab <- do.call(rbind, lapply(unique(x$cell), function(cn) { d <- x[x$cell == cn, ]; a <- pd(g(d, "R2"), g(d, "d3")); b <- pd(g(d, "XF"), g(d, "d3"))
    data.frame(cell = cn, R = sum(!is.na(g(d, "d3"))), d3 = r(g(d, "d3")), T0 = r(g(d, "T0")), OR = r(g(d, "OR")), XF = r(g(d, "XF")), R2 = r(g(d, "R2")),
               R2_1pc = r(g(d, "R2"), .01), d3_1pc = r(g(d, "d3"), .01), R2md3 = a[1], se = a[2], XFmd3 = b[1], se2 = b[2],
               rhoOR = mean(g(d, "rho_OR"), na.rm = TRUE), rhoXF = mean(g(d, "rho_XF"), na.rm = TRUE), rhoR2 = mean(g(d, "rho_R2"), na.rm = TRUE),
               tauOR = mean(g(d, "tau_OR"), na.rm = TRUE), sdOR = sd(g(d, "tau_OR"), na.rm = TRUE), tauR2 = mean(g(d, "tau_R2"), na.rm = TRUE), sdR2 = sd(g(d, "tau_R2"), na.rm = TRUE),
               tauXF = mean(g(d, "tau_XF"), na.rm = TRUE), sdXF = sd(g(d, "tau_XF"), na.rm = TRUE), medD3 = median(g(d, "tau_d3"), na.rm = TRUE), sdD3 = sd(g(d, "tau_d3"), na.rm = TRUE),
               edge = mean(g(d, "tau_d3") < .002 | g(d, "tau_d3") > 19.9, na.rm = TRUE)) }))
  print(format(tab, digits = 3), row.names = FALSE)
  write.csv(tab, file.path(ROOT, "results", paste0("E50_summary_", ix, ".csv")), row.names = FALSE) }
o <- read.csv(file.path(ROOT, "results", "E47c_std_perrep.csv")); m <- merge(x, o, by = c("rep", "cell"))
cat("\nfingerprint: max |d3 here - E47c| I1", max(abs(m$I1_d3 - m$I1), na.rm = TRUE), " full", max(abs(m$full_d3 - m$full), na.rm = TRUE), "\n")
