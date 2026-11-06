HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
x <- read.csv(file.path(ROOT, "results", "E47c_std_perrep.csv")); o <- read.csv(file.path(ROOT, "results", "E43_split_q3_perrep.csv"))
r <- function(p, a = .05) mean(p < a, na.rm = TRUE)
pd <- function(a, b) { d <- (a < .05) - (b < .05); d <- d[!is.na(d)]; if (!length(d)) return(c(NA, NA)); c(mean(d), sd(d) / sqrt(length(d))) }
tab <- do.call(rbind, lapply(unique(x$cell), function(cn) { d <- x[x$cell == cn, ]
  a <- pd(d$I1, d$full); b <- pd(d$I2, d$F1); cc <- pd(d$I2, d$full)
  data.frame(cell = cn, R = nrow(d), I1 = r(d$I1), I2 = r(d$I2), full = r(d$full), F1 = r(d$F1), I1_01 = r(d$I1, .01), I2_01 = r(d$I2, .01), full_01 = r(d$full, .01),
             I1mFull = a[1], se = a[2], I2mF1 = b[1], se2 = b[2], I2mFull = cc[1], se3 = cc[2], NAs = sum(is.na(d$I1))) }))
print(format(tab, digits = 3), row.names = FALSE)
m <- merge(x, o[, c("rep", "cell", "full_q3")], by = c("rep", "cell")); cat("std vs stored default: same decision share", mean((m$full < .05) == (m$full_q3 < .05), na.rm = TRUE), "\n")
write.csv(tab, file.path(ROOT, "results", "E47c_summary.csv"), row.names = FALSE)
