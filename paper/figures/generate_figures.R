## Generate paper figures from precomputed results
## Run from package root: Rscript figures/generate_figures.R

dir.create("figures", showWarnings = FALSE)

# -----------------------------------------------------------------------
# Figure 1: Mechanism comparison boxplot
# -----------------------------------------------------------------------
exp1 <- readRDS("inst/precomputed/paper/exp1_main_comparison.rds")

# Reshape to long format for boxplot
mechs <- c("scheme0", "ordering", "heterogeneous",
           "partial_autopsy", "masked", "periodic")
labels <- c("Scheme 0", "Ordering", "Het. fleet",
            "Partial\nautopsy", "Masking", "Periodic\ninspection")

vals <- lapply(mechs, function(m) exp1[[m]])
names(vals) <- labels

pdf("figures/fig1_mechanism_comparison.pdf", width = 7, height = 4.5)
par(mar = c(5.5, 4.5, 1.5, 1), cex.lab = 1.1)
boxplot(vals,
        ylab = "MAE (sorted rates)",
        col = c("gray85", "gray85", "gray85",
                "steelblue2", "steelblue2", "steelblue2"),
        border = "gray30",
        outline = TRUE,
        las = 2,
        ylim = c(0, max(unlist(vals)) * 1.05))
abline(h = 0.200, lty = 2, col = "gray50", lwd = 0.8)
text(0.7, 0.210, "Symmetric point", cex = 0.75, col = "gray40", adj = 0)
dev.off()
cat("Wrote figures/fig1_mechanism_comparison.pdf\n")

# -----------------------------------------------------------------------
# Figure 2: Autopsy coverage (MAE vs r)
# -----------------------------------------------------------------------
exp6 <- readRDS("inst/precomputed/paper/exp6_autopsy_coverage.rds")

pdf("figures/fig2_autopsy_coverage.pdf", width = 5, height = 4)
par(mar = c(4.5, 4.5, 1.5, 1), cex.lab = 1.1)
plot(exp6$r, exp6$median_mae,
     type = "b", pch = 19, lwd = 2, col = "steelblue3",
     xlab = "Number of inspected components (r)",
     ylab = "Median MAE",
     ylim = c(0, max(exp6$median_mae) * 1.1),
     xaxt = "n")
axis(1, at = exp6$r)
# Highlight the phase transition at r=1
segments(x0 = 0, y0 = exp6$median_mae[exp6$r == 0],
         x1 = 1, y1 = exp6$median_mae[exp6$r == 1],
         col = "firebrick", lwd = 2.5, lty = 1)
points(c(0, 1),
       exp6$median_mae[exp6$r %in% c(0, 1)],
       pch = 19, col = "firebrick", cex = 1.3)
text(0.6, 0.14, expression(paste("91% drop at ", italic(r) == 1)),
     cex = 0.85, col = "firebrick")
dev.off()
cat("Wrote figures/fig2_autopsy_coverage.pdf\n")

# -----------------------------------------------------------------------
# Figure 3: Delta insensitivity
# -----------------------------------------------------------------------
exp5 <- readRDS("inst/precomputed/paper/exp5_inspection_delta.rds")

pdf("figures/fig3_delta_sensitivity.pdf", width = 5, height = 4)
par(mar = c(4.5, 4.5, 1.5, 1), cex.lab = 1.1)
plot(exp5$delta, exp5$median_mae,
     type = "b", pch = 19, lwd = 2, col = "steelblue3",
     xlab = expression(paste("Inspection interval (", delta, ")")),
     ylab = "Median MAE",
     ylim = c(0, max(exp5$mean_mae) * 1.3),
     log = "x",
     xaxt = "n")
axis(1, at = exp5$delta, labels = exp5$delta)
# Add mean MAE line
lines(exp5$delta, exp5$mean_mae,
      type = "b", pch = 17, lwd = 1.5, col = "coral3", lty = 2)
legend("topright", legend = c("Median MAE", "Mean MAE"),
       pch = c(19, 17), col = c("steelblue3", "coral3"),
       lty = c(1, 2), lwd = c(2, 1.5), cex = 0.85, bty = "n")
dev.off()
cat("Wrote figures/fig3_delta_sensitivity.pdf\n")
