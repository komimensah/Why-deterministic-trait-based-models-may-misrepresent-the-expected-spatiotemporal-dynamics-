# ============================================================
# Figure 3. Quarterly spatial comparison of deterministic
#           and stochastic R0 (North Africa)
# Also produces the log-difference distribution figure
# ============================================================

setwd("/Users/kagboka/Desktop/H.manginatrum/")

library(terra)
library(rnaturalearth)
library(sf)

# ------------------------------------------------------------
# 1) LOAD STUDY REGION: NORTH AFRICA
# ------------------------------------------------------------

world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_v  <- vect(world_sf)

north_africa_names <- c(
  "Morocco", "Algeria", "Tunisia",
  "Libya", "Egypt", "Western Sahara"
)

region <- world_v[world_v$name_long %in% north_africa_names, ]
region <- aggregate(region)

# ------------------------------------------------------------
# 2) LOAD QUARTERLY R0 RASTERS
# ------------------------------------------------------------

R0_det_q   <- rast("outputs/R0_det_quarterly.tif")
R0_stoch_q <- rast("outputs/R0_stoch_quarterly.tif")

names(R0_det_q)   <- paste0("Q", 1:4)
names(R0_stoch_q) <- paste0("Q", 1:4)

# ------------------------------------------------------------
# 3) COMPUTE LOG-RATIO DELTA R0
# ------------------------------------------------------------

R0_log_delta_q <- ifel(
  (R0_det_q > 0) & (R0_stoch_q > 0),
  log(R0_stoch_q / R0_det_q),
  NA
)

names(R0_log_delta_q) <- paste0("Q", 1:4)

# ------------------------------------------------------------
# 3B) DISTRIBUTIONAL ANALYSIS OF LOG DIFFERENCES
# ------------------------------------------------------------

all_log <- values(R0_log_delta_q)
all_log <- all_log[is.finite(all_log)]

log_mean   <- mean(all_log)
log_median <- median(all_log)
log_sd     <- sd(all_log)
prop_pos   <- mean(all_log > 0)

cat("\n===== LOG-RATIO DISTRIBUTION SUMMARY =====\n")
cat("Mean log difference: ", log_mean, "\n")
cat("Median log difference: ", log_median, "\n")
cat("SD log difference: ", log_sd, "\n")
cat("Absolute mean: ", abs(log_mean), "\n")
cat("Proportion positive: ", prop_pos, "\n")
cat("\nQuantiles:\n")
print(quantile(all_log, probs = c(0.01, 0.25, 0.5, 0.75, 0.99)))

# ------------------------------------------------------------
# 3C) DISTRIBUTION FIGURE
# ------------------------------------------------------------

tiff(
  "Figure_LogDelta_Distribution.tiff",
  width = 3600,
  height = 1800,
  res = 300,
  compression = "lzw"
)

par(
  mfrow = c(1, 2),
  mar   = c(5, 5, 3, 2),
  cex.axis = 1.5,
  cex.lab  = 1.6
)

# (a) Histogram
hist_obj <- hist(
  all_log,
  breaks = 60,
  col = "#d9d9d9",
  border = "white",
  main = "",
  xlab = expression(log(R[0]^{"stoch"} / R[0]^{"det"})),
  ylab = "Frequency"
)

abline(v = 0, col = "#d7191c", lwd = 3)
abline(v = log_mean, col = "#2c7bb6", lwd = 3)

mtext("(a)", side = 3, line = 0.8, adj = 0, font = 2, cex = 1.8)

text(
  x = log_mean,
  y = max(hist_obj$counts) * 0.88,
  labels = paste0("Mean = ", round(log_mean, 3)),
  col = "#2c7bb6",
  pos = 4,
  cex = 1.3
)

# (b) Density
dens <- density(all_log)

plot(
  dens,
  lwd = 3,
  col = "#2c7bb6",
  main = "",
  xlab = expression(log(R[0]^{"stoch"} / R[0]^{"det"})),
  ylab = "Density"
)

abline(v = 0, col = "#d7191c", lwd = 3)
abline(v = log_mean, col = "#2c7bb6", lwd = 3, lty = 2)

mtext("(b)", side = 3, line = 0.8, adj = 0, font = 2, cex = 1.8)

usr <- par("usr")
text(
  x = usr[1] + 0.03 * (usr[2] - usr[1]),
  y = usr[4] - 0.08 * (usr[4] - usr[3]),
  labels = paste0("Mean = ", round(log_mean, 3),
                  "\nSD = ", round(log_sd, 3)),
  adj = c(0, 1),
  cex = 1.3
)

dev.off()

# ------------------------------------------------------------
# 4) DEFINE CONSISTENT COLOUR SCALES
# ------------------------------------------------------------

r0_lim <- range(values(R0_det_q), values(R0_stoch_q), na.rm = TRUE)

log_vals <- values(R0_log_delta_q)
log_vals <- log_vals[is.finite(log_vals)]

# Use actual observed log-ratio range (not symmetric around zero)
log_min  <- quantile(log_vals, 0.02, na.rm = TRUE)
log_max  <- quantile(log_vals, 0.98, na.rm = TRUE)
log_zlim <- c(log_min, log_max)

# Sequential palette for R0
pal_r0 <- colorRampPalette(c(
  "#ffffe5", "#fee391", "#fec44f",
  "#fe9929", "#d95f0e", "#993404"
))(100)

# Sequential palette for positive log-ratio values
pal_log <- colorRampPalette(c(
  "#f7f7f7", "#f4cccc", "#ea9999", "#e06666", "#cc0000"
))(100)

# ------------------------------------------------------------
# 5) EXPORT FIGURE 3
# ------------------------------------------------------------

tiff(
  filename = "Figure_R0_Quarterly_Comparison_NorthAfrica.tiff",
  width    = 5600,
  height   = 6800,
  res      = 300,
  compression = "lzw"
)

par(
  mfrow = c(4, 3),
  mar   = c(3.5, 3.5, 4, 1.5),
  oma   = c(0, 5, 0, 0),
  cex.axis = 1.4,
  cex.lab  = 1.5
)

row_labels <- letters[1:4]

for (q in 1:4) {
  
  R0_det       <- R0_det_q[[q]]
  R0_stoch     <- R0_stoch_q[[q]]
  R0_log_delta <- R0_log_delta_q[[q]]
  
  # Deterministic
  plot(
    R0_det,
    col  = pal_r0,
    zlim = r0_lim,
    main = paste("Q", q, "Deterministic R0"),
    axes = FALSE,
    cex.main = 1.5,
    legend = FALSE
  )
  plot(region, add = TRUE, border = "black", lwd = 1)
  mtext(row_labels[q], side = 3, line = 0.5, adj = 0, font = 2, cex = 2)
  
  # Stochastic
  plot(
    R0_stoch,
    col  = pal_r0,
    zlim = r0_lim,
    main = paste("Q", q, "Stochastic mean R0"),
    axes = FALSE,
    cex.main = 1.5,
    legend = FALSE
  )
  plot(region, add = TRUE, border = "black", lwd = 1)
  
  # Log-ratio
  plot(
    R0_log_delta,
    col  = pal_log,
    zlim = log_zlim,
    main = expression(log(R[0]^{"stoch"} / R[0]^{"det"})),
    axes = FALSE,
    cex.main = 1.4,
    legend = FALSE
  )
  plot(region, add = TRUE, border = "black", lwd = 1)
}

# ---------------------------
# GLOBAL COLOR LEGENDS
# ---------------------------

# R0 color bar
par(fig = c(0.18, 0.45, 0.01, 0.045), new = TRUE, mar = c(2, 2, 1, 2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

ncol_r0 <- length(pal_r0)
x_r0 <- seq(0, 1, length.out = ncol_r0 + 1)

for (i in 1:ncol_r0) {
  rect(x_r0[i], 0.35, x_r0[i + 1], 0.75, col = pal_r0[i], border = NA)
}

axis(
  1,
  at = seq(0, 1, length.out = 5),
  labels = round(seq(r0_lim[1], r0_lim[2], length.out = 5), 1),
  cex.axis = 1.2,
  lwd = 0,
  lwd.ticks = 1
)
mtext("R0", side = 3, line = 0.2, cex = 1.2)

# log-ratio color bar
par(fig = c(0.55, 0.82, 0.01, 0.045), new = TRUE, mar = c(2, 2, 1, 2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

ncol_log <- length(pal_log)
x_log <- seq(0, 1, length.out = ncol_log + 1)

for (i in 1:ncol_log) {
  rect(x_log[i], 0.35, x_log[i + 1], 0.75, col = pal_log[i], border = NA)
}

axis(
  1,
  at = seq(0, 1, length.out = 5),
  labels = round(seq(log_zlim[1], log_zlim[2], length.out = 5), 2),
  cex.axis = 1.2,
  lwd = 0,
  lwd.ticks = 1
)
mtext(expression(log(R[0]^{"stoch"} / R[0]^{"det"})), side = 3, line = 0.2, cex = 1.1)

dev.off()

cat("✓ Figure 3 exported with global color legends\n")
cat("✓ Figure 3 exported with larger fonts\n")
cat("✓ Figure_LogDelta_Distribution exported with mean and SD annotation\n")
cat("Mean =", round(log_mean, 3), " | SD =", round(log_sd, 3), "\n")