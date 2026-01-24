# ============================================================
# Figure 3. Quarterly spatial comparison of deterministic
#           and stochastic R0 (North Africa)
#
# Description:
# This script generates Figure 3, showing quarterly (Q1–Q4)
# spatial maps of:
#   (a) Deterministic R0
#   (b) Stochastic mean R0
#   (c) Delta R0 (stochastic − deterministic)
#
# All maps use consistent colour scales to allow direct
# spatial comparison across quarters.
#
# Requirements:
# - R >= 4.1
# - terra, rnaturalearth, sf
#
# Inputs:
# - outputs/R0_det_quarterly.tif
# - outputs/R0_stoch_quarterly.tif
#
# Output:
# - Figure_R0_Quarterly_Comparison_NorthAfrica.tiff
# ============================================================

library(terra)
library(rnaturalearth)
library(sf)

# ------------------------------------------------------------
# 1) LOAD STUDY REGION: NORTH AFRICA
# ------------------------------------------------------------

# Download Admin-0 boundaries
world_sf <- ne_countries(scale = "medium", returnclass = "sf")
world_v  <- vect(world_sf)

# Countries defining North Africa
north_africa_names <- c(
  "Morocco", "Algeria", "Tunisia",
  "Libya", "Egypt", "Western Sahara"
)

# Subset region
region <- world_v[world_v$name_long %in% north_africa_names, ]

# ------------------------------------------------------------
# 2) LOAD QUARTERLY R0 RASTERS
# ------------------------------------------------------------

# Deterministic R0 (Q1–Q4)
R0_det_q <- rast(
  "outputs/R0_det_quarterly.tif"
)

# Stochastic mean R0 (Q1–Q4)
R0_stoch_q <- rast(
  "outputs/R0_stoch_quarterly.tif"
)

names(R0_det_q)   <- paste0("Q", 1:4)
names(R0_stoch_q) <- paste0("Q", 1:4)

# ------------------------------------------------------------
# 3) DEFINE CONSISTENT COLOUR SCALES
# ------------------------------------------------------------

# Shared limits for deterministic & stochastic R0
r0_lim <- range(
  values(R0_det_q),
  values(R0_stoch_q),
  na.rm = TRUE
)

# Symmetric limits for delta R0
delta_lim <- max(
  abs(values(R0_stoch_q - R0_det_q)),
  na.rm = TRUE
)

# Sequential palette for R0
pal_r0 <- colorRampPalette(c(
  "#ffffe5", "#fee391", "#fec44f",
  "#fe9929", "#d95f0e", "#993404"
))(100)

# Diverging palette for delta R0
pal_diff <- colorRampPalette(c(
  "#2c7bb6", "#f7f7f7", "#d7191c"
))(100)

# ------------------------------------------------------------
# 4) EXPORT FIGURE 3
# ------------------------------------------------------------

tiff(
  filename = "Figure_R0_Quarterly_Comparison_NorthAfrica.tiff",
  width    = 4200,
  height   = 4800,
  res      = 300
)

# 4 rows (Q1–Q4) × 3 columns (Det / Stoch / Delta)
par(
  mfrow = c(4, 3),
  mar   = c(2, 2, 3, 4),
  oma   = c(0, 4, 0, 0)
)

row_labels <- letters[1:4]  # a–d

for (q in 1:4) {
  
  R0_det   <- R0_det_q[[q]]
  R0_stoch <- R0_stoch_q[[q]]
  R0_delta <- R0_stoch - R0_det
  
  # ---- Deterministic R0 ----
  plot(
    R0_det,
    col   = pal_r0,
    zlim = r0_lim,
    main = paste("Q", q, "Deterministic R0"),
    axes = FALSE
  )
  plot(region, add = TRUE, border = "black", lwd = 0.9)
  
  # Panel label (a–d)
  mtext(
    row_labels[q],
    side = 3,
    line = 0.2,
    adj  = 0,
    font = 2,
    cex  = 1.6
  )
  
  # ---- Stochastic mean R0 ----
  plot(
    R0_stoch,
    col   = pal_r0,
    zlim = r0_lim,
    main = paste("Q", q, "Stochastic mean R0"),
    axes = FALSE
  )
  plot(region, add = TRUE, border = "black", lwd = 0.9)
  
  # ---- Delta R0 ----
  plot(
    R0_delta,
    col   = pal_diff,
    zlim = c(-delta_lim, delta_lim),
    main = expression(Delta*R[0]),
    axes = FALSE
  )
  plot(region, add = TRUE, border = "black", lwd = 0.9)
}

dev.off()

cat("✓ Figure 3 exported: Quarterly deterministic vs stochastic R0 comparison\n")
