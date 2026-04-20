# ============================================================
# Figure 2 (corrected 4-panel version)
# One VPD per panel, uncertainty in fecundity + host-finding
# ============================================================

rm(list = ls()); gc()

# ---------------------------
# 0) CONSTANTS
# ---------------------------
F_mean <- 6867
T_min_activity <- 12

D_host_mean <- 30
D_host_min  <- 15
D_host_max  <- 45

# ---------------------------
# 1) TRIANGULAR SAMPLER
# ---------------------------
rtri <- function(n, a, c, b) {
  u <- runif(n)
  Fc <- (c - a) / (b - a)
  ifelse(
    u < Fc,
    a + sqrt(u * (b - a) * (c - a)),
    b - sqrt((1 - u) * (b - a) * (b - c))
  )
}

# ---------------------------
# 2) TRAIT FUNCTIONS
# ---------------------------
Tp <- function(T, VD) 50.775 - 1.0728*T - 0.229*VD
To <- function(T, VD) 66.19  - 1.5638*T - 0.179*VD
Ti <- function(T, VD) 59.7   - 1.151*T  - 0.1014*VD
Tm <- function(T, VD) 192.23 - 6.054*T + 0.258*VD

Mf_pct <- function(T, VD) 19.32   + 1.212*T - 0.16*VD
Me_pct <- function(T, VD) 108.325 - 3.848*T + 1.414*VD
Mn_pct <- function(T, VD) 51.4786 + 1.525*T - 0.22*VD

Ml_weeks <- function(T, VD) 16.1 + 0.814*T - 0.21*VD
Ma_weeks <- function(T, VD) 14.3 + 0.792*T - 0.14*VD

# keep the value you decided to use
Ql_raw <- function(T) 0.0063 * (T^0.51)
Qa_raw <- function(T) 0.051  * (T^0.49)

# ---------------------------
# 3) HELPERS
# ---------------------------
clip01 <- function(p) pmax(0, pmin(1, p))

mu_from_stage <- function(mort_pct, duration) {
  q <- clip01(mort_pct / 100)
  out <- ifelse(q >= 1, 1e6, -log(1 - q) / duration)
  out[!is.finite(out) | duration <= 0] <- NA
  out
}

mu_from_survival <- function(weeks) {
  out <- 1 / (weeks * 7)
  out[!is.finite(out) | weeks <= 0] <- NA
  out
}

delta_from_time <- function(time) {
  out <- 1 / time
  out[!is.finite(out) | time <= 0] <- NA
  out
}

Q_rule <- function(Qraw, T, D_host = 30) {
  Q_min <- 1 / D_host
  ifelse(T > T_min_activity, pmax(Qraw, Q_min), 0)
}

# ---------------------------
# 4) R0 FUNCTION
# ---------------------------
R0_tick <- function(T, VD, F, D_host = 30) {
  
  tp <- Tp(T, VD)
  to <- To(T, VD)
  ti <- Ti(T, VD)
  tm <- Tm(T, VD)
  
  deltaE <- delta_from_time(ti)
  deltaN <- delta_from_time(tm)
  deltaF <- delta_from_time(tp + to)
  
  muE <- mu_from_stage(Me_pct(T, VD), ti)
  muN <- mu_from_stage(Mn_pct(T, VD), tm)
  muF <- mu_from_stage(Mf_pct(T, VD), tp + to)
  
  muL <- mu_from_survival(Ml_weeks(T, VD))
  muA <- mu_from_survival(Ma_weeks(T, VD))
  
  Ql <- Q_rule(Ql_raw(T), T, D_host = D_host)
  Qa <- Q_rule(Qa_raw(T), T, D_host = D_host)
  
  term <- (deltaF / (deltaF + muF)) *
    (Qa     / (Qa + muA))     *
    (deltaN / (deltaN + muN)) *
    (Ql     / (Ql + muL))     *
    (deltaE / (deltaE + muE))
  
  R0 <- F * term
  R0[!is.finite(R0) | R0 < 0] <- NA
  R0
}

# ---------------------------
# 5) MONTE CARLO WRAPPER
# ---------------------------
R0_MC <- function(T, VD, n = 1000) {
  
  F_samp      <- rtri(n, 3184, F_mean, 13180)
  D_host_samp <- rtri(n, D_host_min, D_host_mean, D_host_max)
  
  R0_mat <- sapply(seq_len(n), function(i) {
    R0_tick(T, VD, F = F_samp[i], D_host = D_host_samp[i])
  })
  
  list(
    mean  = rowMeans(R0_mat, na.rm = TRUE),
    lower = apply(R0_mat, 1, quantile, 0.025, na.rm = TRUE),
    upper = apply(R0_mat, 1, quantile, 0.975, na.rm = TRUE)
  )
}

# ---------------------------
# 6) TEMPERATURE RANGE
# ---------------------------
T_seq <- seq(0, 40, by = 0.5)

# ---------------------------
# 7) SELECT VPD VALUES
# ---------------------------
VD_vals <- c(0, 1.5, 2.5, 3.5)
cols <- c("#2c7bb6", "#9ecae1", "#fdae61", "#d7191c")
labels <- c("(a)", "(b)", "(c)", "(d)")

# ---------------------------
# 8) PRECOMPUTE TO GET SHARED Y-LIMIT
# ---------------------------
det_list <- vector("list", length(VD_vals))
mc_list  <- vector("list", length(VD_vals))

for (i in seq_along(VD_vals)) {
  det_list[[i]] <- R0_tick(T_seq, VD_vals[i], F = F_mean, D_host = D_host_mean)
  mc_list[[i]]  <- R0_MC(T_seq, VD_vals[i], n = 1000)
}

ymax <- max(
  sapply(det_list, max, na.rm = TRUE),
  sapply(mc_list, function(x) max(x$upper, na.rm = TRUE)),
  na.rm = TRUE
)

ymax <- ceiling(ymax / 100) * 100

# ---------------------------
# 9) PLOT
# ---------------------------
tiff(
  "Figure2_clean_4panel.tiff",
  width = 2600,
  height = 2200,
  res = 300,
  compression = "lzw"
)

par(
  mfrow = c(2, 2),
  mar   = c(4.5, 4.5, 3, 1),
  cex.axis = 1.1,
  cex.lab  = 1.3
)

for (i in seq_along(VD_vals)) {
  
  VD  <- VD_vals[i]
  col <- cols[i]
  det <- det_list[[i]]
  mc  <- mc_list[[i]]
  
  plot(
    NA, NA,
    xlim = range(T_seq),
    ylim = c(0, ymax),
    xlab = "Temperature (°C)",
    ylab = expression(R[0]),
    las  = 1,
    bty  = "l"
  )
  
  polygon(
    c(T_seq, rev(T_seq)),
    c(mc$lower, rev(mc$upper)),
    col = adjustcolor(col, alpha.f = 0.18),
    border = NA
  )
  
  lines(T_seq, mc$mean, col = col, lwd = 3)
  lines(T_seq, det, col = "black", lwd = 2.2, lty = 2)
  
  mtext(labels[i], side = 3, adj = 0, font = 2, cex = 1.2)
  mtext(paste("VPD =", VD), side = 3, line = -1, adj = 0.5, cex = 1.1)
  
  if (i == 1) {
    legend(
      "topleft",
      legend = c("Stochastic mean", "Deterministic", "95% CI"),
      col    = c(col, "black", "grey60"),
      lwd    = c(3, 2.2, NA),
      lty    = c(1, 2, NA),
      pch    = c(NA, NA, 15),
      pt.cex = c(NA, NA, 2),
      bty    = "n",
      cex    = 1
    )
  }
}

dev.off()

