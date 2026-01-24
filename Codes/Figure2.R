# ============================================================
# Figure 2
# Temperature-dependent R₀ with propagated biological uncertainty
# Trait-based mechanistic framework for Hyalomma marginatum
#
# This script contrasts deterministic and stochastic expectations
# of population growth by explicitly propagating uncertainty in
# fecundity and host-finding time through a nonlinear R₀ formulation.
# ============================================================
rm(list = ls()); gc()

# ---------------------------
# 0) CONSTANTS
# ---------------------------
F_mean <- 6867  # eggs per female (mean reported)
T_min_activity <- 12  # activity threshold stated in Text S1
# ============================================================
# UNCERTAINTY IN HOST-FINDING TIME (days to host)
# ============================================================

# Literature-informed plausible range (can be justified ecologically)
D_host_mean <- 30   # central assumption (used before)
D_host_min  <- 15   # high host density / favourable landscape
D_host_max  <- 45   # low host density / fragmented landscape

# Triangular sampler (same logic as fecundity)
rtri <- function(n, a, c, b) {
  u <- runif(n)
  Fc <- (c - a) / (b - a)
  ifelse(
    u < Fc,
    a + sqrt(u * (b - a) * (c - a)),
    b - sqrt((1 - u) * (b - a) * (b - c))
  )
}
# Optional: if you want to count only female offspring contributing to future females:
sex_ratio_female <- 1.0  # set to 0.5 if you want a 1:1 sex ratio assumption
F_eff <- F_mean * sex_ratio_female

# ---------------------------
# 1) ESTRADA-PEÑA (2011) TRAIT EQUATIONS
# ---------------------------
# Development times (days)
Tp <- function(T, VD) 50.775 - 1.0728*T - 0.229*VD         # pre-oviposition
To <- function(T, VD) 66.19  - 1.5638*T - 0.179*VD         # oviposition
Ti <- function(T, VD) 59.7   - 1.151*T  - 0.1014*VD        # incubation (egg dev time)
Tm <- function(T, VD) 192.23 - 6.054*T + 0.258*VD          # molt engorged nymphs -> adults

# Mortality (% over the stage period)
Mf_pct <- function(T, VD) 19.32   + 1.212*T - 0.16*VD      # females (fed/repro stage)
Me_pct <- function(T, VD) 108.325 - 3.848*T + 1.414*VD     # eggs (developing)
Mn_pct <- function(T, VD) 51.4786 + 1.525*T - 0.22*VD      # nymphs (molting)

# Survival time of questing stages (weeks)
Ml_weeks <- function(T, VD) 16.1 + 0.814*T - 0.21*VD       # questing larvae survival time
Ma_weeks <- function(T, VD) 14.3 + 0.792*T - 0.14*VD       # questing adults survival time

# Activity/host finding (as given)
Ql_raw <- function(T) 0.063 * (T^0.51)
Qa_raw <- function(T) 0.051 * (T^0.49)

# ---------------------------
# 2) HELPERS: SAFE TRANSFORMS
# ---------------------------
clip01 <- function(p) pmax(0, pmin(1, p))

# Convert "mortality percent over the stage period" into a per-day hazard μ
# If q = probability of dying over the whole stage, survival over stage is (1-q)
# assume constant hazard: survival = exp(-μ * duration) => μ = -ln(1-q)/duration
mu_from_stage_mortality <- function(mort_pct, duration_days) {
  q <- clip01(mort_pct / 100)
  # handle edge cases (q=1 -> infinite hazard; q=0 -> 0 hazard)
  out <- rep(NA_real_, length(q))
  ok <- is.finite(duration_days) & duration_days > 0 & is.finite(q)
  q_ok <- q[ok]
  d_ok <- duration_days[ok]
  # if q==1, set very large hazard
  out[ok] <- ifelse(q_ok >= 1, 1e6, -log(1 - q_ok) / d_ok)
  out
}

# Convert "survival time in weeks" into per-day mortality hazard for questing stages
# If mean survival time is S days and hazard is constant, μ ≈ 1/S
mu_from_survival_time_weeks <- function(S_weeks) {
  S_days <- S_weeks * 7
  out <- rep(NA_real_, length(S_days))
  ok <- is.finite(S_days) & S_days > 0
  out[ok] <- 1 / S_days[ok]
  out
}

# Convert "development time (days)" into per-day development rate δ = 1/time
delta_from_time <- function(time_days) {
  out <- rep(NA_real_, length(time_days))
  ok <- is.finite(time_days) & time_days > 0
  out[ok] <- 1 / time_days[ok]
  out
}

# Enforce activity rule: if T <= 12°C => no host-finding (Q=0)
# if T > 12°C => Q = max(Q_raw, 1/30) to ensure host can be found within 30 days
Q_with_rule <- function(Qraw, T, D_host = 30) {
  Q_min <- 1 / D_host
  Q <- ifelse(T > T_min_activity, pmax(Qraw, Q_min), 0)
  Q
}

# ---------------------------
# 3) BUILD R0(T,VD)
# ---------------------------
R0_tick <- function(T, VD, F = F_eff) {
  # stage durations
  tp <- Tp(T, VD)
  to <- To(T, VD)
  ti <- Ti(T, VD)
  tm <- Tm(T, VD)
  
  # development rates
  # δ_E: egg development (incubation)
  deltaE <- delta_from_time(ti)
  
  # δ_N: molting nymphs to adults (engorged nymph stage)
  deltaN <- delta_from_time(tm)
  
  # δ_f: "fed female reproductive progression" (pre-ovip + ovip)
  # (this matches your intended “female progression” term)
  deltaF <- delta_from_time(tp + to)
  
  # mortalities (per day)
  muE <- mu_from_stage_mortality(Me_pct(T, VD), ti)
  muN <- mu_from_stage_mortality(Mn_pct(T, VD), tm)
  muF <- mu_from_stage_mortality(Mf_pct(T, VD), tp + to)
  
  # questing mortalities (per day)
  muL <- mu_from_survival_time_weeks(Ml_weeks(T, VD))
  muA <- mu_from_survival_time_weeks(Ma_weeks(T, VD))
  
  # activity / host-finding rates (per day, with 30-day rule)
  Ql <- Q_with_rule(Ql_raw(T), T)
  Qa <- Q_with_rule(Qa_raw(T), T)
  
  # build R0 terms (each term is a probability-like factor)
  # δ/(δ+μ)  : survive development/repro vs die during that phase
  # Q/(Q+μ)  : survive questing long enough to find a host
  term_f <- deltaF / (deltaF + muF)
  term_a <- Qa / (Qa + muA)
  term_n <- deltaN / (deltaN + muN)
  term_l <- Ql / (Ql + muL)
  term_e <- deltaE / (deltaE + muE)
  
  R0 <- F * term_f * term_a * term_n * term_l * term_e
  
  # clean impossible values
  R0[!is.finite(R0) | R0 < 0] <- NA_real_
  R0
}


cat("\nDone. R0 surface computed.\n")
cat("F used =", F_eff, "eggs per female\n")
cat("Sex ratio multiplier =", sex_ratio_female, "\n")


# ============================================================
# MONTE CARLO UNCERTAINTY PROPAGATION (fecundity F only)
# Produces mean/low/high R0(T,VD) surfaces
# ============================================================

# --- 1) Fecundity summary from literature ---
F_mean <- 6867
F_min  <- 3184
F_max  <- 13180

# --- 2) Triangular sampler (min, mode, max) ---
rtri <- function(n, a, c, b) {
  # a=min, c=mode, b=max
  u <- runif(n)
  Fc <- (c - a) / (b - a)
  out <- ifelse(
    u < Fc,
    a + sqrt(u * (b - a) * (c - a)),
    b - sqrt((1 - u) * (b - a) * (b - c))
  )
  out
}

# --- 3) Monte Carlo settings ---
set.seed(123)
B <- 2000

# --- 4) Grid (use your same sequences) ---
T_seq  <- seq(0, 40, by = 0.5)
VD_seq <- seq(0, 50, by = 1)
grid   <- expand.grid(T = T_seq, VD = VD_seq)

# --- 5) Run Monte Carlo: sample F, compute R0 surface each time ---
R0_mat <- matrix(NA_real_, nrow = nrow(grid), ncol = B)

F_samp <- rtri(B, a = F_min, c = F_mean, b = F_max)

# Sample host-finding time uncertainty
D_host_samp <- rtri(B, a = D_host_min, c = D_host_mean, b = D_host_max)

for (b in 1:B) {
  
  # Override Q rule internally via closure
  R0_mat[, b] <- {
    
    # temporarily redefine Q rule
    Q_with_rule_mc <- function(Qraw, T) {
      Q_min <- 1 / D_host_samp[b]
      ifelse(T > T_min_activity, pmax(Qraw, Q_min), 0)
    }
    
    # copy of R0 with stochastic Q
    R0_tick(
      grid$T,
      grid$VD,
      F = F_samp[b]
    )
  }
}

# --- 6) Summaries per pixel ---
grid$R0_mean <- rowMeans(R0_mat, na.rm = TRUE)
grid$R0_low  <- apply(R0_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
grid$R0_high <- apply(R0_mat, 1, quantile, probs = 0.975, na.rm = TRUE)

cat("Monte Carlo done: uncertainty propagated over fecundity F only.\n")
summary(F_samp)

# ============================================================
# Monte Carlo wrapper for R0(T, VD)
# Uncertainty in fecundity only (as used in your plot)
# ============================================================

R0_MC <- function(T, VD, F_mean, F_sd = NULL, n = 1000) {
  
  # If SD not provided, infer from min–max range (triangular)
  if (is.null(F_sd)) {
    F_min <- 3184
    F_max <- 13180
    
    F_samp <- rtri(n, a = F_min, c = F_mean, b = F_max)
  } else {
    F_samp <- rnorm(n, mean = F_mean, sd = F_sd)
    F_samp[F_samp < 0] <- 0
  }
  
  # Compute R0 for each Monte Carlo draw
  R0_mat <- sapply(F_samp, function(Fi) {
    R0_tick(T, VD, F = Fi)
  })
  
  list(
    mean  = rowMeans(R0_mat, na.rm = TRUE),
    lower = apply(R0_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(R0_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

# ---------------------------
# SAVE FIGURE (300 dpi)
# ---------------------------
tiff(
  filename = "R0_deterministic_vs_stochastic.tiff",
  width = 2000,
  height = 1600,
  res = 300
)

plot(
  NA, NA,
  xlim = range(T_seq),
  ylim = c(0, 1800),
  xlab = "Temperature (°C)",
  ylab = expression(R[0]),
  main = expression("")
)

for (i in seq_along(VD_vals)) {
  
  VD <- VD_vals[i]
  col <- cols[i]
  
  # Deterministic mean
  det <- R0_tick(T_seq, VD, F = F_mean)
  
  # Monte Carlo
  mc <- R0_MC(T_seq, VD, F_mean, n = 1000)
  
  # Uncertainty band
  polygon(
    c(T_seq, rev(T_seq)),
    c(mc$lower, rev(mc$upper)),
    col = adjustcolor(col, alpha.f = 0.25),
    border = NA
  )
  
  # Stochastic mean
  lines(T_seq, mc$mean, col = col, lwd = 2, lty = 1)
  
  # Deterministic mean (dashed, drawn last)
  lines(T_seq, det, col = col, lwd = 3.5, lty = 2)
}

legend(
  "topleft",
  legend = c(
    "Solid = stochastic mean",
    "Dashed = deterministic mean",
    "Shaded = 95% uncertainty",
    paste("VD =", VD_vals)
  ),
  bty = "n"
)

dev.off()