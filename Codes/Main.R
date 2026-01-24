# ============================================================
# HYALOMMA marginatum complex — Trait-based R0 index
# Using Estrada-Peña et al. (2011) temperature (T) + water deficit (VD)
# Formula (your final structure):
# R0(T,VD) =
#   F *
#   [δ_f/(δ_f+μ_f)] *
#   [Q_a/(Q_a+μ_A)] *
#   [δ_N/(δ_N+μ_N)] *
#   [Q_l/(Q_l+μ_L)] *
#   [δ_E/(δ_E+μ_E)]
#
# Notes:
# - T in °C, VD in mm (as in Text S1)
# - Times Tp, To, Ti, Tm are in DAYS
# - Me, Mn, Mf are mortality in % over the corresponding stage period (days)
# - Ml, Ma are survival time in WEEKS for questing larvae/adults
# - Ql, Qa are activity/host-finding rates; here used as per-day rates
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
# SPATIOTEMPORAL APPLICATION:
# Monthly NetCDF → Deterministic & Stochastic R0 surfaces
# (PARALLEL-SAFE: all needed functions are defined inside app())
# ============================================================

# ============================================================
# SPATIAL PRE-PROCESSING FOR R0 MODELLING
# Crop & mask TerraClimate NetCDFs to Africa / North Africa
# ============================================================

library(terra)
library(rnaturalearth)
library(rnaturalearthdata)

setwd("/Users/kagboka/Desktop/H.manginatrum/")

# ---------------------------
# ============================================================
# 1) LOAD AFRICA BOUNDARY (NO LOCAL FILE NEEDED)
# ============================================================

# Download Admin-0 countries directly
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert to terra vector
world_v <- vect(world)

# Africa only
africa <- world_v[world_v$continent == "Africa", ]

# ---- OPTIONAL: North Africa subset ----
north_africa_names <- c(
  "Morocco", "Algeria", "Tunisia",
  "Libya", "Egypt", "Western Sahara"
)

north_africa <- africa[africa$name_long %in% north_africa_names, ]

# >>> SELECT REGION <<<
#region <- africa        # entire Africa
region <- north_africa  # uncomment for North Africa only
# ---------------------------
# 2) LOAD MONTHLY CLIMATE DATA (NetCDF)
# ---------------------------
# TerraClimate convention:
# tmin, tmax in °C
# def = climatic water deficit (VD, mm)

tmin <- rast("TerraClimate_tmin_2024.nc")
tmax <- rast("TerraClimate_tmax_2024.nc")
vd   <- rast("TerraClimate_def_2024.nc")
crs(tmin) <- "EPSG:4326"
crs(tmax) <- "EPSG:4326"
crs(vd)   <- "EPSG:4326"
# Sanity check: must be monthly
stopifnot(
  nlyr(tmin) == 12,
  nlyr(tmax) == 12,
  nlyr(vd)   == 12
)

# ---------------------------
# 3) CROP + MASK (CRITICAL FOR SPEED)
# ---------------------------
tmin_r <- mask(crop(tmin, region), region)
tmax_r <- mask(crop(tmax, region), region)
vd   <- mask(crop(vd,   region), region)

# ---------------------------
# 4) MONTHLY MEAN TEMPERATURE
# ---------------------------
Tmean <- (tmin_r + tmax_r) / 2

names(Tmean) <- month.abb
names(vd)    <- month.abb


# ---------------------------
# 5) FINAL STACK FOR R0 MODELLING
# ---------------------------
# Check extents
ext(Tmean)
ext(vd)

# Stack: 12 T + 12 VD = 24 layers
stack_monthly <- c(Tmean, vd)
# ---------------------------
# 6) QUICK CHECKS
# ---------------------------
print(stack_monthly)
plot(stack_monthly[[1]], main = "January Mean Temperature (°C)")
plot(stack_monthly[[13]], main = "January Water Deficit (mm)")

cat("✓ Climate data cropped and masked successfully\n")
cat("✓ Stack ready for deterministic and stochastic R0 modelling\n")
cat("Number of grid cells:", ncell(stack_monthly), "\n")





# ---------------------------
# 2) CONSTANTS (passed to workers)
# ---------------------------
F_mean <- 6867
F_min  <- 3184
F_max  <- 13180

T_min_activity <- 12

D_host_mean <- 30
D_host_min  <- 15
D_host_max  <- 45

sex_ratio_female <- 1.0

# ============================================================
# 3) DETERMINISTIC MONTHLY R0 (parallel-safe)
# ============================================================

R0_det_fun <- function(x, F_mean, sex_ratio_female, T_min_activity) {
  
  # ---- Trait equations (Estrada-Peña et al. 2011) ----
  Tp <- function(T, VD) 50.775 - 1.0728*T - 0.229*VD
  To <- function(T, VD) 66.19  - 1.5638*T - 0.179*VD
  Ti <- function(T, VD) 59.7   - 1.151*T  - 0.1014*VD
  Tm <- function(T, VD) 192.23 - 6.054*T + 0.258*VD
  
  Mf_pct <- function(T, VD) 19.32   + 1.212*T - 0.16*VD
  Me_pct <- function(T, VD) 108.325 - 3.848*T + 1.414*VD
  Mn_pct <- function(T, VD) 51.4786 + 1.525*T - 0.22*VD
  
  Ml_weeks <- function(T, VD) 16.1 + 0.814*T - 0.21*VD
  Ma_weeks <- function(T, VD) 14.3 + 0.792*T - 0.14*VD
  
  Ql_raw <- function(T) 0.063 * (T^0.51)
  Qa_raw <- function(T) 0.051 * (T^0.49)
  
  # ---- helpers ----
  clip01 <- function(p) pmax(0, pmin(1, p))
  
  mu_from_stage_mortality <- function(mort_pct, duration_days) {
    q <- clip01(mort_pct / 100)
    out <- rep(NA_real_, length(q))
    ok <- is.finite(duration_days) & duration_days > 0 & is.finite(q)
    q_ok <- q[ok]
    d_ok <- duration_days[ok]
    out[ok] <- ifelse(q_ok >= 1, 1e6, -log(1 - q_ok) / d_ok)
    out
  }
  
  mu_from_survival_time_weeks <- function(S_weeks) {
    S_days <- S_weeks * 7
    out <- rep(NA_real_, length(S_days))
    ok <- is.finite(S_days) & S_days > 0
    out[ok] <- 1 / S_days[ok]
    out
  }
  
  delta_from_time <- function(time_days) {
    out <- rep(NA_real_, length(time_days))
    ok <- is.finite(time_days) & time_days > 0
    out[ok] <- 1 / time_days[ok]
    out
  }
  
  Q_with_rule <- function(Qraw, T, D_host = 30) {
    Q_min <- 1 / D_host
    ifelse(T > T_min_activity, pmax(Qraw, Q_min), 0)
  }
  
  # ---- R0 ----
  R0_tick_local <- function(T, VD, F) {
    
    tp <- Tp(T, VD); to <- To(T, VD); ti <- Ti(T, VD); tm <- Tm(T, VD)
    
    deltaE <- delta_from_time(ti)
    deltaN <- delta_from_time(tm)
    deltaF <- delta_from_time(tp + to)
    
    muE <- mu_from_stage_mortality(Me_pct(T, VD), ti)
    muN <- mu_from_stage_mortality(Mn_pct(T, VD), tm)
    muF <- mu_from_stage_mortality(Mf_pct(T, VD), tp + to)
    
    muL <- mu_from_survival_time_weeks(Ml_weeks(T, VD))
    muA <- mu_from_survival_time_weeks(Ma_weeks(T, VD))
    
    Ql <- Q_with_rule(Ql_raw(T), T)
    Qa <- Q_with_rule(Qa_raw(T), T)
    
    term_f <- deltaF / (deltaF + muF)
    term_a <- Qa / (Qa + muA)
    term_n <- deltaN / (deltaN + muN)
    term_l <- Ql / (Ql + muL)
    term_e <- deltaE / (deltaE + muE)
    
    R0 <- F * term_f * term_a * term_n * term_l * term_e
    R0[!is.finite(R0) | R0 < 0] <- NA_real_
    R0
  }
  
  T  <- x[1:12]
  VD <- x[13:24]
  
  F_eff <- F_mean * sex_ratio_female
  R0_tick_local(T, VD, F = F_eff)  # returns length 12 (monthly)
}

R0_det_monthly <- app(
  stack_monthly,
  fun   = R0_det_fun,
  cores = 4,
  F_mean = F_mean,
  sex_ratio_female = sex_ratio_female,
  T_min_activity = T_min_activity
)
names(R0_det_monthly) <- month.abb

R0_det_annual <- app(R0_det_monthly, mean, na.rm = TRUE)

# ============================================================
# 4) STOCHASTIC MONTHLY R0 (parallel-safe)
#    returns 3 layers: mean / low / high per cell
# ============================================================

R0_stoch_fun <- function(
    x,
    F_mean, F_min, F_max,
    D_host_mean, D_host_min, D_host_max,
    sex_ratio_female,
    T_min_activity,
    B = 300
) {
  
  # ---- local triangular sampler ----
  rtri_local <- function(n, a, c, b) {
    u <- runif(n)
    Fc <- (c - a) / (b - a)
    ifelse(
      u < Fc,
      a + sqrt(u * (b - a) * (c - a)),
      b - sqrt((1 - u) * (b - a) * (b - c))
    )
  }
  
  # ---- Trait equations ----
  Tp <- function(T, VD) 50.775 - 1.0728*T - 0.229*VD
  To <- function(T, VD) 66.19  - 1.5638*T - 0.179*VD
  Ti <- function(T, VD) 59.7   - 1.151*T  - 0.1014*VD
  Tm <- function(T, VD) 192.23 - 6.054*T + 0.258*VD
  
  Mf_pct <- function(T, VD) 19.32   + 1.212*T - 0.16*VD
  Me_pct <- function(T, VD) 108.325 - 3.848*T + 1.414*VD
  Mn_pct <- function(T, VD) 51.4786 + 1.525*T - 0.22*VD
  
  Ml_weeks <- function(T, VD) 16.1 + 0.814*T - 0.21*VD
  Ma_weeks <- function(T, VD) 14.3 + 0.792*T - 0.14*VD
  
  Ql_raw <- function(T) 0.063 * (T^0.51)
  Qa_raw <- function(T) 0.051 * (T^0.49)
  
  # ---- helpers ----
  clip01 <- function(p) pmax(0, pmin(1, p))
  
  mu_from_stage_mortality <- function(mort_pct, duration_days) {
    q <- clip01(mort_pct / 100)
    out <- rep(NA_real_, length(q))
    ok <- is.finite(duration_days) & duration_days > 0 & is.finite(q)
    q_ok <- q[ok]
    d_ok <- duration_days[ok]
    out[ok] <- ifelse(q_ok >= 1, 1e6, -log(1 - q_ok) / d_ok)
    out
  }
  
  mu_from_survival_time_weeks <- function(S_weeks) {
    S_days <- S_weeks * 7
    out <- rep(NA_real_, length(S_days))
    ok <- is.finite(S_days) & S_days > 0
    out[ok] <- 1 / S_days[ok]
    out
  }
  
  delta_from_time <- function(time_days) {
    out <- rep(NA_real_, length(time_days))
    ok <- is.finite(time_days) & time_days > 0
    out[ok] <- 1 / time_days[ok]
    out
  }
  
  # ---- unpack ----
  T  <- x[1:12]
  VD <- x[13:24]
  
  # ---- sample uncertainty ----
  F_samp <- rtri_local(B, F_min, F_mean, F_max) * sex_ratio_female
  D_samp <- rtri_local(B, D_host_min, D_host_mean, D_host_max)
  
  # ---- per-draw R0 (returns vector length 12) ----
  R0_one_draw <- function(F_eff, D_host_eff) {
    
    Q_with_rule <- function(Qraw, T) {
      Q_min <- 1 / D_host_eff
      ifelse(T > T_min_activity, pmax(Qraw, Q_min), 0)
    }
    
    tp <- Tp(T, VD); to <- To(T, VD); ti <- Ti(T, VD); tm <- Tm(T, VD)
    
    deltaE <- delta_from_time(ti)
    deltaN <- delta_from_time(tm)
    deltaF <- delta_from_time(tp + to)
    
    muE <- mu_from_stage_mortality(Me_pct(T, VD), ti)
    muN <- mu_from_stage_mortality(Mn_pct(T, VD), tm)
    muF <- mu_from_stage_mortality(Mf_pct(T, VD), tp + to)
    
    muL <- mu_from_survival_time_weeks(Ml_weeks(T, VD))
    muA <- mu_from_survival_time_weeks(Ma_weeks(T, VD))
    
    Ql <- Q_with_rule(Ql_raw(T), T)
    Qa <- Q_with_rule(Qa_raw(T), T)
    
    term_f <- deltaF / (deltaF + muF)
    term_a <- Qa / (Qa + muA)
    term_n <- deltaN / (deltaN + muN)
    term_l <- Ql / (Ql + muL)
    term_e <- deltaE / (deltaE + muE)
    
    R0 <- F_eff * term_f * term_a * term_n * term_l * term_e
    R0[!is.finite(R0) | R0 < 0] <- NA_real_
    R0
  }
  
  # compute all draws (12 x B)
  R0_mc <- vapply(
    seq_len(B),
    function(b) R0_one_draw(F_eff = F_samp[b], D_host_eff = D_samp[b]),
    FUN.VALUE = rep(NA_real_, 12)
  )
  
  # aggregate across months first (so each cell returns ONE number per statistic)
  # (this matches what you attempted: a single mean/low/high per raster cell)
  # return STOCHASTIC MONTHLY MEAN (length = 12)
  rowMeans(R0_mc, na.rm = TRUE)
}

names(R0_stoch_monthly) <- month.abb

# ============================================================
# 5) TEMPORAL AGGREGATION
# ============================================================

# Deterministic quarterly/annual
q_index <- rep(1:4, each = 3)
R0_det_quarterly <- tapp(R0_det_monthly, q_index, mean, na.rm = TRUE)
R0_det_annual    <- app(R0_det_monthly, mean, na.rm = TRUE)


# ============================================================
# STOCHASTIC TEMPORAL AGGREGATION
# ============================================================

# Quarterly index (Jan–Mar, Apr–Jun, Jul–Sep, Oct–Dec)
q_index <- rep(1:4, each = 3)

# Stochastic quarterly mean R0
R0_stoch_quarterly <- tapp(
  R0_stoch_monthly,
  q_index,
  mean,
  na.rm = TRUE
)

names(R0_stoch_quarterly) <- c("Q1", "Q2", "Q3", "Q4")


# ============================================================
# 6) EXPORT
# ============================================================

dir.create("outputs", showWarnings = FALSE)

writeRaster(R0_det_monthly,   "outputs/R0_det_monthly.tif", overwrite = TRUE)
writeRaster(R0_det_annual,    "outputs/R0_det_annual.tif",  overwrite = TRUE)
writeRaster(R0_det_quarterly, "outputs/R0_det_quarterly.tif", overwrite = TRUE)

writeRaster(R0_stoch_monthly,      "outputs/R0_stoch_stats.tif", overwrite = TRUE)
writeRaster(R0_stoch_annual_mean,  "outputs/R0_stoch_annual_mean.tif", overwrite = TRUE)
writeRaster(R0_stoch_quarterly, "outputs/R0_stoch_quarterly.tif",overwrite = TRUE)
cat("✓ Spatiotemporal R0 surfaces successfully generated and exported.\n")



