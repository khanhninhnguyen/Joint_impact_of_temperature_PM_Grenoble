# ==============================================================================
# SCRIPT: convert_temp_to_RR.R
#
# PURPOSE
#   Self-contained example that:
#     1. Loads the real fitted DLNM model (ANAS_with_UHI, 2000-2016) and its
#        crossbasis from local RDS files.
#     2. Simulates 100 years of synthetic daily temperature for multiple grid
#        points (seasonal sine-wave + noise; NOT the real-year resampling done
#        by 02_prepare_temp_series.R — this is illustrative only).
#     3. Applies the model's point-estimate coefficients to every grid point
#        and assembles a 3D cumulative-lag risk array
#        [n_grids × n_years × n_day_max] — identical Phase 03 logic
#        (03_compute_risk_3d.R / 00_compute_risk_3d_raw.R).
#     4. Quantizes log(RR) to int16 and saves as an xz-compressed RDS — same
#        storage format as Phase 03 output.
#
# RECONSTRUCT SAVED OUTPUT DOWNSTREAM
#   obj        <- readRDS("risk_3d_example.rds")
#   qi         <- readBin(obj, integer(), n = prod(attr(obj, "array_dim")),
#                         size = 2L, signed = TRUE, endian = "little")
#   qi[qi == attr(obj, "na_sentinel")] <- NA_integer_
#   risk_array <- array(exp(qi / attr(obj, "quant_scale")),
#                       dim = attr(obj, "array_dim"))
#   # risk_array[grid, year, day_of_year]
#
# REQUIRED PACKAGES
#   install.packages(c("dlnm", "parallel"))
# ==============================================================================

library(dlnm)
library(parallel)

# ==============================================================================
# SECTION 0 — CONFIGURATION
# Edit these values to match your use case.
# ==============================================================================

# ---- Model files (real fitted model on temp_with_uhi, 2000-2016) -------------
ERF_DIR    <- "/Users/knguyen04/Documents/PACC-MACS/Results/ERFs"
model_name <- "ANAS_with_UHI"
period_fit <- "2000-2016"

# ---- Simulation settings -----------------------------------------------------
n_grids    <- 3L    # number of spatial grid points to simulate
n_years    <- 100L  # simulation length (years)
year_start <- 2001L # first calendar year

# Per-grid mean temperature (°C); values should stay within the model's
# training range (~temp_with_uhi observed in French cities 2000-2016).
grid_mean_temp <- c(10.0, 14.0, 18.0)  # length must equal n_grids
temp_amplitude <- 10.0   # seasonal swing (°C)
temp_noise_sd  <- 2.0    # day-to-day random noise (°C)

seed_temp <- 42L   # RNG seed for temperature simulation

out_file <- "risk_3d_example.rds"   # output (quantized int16, xz-compressed)

# ==============================================================================
# VALIDATION
# ==============================================================================
stopifnot(length(grid_mean_temp) == n_grids, n_years >= 1L, n_grids >= 1L)

cat("================================================================\n")
cat("convert_temp_to_RR.R  — Phase 03 worked example (point estimate)\n")
cat("  model      :", model_name, "/", period_fit, "\n")
cat("  n_grids    :", n_grids,  "\n")
cat("  n_years    :", n_years,  "\n")
cat("  out_file   :", out_file, "\n")
cat("================================================================\n\n")

# ==============================================================================
# SECTION 1 — LOAD FITTED MODEL AND CROSSBASIS
#
# Mirrors 00_compute_risk_3d_raw.R:130-137.
# ==============================================================================
mod_file <- file.path(ERF_DIR, paste0("fitted_ERF_82_", period_fit, model_name, ".rds"))
cb_file  <- file.path(ERF_DIR, paste0("crossbasis_82_", period_fit, model_name, ".rds"))

if (!file.exists(mod_file)) stop("Model RDS not found:\n  ", mod_file)
if (!file.exists(cb_file))  stop("Crossbasis RDS not found:\n  ", cb_file)

cat("Loading model and crossbasis...\n")
mod <- readRDS(mod_file)
cb  <- readRDS(cb_file)

max_lag <- attr(cb, "lag")[2]

# Centering temperature = temperature minimising cumulative RR in the
# training data range.  Identical to 00_compute_risk_3d_raw.R:136-137.
pred0          <- crosspred(cb, mod, at = mod$data$temp, cumul = TRUE)
centering_temp <- as.numeric(names(pred0$allRRfit)[which.min(pred0$allRRfit)])
rm(pred0)

cat(sprintf("  Model loaded  : %s / %s\n", model_name, period_fit))
cat(sprintf("  max_lag       : %d days\n", max_lag))
cat(sprintf("  centering_temp: %.2f °C\n\n", centering_temp))

# ==============================================================================
# SECTION 2 — SIMULATE TEMPERATURE (3D array: n_grids × n_years × 366)
#
# Layout: temp_3d[g, y, d] = temperature at grid g, year y, day-of-year d.
# Day 366 slot is NA for non-leap years (same convention as Phase 02/03).
#
# NOTE: Phase 02 (02_prepare_temp_series.R) builds this array by resampling
# real historical years with replacement plus a scenario trend offset.
# Here we use a simple sine-wave generator for illustration only.
# ==============================================================================
is_leap_year <- function(y) (y %% 4 == 0 & y %% 100 != 0) | (y %% 400 == 0)

years_vec  <- year_start + seq_len(n_years) - 1L
leap_flags <- is_leap_year(years_vec)
n_day_max  <- 366L

cat("Simulating temperature array...\n")
set.seed(seed_temp)

temp_3d <- array(NA_real_, dim = c(n_grids, n_years, n_day_max))

for (g in seq_len(n_grids)) {
  for (yi in seq_len(n_years)) {
    n_days <- if (leap_flags[yi]) 366L else 365L
    doy    <- seq_len(n_days)
    temp_day <- grid_mean_temp[g] +
                temp_amplitude * sin(2 * pi * (doy - 80) / n_days) +
                rnorm(n_days, sd = temp_noise_sd)
    temp_3d[g, yi, seq_len(n_days)] <- temp_day
    # Slot [g, yi, 366] stays NA for non-leap years
  }
}

cat(sprintf("  Dimensions : [%d grids × %d years × %d day slots]\n",
            n_grids, n_years, n_day_max))
cat(sprintf("  Temp range : [%.1f, %.1f] °C\n",
            min(temp_3d, na.rm = TRUE), max(temp_3d, na.rm = TRUE)))
cat(sprintf("  Leap years : %d / %d\n\n", sum(leap_flags), n_years))

# ==============================================================================
# SECTION 3 — HELPER: CUMULATIVE RISK VIA MANUAL LAG RECONSTRUCTION
#
# Copied verbatim from 03_compute_risk_3d.R:312-331.
#
# For each day t:
#   - lag 0 contribution : RR from today's temperature at lag 0
#   - lag l contribution : RR from temperature on day t-l at lag l
# Final daily RR = product of contributions across all lags 0 … max_lag.
# ==============================================================================
compute_grid_risk <- function(temp_series, cb, mod, cen_temp, max_lag_val) {
  pred_exact   <- crosspred(cb, mod, at = temp_series, cen = cen_temp)
  time_indices <- match(temp_series, pred_exact$predvar)
  RR_all       <- pred_exact$matRRfit[time_indices, , drop = FALSE]

  n_days <- length(temp_series)
  RR_mat <- matrix(1, nrow = n_days, ncol = max_lag_val + 1L)

  for (l in 0:max_lag_val) {
    col_idx <- l + 1L
    if (l == 0L) {
      RR_mat[, col_idx] <- RR_all[, col_idx]
    } else {
      idx_curr                  <- (l + 1L):n_days
      idx_exp                   <- 1L:(n_days - l)
      RR_mat[idx_curr, col_idx] <- RR_all[idx_exp, col_idx]
    }
  }
  apply(RR_mat, 1L, prod)
}

# ==============================================================================
# SECTION 4 — APPLY ERF TO ONE GRID → [n_years × n_day_max] RISK MATRIX
#
# Mirrors apply_erf_to_grid() in 00_compute_risk_3d_raw.R:172-193.
#
# Key ordering trick: transpose temp_3d[g, , ] to [n_day_max × n_years] so
# column-major traversal of which() gives year-by-year, day-within-year order —
# the correct temporal sequence for lag reconstruction.
# NAs (day 366 in non-leap years) are excluded and mapped back afterwards.
# ==============================================================================
apply_erf_to_grid <- function(g, n_yrs, n_dmax) {
  grid_temps   <- temp_3d[g, , ]      # [n_years × n_day_max]
  grid_temps_t <- t(grid_temps)       # [n_day_max × n_years]

  idx <- which(!is.na(grid_temps_t), arr.ind = TRUE)
  # idx[, 1] = day-of-year index, idx[, 2] = year index

  if (nrow(idx) == 0L) return(matrix(NA_real_, nrow = n_yrs, ncol = n_dmax))

  d_idx     <- idx[, 1L]
  y_idx     <- idx[, 2L]
  temp_vals <- grid_temps_t[idx]

  risks <- compute_grid_risk(temp_vals, cb, mod, centering_temp, max_lag)

  risk_mat <- matrix(NA_real_, nrow = n_yrs, ncol = n_dmax)
  risk_mat[cbind(y_idx, d_idx)] <- risks
  risk_mat
}

# ==============================================================================
# SECTION 5 — APPLY MODEL TO ALL GRIDS (Phase 00/03 core, single ERF)
# ==============================================================================
n_cores_avail <- max(1L, detectCores(logical = FALSE) - 1L)
use_parallel  <- n_cores_avail > 1L && .Platform$OS.type == "unix"

cat(sprintf("Applying model to %d grids (%s)...\n", n_grids,
            if (use_parallel) sprintf("parallel, %d cores", n_cores_avail)
            else "sequential"))

t0 <- proc.time()

if (use_parallel) {
  grid_risks <- mclapply(seq_len(n_grids), apply_erf_to_grid,
                         n_yrs = n_years, n_dmax = n_day_max,
                         mc.cores = n_cores_avail)
} else {
  grid_risks <- lapply(seq_len(n_grids), apply_erf_to_grid,
                       n_yrs = n_years, n_dmax = n_day_max)
}

elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Done in %.1f s\n\n", elapsed))

# Assemble 3D array [n_grids × n_years × n_day_max]
risk_3d <- array(NA_real_, dim = c(n_grids, n_years, n_day_max))
for (g in seq_len(n_grids)) risk_3d[g, , ] <- grid_risks[[g]]
rm(grid_risks)

cat(sprintf("risk_3d dims     : [%d × %d × %d]\n", dim(risk_3d)[1], dim(risk_3d)[2], dim(risk_3d)[3]))
cat(sprintf("RR range         : [%.3f, %.3f]\n\n",
            min(risk_3d, na.rm = TRUE), max(risk_3d, na.rm = TRUE)))

# ==============================================================================
# SECTION 6 — QUANTIZE AND SAVE (Phase 03 output format)
#
# log(RR) × QUANT_SCALE → int16 raw bytes, little-endian.
# NA cells (day 366 in non-leap years) use sentinel -32768.
# ==============================================================================
QUANT_SCALE <- 10000L
NA_SENTINEL <- -32768L
INT16_MIN   <- -32767L
INT16_MAX   <-  32767L

cat(sprintf("Quantizing to int16 and saving -> '%s'...\n", out_file))

qi      <- as.integer(round(log(risk_3d) * QUANT_SCALE))
na_mask <- is.na(qi)
if (any(!na_mask)) {
  qi[!na_mask & qi >  INT16_MAX] <- INT16_MAX
  qi[!na_mask & qi <  INT16_MIN] <- INT16_MIN
}
qi[na_mask] <- NA_SENTINEL

con       <- rawConnection(raw(0L), "wb")
writeBin(qi, con, size = 2L, endian = "little")
raw_bytes <- rawConnectionValue(con)
close(con)
rm(qi, con)

attr(raw_bytes, "array_dim")      <- dim(risk_3d)
attr(raw_bytes, "quant_scale")    <- QUANT_SCALE
attr(raw_bytes, "na_sentinel")    <- NA_SENTINEL
attr(raw_bytes, "encoding")       <- "log_rr_int16_le"
attr(raw_bytes, "dim_names")      <- c("grid", "year", "day_of_year")
attr(raw_bytes, "years")          <- years_vec
attr(raw_bytes, "model_name")     <- model_name
attr(raw_bytes, "period_fit")     <- period_fit
attr(raw_bytes, "centering_temp") <- centering_temp
attr(raw_bytes, "max_lag")        <- max_lag

saveRDS(raw_bytes, out_file, compress = "xz")
cat(sprintf("  Saved: %.1f KB\n\n", file.info(out_file)$size / 1024))

# ==============================================================================
# SECTION 7 — SUMMARY
# ==============================================================================
cat("=== Summary ===\n")
cat(sprintf("  model          : %s / %s\n", model_name, period_fit))
cat(sprintf("  centering_temp : %.2f °C\n", centering_temp))
cat(sprintf("  max_lag        : %d days\n", max_lag))
cat(sprintf("  Output dims    : [%d grids × %d years × %d days]\n",
            n_grids, n_years, n_day_max))

cat("\nObjects in workspace:\n")
cat("  risk_3d       — 3D RR array [n_grids × n_years × n_day_max]\n")
cat("  temp_3d       — simulated temperature array [n_grids × n_years × 366]\n")
cat("  cb, mod       — fitted DLNM crossbasis and GLM objects\n")
cat("  centering_temp— temperature minimising cumulative RR\n")

cat(sprintf("\nSaved -> '%s'\n", out_file))
cat("Reconstruct with:\n")
cat("  obj  <- readRDS('", out_file, "')\n", sep = "")

cat("\n================================================================\n")
cat("DONE — convert_temp_to_RR.R\n")
cat("================================================================\n")
