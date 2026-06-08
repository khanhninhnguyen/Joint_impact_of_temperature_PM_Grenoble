# ==============================================================================
# convert_temp_to_RR.R
#
#   1. Load the fitted DLNM model (ANAS_with_UHI, 2000-2016) and its crossbasis.
#   2. Simulate daily temperature for several grid points.
#   3. Apply the model's point-estimate coefficients and assemble a 3D
#      cumulative-lag risk array [n_grids × n_years × n_day_max].
#   4. Quantize log(RR) to int16 and save as an xz-compressed RDS.
#
# Reconstruct the saved output:
#   obj        <- readRDS("risk_3d_example.rds")
#   qi         <- readBin(obj, integer(), n = prod(attr(obj, "array_dim")),
#                         size = 2L, signed = TRUE, endian = "little")
#   qi[qi == attr(obj, "na_sentinel")] <- NA_integer_
#   risk_array <- array(exp(qi / attr(obj, "quant_scale")),
#                       dim = attr(obj, "array_dim"))
#
# Requires: install.packages(c("dlnm", "parallel"))
# ==============================================================================

library(dlnm)
library(parallel)

# ---- Configuration -----------------------------------------------------------

# ---- Model files -------------------------------------------------------------
ERF_DIR    <- "/Users/knguyen04/Documents/PACC-MACS/Results/ERFs"
model_name <- "ANAS_with_UHI"
period_fit <- "2000-2016"

# ---- Simulation settings -----------------------------------------------------
n_grids    <- 3L    # number of spatial grid points to simulate
n_years    <- 100L  # simulation length (years)
year_start <- 2001L # first calendar year

# Per-grid mean temperature (°C); keep within the model's training range.
grid_mean_temp <- c(10.0, 14.0, 18.0)  # length must equal n_grids
temp_amplitude <- 10.0   # seasonal swing (°C)
temp_noise_sd  <- 2.0    # day-to-day random noise (°C)

seed_temp <- 42L   # RNG seed for temperature simulation

out_file <- "risk_3d_example.rds"   # output (quantized int16, xz-compressed)

stopifnot(length(grid_mean_temp) == n_grids, n_years >= 1L, n_grids >= 1L)

cat("================================================================\n")
cat("convert_temp_to_RR.R\n")
cat("  model      :", model_name, "/", period_fit, "\n")
cat("  n_grids    :", n_grids,  "\n")
cat("  n_years    :", n_years,  "\n")
cat("  out_file   :", out_file, "\n")
cat("================================================================\n\n")

# ---- Load fitted model and crossbasis ----------------------------------------
mod_file <- file.path(ERF_DIR, paste0("fitted_ERF_82_", period_fit, model_name, ".rds"))
cb_file  <- file.path(ERF_DIR, paste0("crossbasis_82_", period_fit, model_name, ".rds"))

if (!file.exists(mod_file)) stop("Model RDS not found:\n  ", mod_file)
if (!file.exists(cb_file))  stop("Crossbasis RDS not found:\n  ", cb_file)

cat("Loading model and crossbasis...\n")
mod <- readRDS(mod_file)
cb  <- readRDS(cb_file)

max_lag <- attr(cb, "lag")[2]

# Centering temperature = temperature minimising cumulative RR over the
# training data range.
pred0          <- crosspred(cb, mod, at = mod$data$temp, cumul = TRUE)
centering_temp <- as.numeric(names(pred0$allRRfit)[which.min(pred0$allRRfit)])
rm(pred0)

cat(sprintf("  Model loaded  : %s / %s\n", model_name, period_fit))
cat(sprintf("  max_lag       : %d days\n", max_lag))
cat(sprintf("  centering_temp: %.2f °C\n\n", centering_temp))

# ---- Simulate temperature (3D array: n_grids × n_years × 366) -----------------
# temp_3d[g, y, d] = temperature at grid g, year y, day-of-year d.
# Day 366 slot is NA for non-leap years.
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
  }
}

cat(sprintf("  Dimensions : [%d grids × %d years × %d day slots]\n",
            n_grids, n_years, n_day_max))
cat(sprintf("  Temp range : [%.1f, %.1f] °C\n",
            min(temp_3d, na.rm = TRUE), max(temp_3d, na.rm = TRUE)))
cat(sprintf("  Leap years : %d / %d\n\n", sum(leap_flags), n_years))

# ---- Cumulative risk via manual lag reconstruction ---------------------------
# For each day t, the daily RR is the product over lags 0 … max_lag of the
# RR contributed by the temperature on day t-l at lag l.
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

# ---- Apply ERF to one grid → [n_years × n_day_max] risk matrix ----------------
# Transpose to [n_day_max × n_years] so column-major traversal of which() yields
# year-by-year, day-within-year order. NAs are excluded and mapped back after.
apply_erf_to_grid <- function(g, n_yrs, n_dmax) {
  grid_temps   <- temp_3d[g, , ]      # [n_years × n_day_max]
  grid_temps_t <- t(grid_temps)       # [n_day_max × n_years]

  idx <- which(!is.na(grid_temps_t), arr.ind = TRUE)

  if (nrow(idx) == 0L) return(matrix(NA_real_, nrow = n_yrs, ncol = n_dmax))

  d_idx     <- idx[, 1L]
  y_idx     <- idx[, 2L]
  temp_vals <- grid_temps_t[idx]

  risks <- compute_grid_risk(temp_vals, cb, mod, centering_temp, max_lag)

  risk_mat <- matrix(NA_real_, nrow = n_yrs, ncol = n_dmax)
  risk_mat[cbind(y_idx, d_idx)] <- risks
  risk_mat
}

# ---- Apply model to all grids ------------------------------------------------
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

# ---- Quantize and save -------------------------------------------------------
# log(RR) × QUANT_SCALE → int16 raw bytes, little-endian. NA cells use -32768.
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

# ---- Summary -----------------------------------------------------------------
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
