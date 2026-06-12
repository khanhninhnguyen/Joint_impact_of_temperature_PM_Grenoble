# ==============================================================================
# run_microsim_from_RR.R
#
#   Survival microsimulation driven by the 3D temperature relative-risk array
#   from convert_temp_to_RR.R:
#     * simulate a synthetic population (ID, start age, grid),
#     * simulate a age–mortality baseline hazard,
#     * combine baseline daily risk with the temperature RR into a daily
#       p(death) and build the survival curve S(t) = prod(1 - p),
#     * draw one Uniform(0,1) per individual and read off the death day as the
#       first day where S(t) < U.
#
#   Every individual is followed until certain death at MAX_AGE (= 140), so the
#   population average life expectancy is well defined (no censoring). The
#   temperature RR series is padded with neutral RR = 1 beyond the array.
#
#   The daily survival curve depends only on (grid, start_age), so each of the
#   n_grids × n_start_ages distinct curves is computed once and every death day
#   is read off with findInterval().
#
# Input : risk_3d_example.rds  (output of convert_temp_to_RR.R)
# Output: microsim_death_ages.csv  (ID, start_age, grid, death_day, death_age)
# ==============================================================================

# ---- Configuration -----------------------------------------------------------

rr_file   <- "risk_3d_example.rds"        # input from convert_temp_to_RR.R
out_file  <- "microsim_death_ages.csv"    # individual-level output

N_POP     <- 100000L   # synthetic population size
MIN_AGE   <- 30L       # youngest start age (also baseline age-specific mortality)
MAX_START <- 90L       # oldest start age
MAX_AGE   <- 140L      # certain death at this age (removes censoring)
DAYS_PER_YEAR <- 365L  # day 366 of leap years dropped

# ---- annual baseline mortality hazard -----------------------
#   mu(age) = GM_A + GM_B * exp(GM_b * (age - MIN_AGE))    # Gompertz rate (mortality ~doubles every ~8 yr)

seed_pop <- 123L   # RNG seed for population age / grid assignment
seed_sim <- 1234L  # RNG seed for the single Uniform(0,1) death draw

cat("================================================================\n")
cat("run_microsim_from_RR.R — fast survival microsimulation (base R)\n")
cat("================================================================\n\n")

# ---- Load and reconstruct the 3D temperature RR array ------------------------
if (!file.exists(rr_file)) {
  stop("RR file not found: ", rr_file,
       "\nRun convert_temp_to_RR.R first to produce it.")
}

cat(sprintf("Loading RR array from '%s'...\n", rr_file))
obj <- readRDS(rr_file)
qi  <- readBin(obj, integer(), n = prod(attr(obj, "array_dim")),
               size = 2L, signed = TRUE, endian = "little")
qi[qi == attr(obj, "na_sentinel")] <- NA_integer_
risk_array <- array(exp(qi / attr(obj, "quant_scale")),
                    dim = attr(obj, "array_dim"))   # [grid, year, day_of_year]

n_grids <- dim(risk_array)[1]
n_years <- dim(risk_array)[2]
cat(sprintf("  Dimensions : [%d grids x %d years x %d day slots]\n",
            dim(risk_array)[1], dim(risk_array)[2], dim(risk_array)[3]))
cat(sprintf("  RR range   : [%.3f, %.3f]\n\n",
            min(risk_array, na.rm = TRUE), max(risk_array, na.rm = TRUE)))

# ---- Flatten temperature RR to a daily series per grid -----------------------
# First DAYS_PER_YEAR days of each year, flattened year-by-year, padded with
# neutral RR = 1 to cover the youngest person's span to certain death.
horizon_days <- (MAX_AGE - MIN_AGE) * DAYS_PER_YEAR   # youngest person's max span

env_by_grid <- vector("list", n_grids)
for (g in seq_len(n_grids)) {
  series <- as.vector(t(risk_array[g, , seq_len(DAYS_PER_YEAR)]))
  if (length(series) < horizon_days) {
    series <- c(series, rep(1, horizon_days - length(series)))     # neutral tail
  }
  env_by_grid[[g]] <- series[seq_len(horizon_days)]
}
cat(sprintf("Built daily temperature RR series per grid (length %d days each;\n",
            horizon_days))
cat(sprintf("  %d real RR years + %d neutral-RR tail years).\n\n",
            n_years, (MAX_AGE - MIN_AGE) - n_years))

# ---- Age–mortality baseline -------------------------------
# Daily baseline p(death) for ages MIN_AGE .. MAX_AGE; at MAX_AGE it is forced
# to 1.0 so every individual eventually dies.
ages_all   <- MIN_AGE:MAX_AGE
mu_annual  <- GM_A + GM_B * exp(GM_b * (ages_all - MIN_AGE))   # annual hazard
daily_base <- pmin(mu_annual / DAYS_PER_YEAR, 1.0)             # daily p(death)
daily_base[length(daily_base)] <- 1.0                          # MAX_AGE -> certain
names(daily_base) <- ages_all

cat("Synthetic baseline annual mortality (per 1,000):\n")
for (a in c(30L, 50L, 70L, 90L)) {
  cat(sprintf("  age %3d : %6.2f\n", a, 1000 * mu_annual[ages_all == a]))
}
cat("\n")

# ---- Simulate the population (ID, start age, grid) ---------------------------
# Start ages are drawn proportional to the survivorship l(a) implied by the
# baseline hazard; each individual is assigned a grid uniformly at random.
set.seed(seed_pop)

start_ages_pool <- MIN_AGE:MAX_START
mu_pool   <- GM_A + GM_B * exp(GM_b * (start_ages_pool - MIN_AGE))
surv_l    <- exp(-cumsum(mu_pool))                  # survivorship from MIN_AGE
age_w     <- surv_l / sum(surv_l)

start_age <- sample(start_ages_pool, size = N_POP, replace = TRUE, prob = age_w)
grid_id   <- sample(seq_len(n_grids),  size = N_POP, replace = TRUE)
ID        <- seq_len(N_POP)

cat(sprintf("Simulated population: %d individuals, ages %d-%d, %d grids.\n",
            N_POP, MIN_AGE, MAX_START, n_grids))
cat(sprintf("  Mean start age: %.1f\n\n", mean(start_age)))

# ---- Survival simulation -----------------------------------------------------
# For each distinct (grid, start_age) compute the daily survival curve once:
#   p(t)       = clip(baseline_daily(age at t) * temp_RR(grid, t), 0, 1)
#   negLogS(t) = -cumsum(log(1 - p(t)))
# then read off each death day as findInterval(-log(U), negLogS).
set.seed(seed_sim)
U <- runif(N_POP)
E <- -log(U)                          # Exponential(1) target cumulative hazard

death_day <- integer(N_POP)

combos <- unique(data.frame(grid = grid_id, age = start_age))
cat(sprintf("Computing %d distinct (grid x start_age) survival curves...\n",
            nrow(combos)))

for (k in seq_len(nrow(combos))) {
  g <- combos$grid[k]
  a <- combos$age[k]

  n_days_c <- (MAX_AGE - a) * DAYS_PER_YEAR        # days until certain death

  # Daily baseline aligned to this person's age track, and the grid RR.
  base_idx  <- (a - MIN_AGE) + rep(seq_len(MAX_AGE - a) - 1L, each = DAYS_PER_YEAR) + 1L
  base_days <- daily_base[base_idx][seq_len(n_days_c)]
  env_days  <- env_by_grid[[g]][seq_len(n_days_c)]

  p_death <- pmin(base_days * env_days, 1.0)
  negLogS <- -cumsum(log1p(-p_death))
  negLogS[!is.finite(negLogS)] <- Inf              # guarantee a crossing

  sel <- which(grid_id == g & start_age == a)
  death_day[sel] <- findInterval(E[sel], negLogS)
}

death_age <- start_age + death_day / DAYS_PER_YEAR

# ---- Output and summary ------------------------------------------------------
out <- data.frame(
  ID        = ID,
  start_age = start_age,
  grid      = grid_id,
  death_day = death_day,
  death_age = round(death_age, 4)
)
write.csv(out, out_file, row.names = FALSE)
cat(sprintf("\nSaved individual results -> '%s'  (%d rows)\n\n", out_file, nrow(out)))

deaths_first_year <- sum(death_day < DAYS_PER_YEAR)
le_avg            <- mean(death_age)

cat("=== RESULTS ===\n")
cat(sprintf("  Population size                 : %d\n", N_POP))
cat(sprintf("  Deaths in first simulation year : %d  (%.3f%% of population)\n",
            deaths_first_year, 100 * deaths_first_year / N_POP))
cat(sprintf("  Average life expectancy (age)   : %.2f years\n", le_avg))
cat(sprintf("  Mean remaining years from start : %.2f years\n",
            mean(death_age - start_age)))
cat("================================================================\n")
