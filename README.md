# Joint Impact of Temperature and PM — Grenoble

Tools for translating a fitted temperature–mortality exposure–response function (ERF)
into a daily relative-risk (RR) surface, and for running a survival microsimulation
driven by that surface.

The pipeline has two stages:

1. **`convert_temp_to_RR.R`** — turn a fitted DLNM model into a 3D RR array.
2. **`run_microsim_from_RR.R`** — run a population survival microsimulation on that array.

```
fitted DLNM model ──► convert_temp_to_RR.R ──► risk_3d_example.rds ──► run_microsim_from_RR.R ──► microsim_death_ages.csv
```

## Requirements

- R (base R is enough for the microsimulation)
- R packages: `dlnm`, `parallel`

```r
install.packages(c("dlnm", "parallel"))
```

The parallel path uses `mclapply` and therefore only kicks in on Unix-like
systems (macOS, Linux). On Windows it falls back to sequential execution.

## Stage 1 — `convert_temp_to_RR.R`

Loads a fitted distributed-lag non-linear model (DLNM) and its crossbasis,
simulates daily temperature for a handful of grid points, applies the model's
point-estimate coefficients, and assembles a cumulative-lag RR array of shape
`[n_grids × n_years × 366]` (the 366th day-of-year slot is `NA` in non-leap years).

**Input.** Two RDS files expected under `ERF_DIR` (set near the top of the script,
currently `…/PACC-MACS/Results/ERFs`):

- `fitted_ERF_82_<period_fit><model_name>.rds` — the fitted GLM
- `crossbasis_82_<period_fit><model_name>.rds` — the matching DLNM crossbasis

Defaults: `model_name = "ANAS_with_UHI"`, `period_fit = "2000-2016"`.

**Key settings** (top of the script):

| Setting | Default | Meaning |
|---|---|---|
| `n_grids` | `3` | number of spatial grid points |
| `n_years` | `100` | simulation length in years |
| `year_start` | `2001` | first calendar year |
| `grid_mean_temp` | `c(10, 14, 18)` | per-grid mean temperature (°C); length must equal `n_grids` |
| `temp_amplitude` | `10` | seasonal swing (°C) |
| `temp_noise_sd` | `2` | day-to-day random noise (°C) |
| `seed_temp` | `42` | RNG seed for temperature simulation |

The centering temperature (the temperature that minimises cumulative RR over the
training range) is derived automatically from the model.

**Output.** `risk_3d_example.rds` — `log(RR)` quantized to `int16`, stored as
xz-compressed raw bytes. NA cells use the sentinel `-32768`. Metadata
(array dims, quant scale, model name, centering temperature, max lag, etc.) is
attached as attributes.

Reconstruct the RR array:

```r
obj        <- readRDS("risk_3d_example.rds")
qi         <- readBin(obj, integer(), n = prod(attr(obj, "array_dim")),
                      size = 2L, signed = TRUE, endian = "little")
qi[qi == attr(obj, "na_sentinel")] <- NA_integer_
risk_array <- array(exp(qi / attr(obj, "quant_scale")),
                    dim = attr(obj, "array_dim"))   # [grid, year, day_of_year]
```

Run it:

```sh
Rscript convert_temp_to_RR.R
```

## Stage 2 — `run_microsim_from_RR.R`

A fast base-R survival microsimulation driven by the RR array:

- simulates a synthetic population (`ID`, start age, grid);
- builds a Gompertz–Makeham age–mortality baseline hazard;
- combines the daily baseline risk with the temperature RR into a daily
  `p(death)` and forms the survival curve `S(t) = ∏ (1 − p)`;
- draws one `Uniform(0,1)` per individual and reads off the death day as the
  first day where `S(t)` drops below the draw.

Everyone is followed to certain death at `MAX_AGE` (140), so there is no
censoring and population life expectancy is well defined. The temperature RR
series is padded with neutral RR = 1 beyond the array horizon. The daily
survival curve depends only on `(grid, start_age)`, so each distinct curve is
computed once and reused via `findInterval()`.

**Key settings** (top of the script):

| Setting | Default | Meaning |
|---|---|---|
| `N_POP` | `100000` | synthetic population size |
| `MIN_AGE` | `30` | youngest start age (and baseline-hazard anchor) |
| `MAX_START` | `90` | oldest start age |
| `MAX_AGE` | `140` | certain death age (removes censoring) |
| `DAYS_PER_YEAR` | `365` | leap-year day 366 is dropped |
| `GM_A`, `GM_B`, `GM_b` | `0.0005`, `0.0008`, `0.085` | Gompertz–Makeham hazard parameters |
| `seed_pop`, `seed_sim` | `123`, `1234` | RNG seeds for population and death draw |

**Input.** `risk_3d_example.rds` (from stage 1).

**Output.** `microsim_death_ages.csv` with one row per individual:
`ID, start_age, grid, death_day, death_age`. A summary (population size, deaths
in the first simulation year, average life expectancy) is printed to the console.

Run it:

```sh
Rscript run_microsim_from_RR.R
```

## Files

| File | Description |
|---|---|
| `convert_temp_to_RR.R` | Stage 1 — DLNM model → 3D RR array |
| `run_microsim_from_RR.R` | Stage 2 — survival microsimulation on the RR array |
| `risk_3d_example.rds` | Quantized 3D RR array (output of stage 1) |
| `microsim_death_ages.csv` | Individual-level death ages (output of stage 2) |
