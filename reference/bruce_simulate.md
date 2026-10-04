# Simulate a synthetic city

Generates data from the BRUCE data-generating process on an invented
city, for testing the pipeline and checking parameter recovery: two
covariates (`elev`, a pair of hills; `urb`, an urban core), a reference
network, a biased low-cost network, an optional coarse background
product, a Matern field with AR(1) dependence between days and a
city-wide daily AR(1) term.

## Usage

``` r
bruce_simulate(
  n_reference = 12,
  n_biased = 40,
  n_days = 10,
  extent = 30000,
  range = 8000,
  sigma = 1,
  group_rho = 0.6,
  intercept = 25,
  beta_elev = -0.0065,
  beta_urb = 1.5,
  day_sd = 1.5,
  day_rho = 0.7,
  bias = 1.5,
  sd_reference = 0.3,
  sd_biased = 1,
  sd_station = 0.5,
  background = TRUE,
  background_res = 9000,
  background_sd = 0.7,
  res = 250,
  crs = 31983,
  origin = c(650000, 7450000),
  clustered = FALSE,
  seed = NULL
)
```

## Arguments

- n_reference, n_biased:

  Number of reference and low-cost stations.

- n_days:

  Number of days.

- extent:

  Side of the square domain (metres).

- range, sigma:

  Matern range (metres) and standard deviation of the field (smoothness
  1).

- group_rho:

  Day-to-day correlation of the field.

- intercept, beta_elev, beta_urb:

  Mean level and covariate effects (`beta_elev` per metre of elevation,
  `beta_urb` per unit urban fraction).

- day_sd, day_rho:

  Standard deviation and autocorrelation of the daily term.

- bias:

  Offset of the low-cost network.

- sd_reference, sd_biased:

  Observation noise of each network.

- sd_station:

  Standard deviation of persistent per-station offsets in the low-cost
  network.

- background:

  Simulate a coarse background product.

- background_res, background_sd:

  Its cell size (metres) and noise.

- res:

  Covariate raster resolution (metres).

- crs:

  EPSG code used for the coordinates.

- origin:

  Lower-left corner of the domain in `crs`.

- clustered:

  If `TRUE`, low-cost stations concentrate in the urban core, mimicking
  the deployment pattern that drove field collapse in the original
  application.

- seed:

  Random seed.

## Value

A list with `obs`, `stations`, `sources`, `background`, `covariates`
(SpatRaster), `domain` (sf polygon), `crs`, `truth` (SpatRaster, one
layer per day, latent temperature) and `params`. The first four plug
straight into
[`bruce_data()`](https://asgouveiaa.github.io/bruce/reference/bruce_data.md).
