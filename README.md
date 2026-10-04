# bruce <img src="https://img.shields.io/badge/status-experimental-orange" align="right"/>

**B**ayesian **R**io **U**rban **C**limate **E**stimator, generalised to any city.

`bruce` interpolates daily urban air temperature by fusing observation
networks of different quality — official stations, airports, low-cost IoT
sensors, citizen weather stations — with an optional gridded background
(e.g. ERA5-Land). It was developed for Rio de Janeiro's heat surveillance and
is written so that any city with comparable data can use it: everything
city-specific enters through a small data contract.

## Model

For an observation of network *k* at station *i*, location *s*, day *t*:

η = β₀ + δₖ + x(s)ᵀβ + uᵢ + dₜ + w(s, t)

* **δₖ** — offset of a biased network (estimated only for networks declared `bias = TRUE`)
* **x(s)ᵀβ** — spatial covariates (elevation, urban fraction, …)
* **uᵢ** — iid station effect (persistent siting bias)
* **dₜ** — city-wide daily AR(1) term
* **w(s, t)** — Matérn SPDE field with AR(1) dependence between days

Each network (and the background) has its own Gaussian likelihood with a
fixed or PC-prior observation standard deviation. Windows (one month is the
tested scale) are fitted independently with R-INLA.

## Installation

```r
install.packages("INLA", repos = c(getOption("repos"),
                 INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
remotes::install_github("asgouveiaa/bruce")
```

## Data contract

| object | columns |
|---|---|
| `obs` | `station_id`, `date`, `value` (one value per station and day) |
| `stations` | `station_id`, `source`, and `x`/`y` (projected) or `lon`/`lat` |
| `sources` | `source`, `bias` (`FALSE` = reference), `obs_sd` (fixed, or `NA` to estimate) |
| `background` *(optional)* | `date`, `value`, `x`/`y` or `lon`/`lat` of cell centres |
| covariates | a `terra::SpatRaster`, one named layer per covariate |

At least one reference network is required.

## Workflow

```r
library(bruce)
sim  <- bruce_simulate(seed = 1)          # or your own city's data

d    <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs,
                   background = sim$background)
d    <- bruce_covariates(d, sim$covariates)
mesh <- bruce_mesh(d)
fit  <- bruce_fit(d, mesh, bruce_spec(), keep_config = TRUE)
bruce_check(fit)                          # every attempt, collapse diagnostics

pred <- bruce_predict(fit, n_samples = 200, keep_samples = TRUE)
expo <- bruce_zonal(pred, zones, id_col = "zone", weights = population)
cv   <- bruce_cv(d, mesh, scheme = "lso", checkpoint_dir = "cv_lso")
```

## Degenerate fits

With clustered networks the day-to-day correlation of the field can drift to
1 and its variance explode. `bruce_fit()` checks every fit and, by default,
restarts the mode search from a reset point and, if that also collapses,
refits with the correlation fixed (`fallback_group_rho`, default 0.55). All
attempts are kept in `fit$attempts`; `fit$status` records the outcome. See
`?bruce_spec` (`on_collapse`).

## Status

Version 0.1: daily temperature, independent windows. Planned: sub-daily
resolution, humidity and heat index, bivariate fields.

`inst/examples/rio_adapter.R` maps the original Rio de Janeiro inputs onto
the contract.
