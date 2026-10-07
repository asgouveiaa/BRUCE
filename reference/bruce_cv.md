# Cross-validation by station or by region

Refits the model with stations held out and scores predictions at the
held-out reference stations.

## Usage

``` r
bruce_cv(
  d,
  mesh = NULL,
  spec = NULL,
  window = NULL,
  scheme = c("lso", "lro"),
  regions = NULL,
  region_col = NULL,
  test_sources = NULL,
  n_samples = 200,
  level = 0.9,
  scale = c("observation", "field"),
  noise_sd = 0,
  checkpoint_dir = NULL,
  theta_init = NULL,
  fix_theta = FALSE,
  fold_ids = NULL,
  seed = 1,
  verbose = TRUE
)
```

## Arguments

- d:

  A `bruce_fit` (recommended) or a `bruce_data`. Given a fit, its data,
  mesh, specification and window are reused (each can still be
  overridden), and every fold starts the hyperparameter search from the
  fit's mode (`theta_init = fit`): the result is the same as starting
  from scratch, only faster. If the fit ended with the day-to-day
  correlation fixed after a collapse, the folds keep it fixed at the
  same value.

- mesh, spec, window:

  As in
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md).
  `mesh` is required when `d` is a `bruce_data`; `spec` defaults to
  [`bruce_spec()`](https://asgouveiaa.github.io/bruce/reference/bruce_spec.md).

- scheme:

  `"lso"` (leave-station-out: one fold per test station, removing only
  that station) or `"lro"` (leave-region-out: one fold per region,
  removing every station of every network in the region; the background
  stays).

- regions:

  For `"lro"`: polygons (`sf`) with an id column.

- region_col:

  Name of the id column in `regions`.

- test_sources:

  Networks scored at held-out locations; default all reference networks
  (`bias = FALSE`).

- n_samples:

  Posterior samples per fold. With `0`, scores use posterior means and
  CRPS/coverage are not computed.

- level:

  Nominal level of the predictive interval for coverage.

- scale:

  What the predictive distribution describes, and therefore what RMSE,
  CRPS and coverage evaluate:

  - `"observation"` (default): the reading of a new station of the test
    network at the held-out location: the spatial surface plus the
    network's offset, a new station effect (drawn with the estimated
    station sd) and the network's measurement noise. This is the
    distribution the held-out readings come from, so coverage near
    `level` means calibrated intervals.

  - `"field"`: the reference-scale surface only, as mapped by
    [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)
    and as scored in the original Rio de Janeiro application. Readings
    vary around it by station effect and noise, so its intervals are
    expected to cover held-out readings less often than `level`. With
    `"observation"`, field-scale coverage and CRPS are also reported
    (`coverage_field`, `crps_field`).

- noise_sd:

  With `scale = "field"` only: extra Gaussian noise added to the samples
  before scoring (kept for compatibility; `scale = "observation"`
  accounts for noise and station effects properly).

- checkpoint_dir:

  Optional directory; each finished fold is saved there and skipped on
  rerun, so long runs can be resumed.

- theta_init:

  Passed to
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md)
  (a `bruce_fit` or a vector). Starting each fold from the mode of the
  full-data fit usually saves most of the mode search, since removing a
  few stations barely moves it; this is the default when `d` is a fit.

- fix_theta:

  Keep the hyperparameters fixed at `theta_init` in every fold (see
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md)):
  each fold becomes a single Gaussian approximation, much faster, at the
  cost of ignoring how hyperparameters would change without the held-out
  data. With `theta_init = NULL` the full data are fitted once first and
  their mode is used.

- fold_ids:

  Run only these folds: station ids for `"lso"`, region ids for `"lro"`.
  Useful for quick checks on a subset.

- seed:

  Seed for posterior sampling (`NULL` for parallel, non-reproducible
  sampling; see
  [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)).

- verbose:

  Print progress.

## Value

An object of class `bruce_cv` with per-observation `predictions`,
per-fold `folds` (scores and fit status) and `overall` scores, all as
tibbles, plus the `scale` that was scored.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- bruce_fit(d, mesh, spec, keep_config = TRUE)
cv  <- bruce_cv(fit)                          # leave-station-out, warm start
cv  <- bruce_cv(fit, fold_ids = c("A621", "SBRJ"), n_samples = 100)
cv  <- bruce_cv(fit, scale = "field")         # as in the original application
} # }
```
