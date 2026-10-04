# Cross-validation by station or by region

Refits the model with stations held out and scores predictions at the
held-out reference stations.

## Usage

``` r
bruce_cv(
  d,
  mesh,
  spec = bruce_spec(),
  window = NULL,
  scheme = c("lso", "lro"),
  regions = NULL,
  region_col = NULL,
  test_sources = NULL,
  n_samples = 200,
  level = 0.9,
  noise_sd = 0,
  checkpoint_dir = NULL,
  theta_init = NULL,
  seed = 1,
  verbose = TRUE
)
```

## Arguments

- d, mesh, spec, window:

  As in
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md).

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

- noise_sd:

  Observation noise added to samples before scoring (`0` scores the
  latent field, as in the original application; use the reference
  network's `obs_sd` to score the predictive distribution of an
  observation).

- checkpoint_dir:

  Optional directory; each finished fold is saved there and skipped on
  rerun, so long runs can be resumed.

- theta_init:

  Passed to
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md).

- seed:

  Seed for posterior sampling (`NULL` for parallel, non-reproducible
  sampling; see
  [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)).

- verbose:

  Print progress.

## Value

An object of class `bruce_cv` with per-observation `predictions`,
per-fold `folds` (scores and fit status) and `overall` scores, all as
tibbles.
