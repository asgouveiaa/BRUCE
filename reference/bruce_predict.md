# Predict on a raster grid

Predicts the variable at the reference-network scale (no network offset,
no station effect) for each day of the fitted window.

## Usage

``` r
bruce_predict(
  fit,
  grid = NULL,
  domain = NULL,
  dates = NULL,
  n_samples = 0,
  probs = c(0.05, 0.95),
  keep_samples = FALSE,
  cap = TRUE,
  resample = "average",
  seed = 1
)
```

## Arguments

- fit:

  A `bruce_fit`.

- grid:

  A
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  defining the prediction grid. Cells where its first layer is `NA` are
  skipped, so a masked raster works as a domain. `NULL` uses the
  covariate rasters' own grid.

- domain:

  Optional polygon (`sf` or `SpatVector`) to mask the grid.

- dates:

  Subset of the fitted days; `NULL` predicts all.

- n_samples:

  Number of joint posterior samples. `0` returns posterior means only
  (fast, from INLA's marginal summaries). `> 0` returns mean, standard
  deviation and quantiles from samples; requires the fit to have been
  made with `keep_config = TRUE`.

- probs:

  Quantiles to return when sampling.

- keep_samples:

  Keep the per-cell samples (needed by
  [`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md)
  for uncertainty of zonal means). Memory grows as cells x samples x
  days.

- cap:

  Cap covariates to the range observed at stations and background cells
  (avoids extrapolating linear effects).

- resample:

  How each covariate is brought to the grid: `"average"` (block mean of
  the native raster over each grid cell, via
  [`terra::project()`](https://rspatial.github.io/terra/reference/project.html)),
  `"bilinear"`, or `"near"` (the native raster's value at the grid cell
  centre, read directly from that raster in its own projection). A
  single value or a named vector, e.g.
  `c(elev = "average", svf = "near")`. Rasters that already share the
  grid's geometry are used as they are.

- seed:

  Seed for posterior sampling. A number makes samples reproducible (INLA
  then samples on one thread); `NULL` samples in parallel, faster but
  not reproducible.

## Value

An object of class `bruce_prediction` with `mean`, `sd` and `quantiles`
(SpatRasters with one layer per day) and, optionally, `samples`.
