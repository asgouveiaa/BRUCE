# Changelog

## bruce 0.1.1

- [`bruce_covariates()`](https://asgouveiaa.github.io/bruce/reference/bruce_covariates.md)
  accepts a named list of rasters, each kept in its own projection and
  resolution and read directly at stations and background cells, with no
  resampling to a common grid. New `na_value` argument sets the raw
  value used where a raster has no data (e.g. `c(urb = 0)`).
- [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)
  and
  [`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md)
  read each covariate from its native raster: `resample = "near"` takes
  the raster’s value at each grid cell centre; `"average"` and
  `"bilinear"` project onto the grid.
- Validated against the Rio de Janeiro article fit (January 2025, daily
  median): with identical inputs, all fixed effects within 0.03
  posterior sd and spatial range 7682 m vs 7669 m.

## bruce 0.1.0

- First version: data contract
  ([`bruce_data()`](https://asgouveiaa.github.io/bruce/reference/bruce_data.md),
  [`bruce_covariates()`](https://asgouveiaa.github.io/bruce/reference/bruce_covariates.md)),
  mesh, model specification and fitting with collapse detection
  ([`bruce_mesh()`](https://asgouveiaa.github.io/bruce/reference/bruce_mesh.md),
  [`bruce_spec()`](https://asgouveiaa.github.io/bruce/reference/bruce_spec.md),
  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md),
  [`bruce_check()`](https://asgouveiaa.github.io/bruce/reference/bruce_check.md)),
  prediction, decomposition and zonal aggregation
  ([`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md),
  [`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md),
  [`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md)),
  leave-station-out and leave-region-out validation with CRPS
  ([`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md))
  and a synthetic city simulator
  ([`bruce_simulate()`](https://asgouveiaa.github.io/bruce/reference/bruce_simulate.md)).
