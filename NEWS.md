# bruce 0.1.2

* `bruce_mesh()` builds the domain boundary with INLA's non-convex hull
  algorithm (as `INLA::inla.nonconvex.hull()`), so meshes match those of
  `INLA::inla.mesh.2d()` from the same points. Before, older fmesher versions
  produced a different (sf buffer) boundary, which changed the mesh and could
  shift the estimated spatial range.

# bruce 0.1.1

* `bruce_covariates()` accepts a named list of rasters, each kept in its own
  projection and resolution and read directly at stations and background
  cells, with no resampling to a common grid. New `na_value` argument sets
  the raw value used where a raster has no data (e.g. `c(urb = 0)`).
* `bruce_predict()` and `bruce_decompose()` read each covariate from its
  native raster: `resample = "near"` takes the raster's value at each grid
  cell centre; `"average"` and `"bilinear"` project onto the grid.
* Validated against the Rio de Janeiro article fit (January 2025, daily
  median): with identical inputs, all fixed effects within 0.03 posterior sd
  and spatial range 7682 m vs 7669 m.

# bruce 0.1.0

* First version: data contract (`bruce_data()`, `bruce_covariates()`),
  mesh, model specification and fitting with collapse detection
  (`bruce_mesh()`, `bruce_spec()`, `bruce_fit()`, `bruce_check()`),
  prediction, decomposition and zonal aggregation (`bruce_predict()`,
  `bruce_decompose()`, `bruce_zonal()`), leave-station-out and
  leave-region-out validation with CRPS (`bruce_cv()`) and a synthetic
  city simulator (`bruce_simulate()`).
