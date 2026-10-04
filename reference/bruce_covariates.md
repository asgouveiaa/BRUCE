# Attach spatial covariates

Extracts covariates at station locations and at background cells,
standardises them, and keeps the rasters for prediction.

## Usage

``` r
bruce_covariates(
  d,
  rasters,
  background_method = "buffer_mean",
  background_buffer = 5500,
  scale_on = NULL
)
```

## Arguments

- d:

  A `bruce_data` object.

- rasters:

  A
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  whose layer names are the covariate names (syntactic R names), e.g.
  `elev`, `urb`, `svf`.

- background_method:

  How each covariate is represented at background cells: `"buffer_mean"`
  (mean within `background_buffer` metres, matching the footprint of a
  coarse grid cell), `"point"` (value at the cell centre), or `"zero"`
  (set to the covariate mean; for covariates that have no meaning at the
  background's scale, such as sky view factor). A single value or a
  named vector per covariate.

- background_buffer:

  Buffer radius in metres for `"buffer_mean"`.

- scale_on:

  Which locations define the standardisation of each covariate:
  `"background"` or `"stations"`. Defaults to `"background"` when a
  background exists and the covariate is not `"zero"`, otherwise
  `"stations"`.

## Value

`d` with a `covariates` element.

## Details

Covariate values that are missing after extraction are set to the
covariate mean (0 on the standardised scale). The observed range at
stations and background cells is stored and used by
[`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)
to cap covariates on the prediction grid, so the model never
extrapolates a linear effect beyond the data.
