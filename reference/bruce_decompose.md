# Decompose the predicted surface into its components

Splits the posterior-mean prediction on a grid into the part explained
by covariates, the city-wide daily term and the spatial field, and
measures how much of the spatial variance of the surface comes from
each.

## Usage

``` r
bruce_decompose(
  fit,
  grid = NULL,
  domain = NULL,
  dates = NULL,
  cap = TRUE,
  resample = "average"
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

## Value

An object of class `bruce_decomposition`: `total`, `covariates`, `field`
(SpatRasters, one layer per day), `by_covariate` (list of SpatRasters,
one per covariate), `day_effect` (one value per day), `variance` (one
row per day) and `stations` (locations used in the fit, for plotting).

## Details

For each day the prediction at cell \\s\\ is \$\$\hat y(s) = \beta_0 +
\underbrace{x(s)^\top\beta}\_{\text{covariates}} +
\underbrace{d_t}\_{\text{day}} + \underbrace{w(s,t)}\_{\text{field}}\$\$
using posterior means. The intercept and the daily term are constant in
space, so the spatial variance across cells is \$\$\mathrm{Var}(\hat y)
= \mathrm{Var}(x^\top\beta) + \mathrm{Var}(w) +
2\\\mathrm{Cov}(x^\top\beta, w).\$\$ `share_covariates` and
`share_field` are each component's variance over
`var_covariates + var_field`; the covariance term is reported separately
so that a large overlap between the two parts is visible rather than
hidden in the shares.
