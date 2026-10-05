# Prediction grid at a chosen resolution

Builds a raster covering a polygon at the requested cell size, with
cells outside the polygon set to `NA` (cell centre inside the polygon),
ready for the `grid` argument of
[`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md),
[`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md)
and
[`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md).

## Usage

``` r
bruce_grid(domain, resolution, crs = NULL)
```

## Arguments

- domain:

  Polygon(s) (`sf`, `sfc` or `SpatVector`).

- resolution:

  Cell size in the units of the projection (metres).

- crs:

  Projection of the grid: anything
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  accepts, or a `bruce_data` / `bruce_fit` to use the model's
  projection. `NULL` keeps the domain's projection, which must then be
  projected.

## Value

A single-layer
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
named `grid`, 1 inside the domain and `NA` outside.

## Details

The grid resolution sets where the surface is evaluated, not how much
spatial detail the model has. The spatial field varies on the scale of
the mesh and of the station network (kilometres), so a finer grid adds
detail only through the covariates. Two consequences:

- covariates need native data at least as fine as the grid, or the map
  shows the blocks of the coarser raster;

- a coefficient refers to the covariate as it was measured at the
  stations. If urban fraction was the mean over a 500 m block,
  predicting at 100 m with urban fraction averaged over 100 m cells
  gives larger contrasts than the model learned. Keep the covariate's
  support (e.g. a 500 m moving-window mean computed at 100 m) when
  refining the grid.

Prediction time and memory grow with the number of cells: halving the
cell size multiplies them by four, and posterior samples kept with
`keep_samples = TRUE` take cells x samples x days x 8 bytes.

## Examples

``` r
sim <- bruce_simulate(n_days = 2, seed = 1)
g <- bruce_grid(sim$domain, resolution = 1000)
g
#> class       : SpatRaster
#> size        : 30, 30, 1  (nrow, ncol, nlyr)
#> resolution  : 1000, 1000  (x, y)
#> extent      : 650000, 680000, 7450000, 7480000  (xmin, xmax, ymin, ymax)
#> coord. ref. : SIRGAS 2000 / UTM zone 23S (EPSG:31983)
#> source(s)   : memory
#> name        : grid
#> min value   :    1
#> max value   :    1
```
