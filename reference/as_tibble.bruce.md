# Prediction or decomposition as a long table

One row per grid cell and day, ready for `dplyr` and `ggplot2` (e.g.
`geom_raster(aes(x, y, fill = mean))`).

## Usage

``` r
# S3 method for class 'bruce_prediction'
as_tibble(x, ...)

# S3 method for class 'bruce_decomposition'
as_tibble(x, ...)
```

## Arguments

- x:

  A `bruce_prediction` (from
  [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md))
  or a `bruce_decomposition` (from
  [`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md)).

- ...:

  Unused.

## Value

A tibble with the cell centre coordinates `x`, `y` (in the grid's
projection) and `date`, plus

- for a prediction: `mean` and, when posterior samples were drawn, `sd`
  and one column per quantile (`q0.05`, `q0.95`, ...);

- for a decomposition: `total`, `covariates`, `field` and one column per
  covariate's contribution (`cov_elev`, ...).

Cells outside the domain are dropped.

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- bruce_predict(fit, n_samples = 200)
as_tibble(pred) |>
  dplyr::filter(date == as.Date("2025-01-15")) |>
  ggplot2::ggplot(ggplot2::aes(x, y, fill = mean)) +
  ggplot2::geom_raster() + ggplot2::coord_equal()
} # }
```
