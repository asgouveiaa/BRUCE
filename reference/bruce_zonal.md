# Aggregate predictions over zones

Weighted zonal means (e.g. population-weighted exposure by neighbourhood
or health territory), propagating posterior uncertainty when samples
were kept.

## Usage

``` r
bruce_zonal(pred, zones, id_col, weights = NULL, probs = c(0.05, 0.95))
```

## Arguments

- pred:

  A `bruce_prediction`.

- zones:

  Polygons (`sf`).

- id_col:

  Name of the zone id column.

- weights:

  Optional
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  of weights (e.g. population). Resampled to the prediction grid by
  averaging, so only relative weights within a zone matter.

- probs:

  Quantiles of the zonal mean (only with samples).

## Value

A data frame with one row per zone and date: `mean` and, with samples,
`sd` and the requested quantiles of the zonal mean.
