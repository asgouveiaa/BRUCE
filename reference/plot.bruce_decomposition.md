# Plot a decomposition

Three maps for one day: total prediction, covariate part (centred, so it
reads as a deviation from the city mean) and spatial field with the
stations used in the fit. The covariate and field panels share a
symmetric colour scale so their magnitudes can be compared directly.

## Usage

``` r
# S3 method for class 'bruce_decomposition'
plot(x, date = NULL, stations = TRUE, ...)
```

## Arguments

- x:

  A `bruce_decomposition`.

- date:

  Day to plot (default: the first).

- stations:

  Overlay station locations on the field panel.

- ...:

  Unused.
