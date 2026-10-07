# City series: the map against the reference stations

For each day, the median of the mapped surface over the grid and its
5-95% spread across cells, with the daily mean of the reference
stations' readings. A quick check that the map follows the stations
through the period (level and day-to-day changes).

## Usage

``` r
bruce_plot_series(prediction, data, sources = NULL, ylab = "Value")
```

## Arguments

- prediction:

  A `bruce_prediction`.

- data:

  The `bruce_data` (or `bruce_fit`) with the observations.

- sources:

  Networks averaged as points (default: the reference networks).

- ylab:

  Axis label.

## Value

A data frame with the daily values, invisibly.
