# Plot a prediction: one map per day

One small map per day. With `scale = "common"` all days share one colour
scale, so days can be compared directly (from the 0.5% to the 99.5%
quantile of all days, so a few extreme cells do not wash out the rest).
With `scale = "free"` each day gets its own scale and legend, which
shows the pattern within the city on every day, even on days much cooler
or warmer than the rest. With `scale = "anomaly"` each day is shown
minus its own city mean, on one symmetric scale: how much hotter or
cooler each place is than the city that day, comparable across days.

## Usage

``` r
# S3 method for class 'bruce_prediction'
plot(
  x,
  what = c("mean", "sd"),
  scale = c("common", "free", "anomaly"),
  dates = NULL,
  ncol = NULL,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A `bruce_prediction`.

- what:

  `"mean"` (default) or `"sd"` (uncertainty; needs a prediction made
  with `n_samples > 0`).

- scale:

  `"common"` (default), `"free"` or `"anomaly"` (see Details;
  `"anomaly"` only for `what = "mean"`).

- dates:

  Days to plot (default: all).

- ncol:

  Maps per row (default: up to 6).

- main:

  Overall title.

- ...:

  Unused.
