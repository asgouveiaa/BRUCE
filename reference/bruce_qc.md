# Quality control of stations against the background

Compares every station with the background product (the nearest
background cell, same day) and flags stations that behave unlike the
rest of their own network. A network's typical difference to the
background is expected (that is its bias, which the model estimates);
what is suspect is a station far from its own network, or one that
swings much more than the others from day to day. A broken humidity
sensor, a station with wrong coordinates or a sensor in direct sun show
up this way.

## Usage

``` r
bruce_qc(d, threshold = 3, min_days = 5, drop = FALSE)
```

## Arguments

- d:

  A `bruce_data`.

- threshold:

  Flag a station whose mean difference is more than `threshold` (in the
  variable's units) away from its network's median difference, or whose
  day-to-day standard deviation of the difference is above `threshold`.

- min_days:

  Flag stations with fewer days of data than this.

- drop:

  `FALSE` (default) returns the report. `TRUE` returns `d` with the
  flagged stations of biased networks removed, and the report in `d$qc`.

## Value

With `drop = FALSE`, a tibble (class `bruce_qc`) with one row per
station: `n_days`, `n_compared` (days with a background value),
`mean_diff` and `sd_diff` (station minus background), `deviation` (mean
difference minus the network's median), `flag` (reason, empty if none)
and `drop`. With `drop = TRUE`, the filtered `bruce_data`.

## Details

Without a background, each station is compared with the daily mean of
the reference stations instead.

Stations of reference networks (`bias = FALSE`) are reported but never
dropped: a reference station far from the background is more likely real
local climate (a valley, the waterfront) than a fault.

## Examples

``` r
sim <- bruce_simulate(n_days = 5, seed = 1)
d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
bruce_qc(d)
#> <bruce_qc> 52 stations compared with the background (threshold 3)
#> Station minus background, by network:
#>     source      role stations median_diff flagged dropped
#>    lowcost    biased       40        1.65       1       1
#>  reference reference       12       -0.10       0       0
#> 
#> Flagged stations:
#>  station_id  source n_days mean_diff sd_diff                  flag action
#>      IOT022 lowcost      5     -1.87    1.44 -3.5 from its network   drop
```
