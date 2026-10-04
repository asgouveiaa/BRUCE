# Assemble and validate input data

The data contract of the package. Everything specific to a city enters
here: which networks exist, which ones are trusted, and where the
stations are.

## Usage

``` r
bruce_data(obs, stations, sources, crs, background = NULL, background_sd = 0.7)
```

## Arguments

- obs:

  Data frame of station observations with columns `station_id`, `date`
  (coercible to `Date`) and `value` (e.g. daily median or maximum
  temperature).

- stations:

  Data frame with one row per station: `station_id`, `source`, and
  either `x`/`y` in `crs` or `lon`/`lat` in WGS84.

- sources:

  Data frame describing each network, with columns

  - `source`: network name, matching `stations$source`;

  - `bias`: `FALSE` for reference networks (e.g. official stations),
    `TRUE` for networks that get an estimated offset (low-cost sensors,
    citizen stations, airports with a known systematic difference);

  - `obs_sd`: fixed observation standard deviation, or `NA` to estimate
    it with a PC prior. At least one reference network (`bias = FALSE`)
    is needed to identify the intercept.

- crs:

  Projected coordinate reference system in metres (EPSG code or anything
  accepted by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)).

- background:

  Optional gridded product on the same variable (e.g. ERA5-Land daily
  statistic), as a data frame with `date`, `value` and `x`/`y` or
  `lon`/`lat` of cell centres. Enters as an additional unbiased
  likelihood with fixed standard deviation `background_sd`.

- background_sd:

  Fixed observation standard deviation of the background.

## Value

An object of class `bruce_data`.
