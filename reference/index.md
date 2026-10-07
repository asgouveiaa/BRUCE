# Package index

## Data

Declare what the city has: observations, stations, networks, background
and covariates.

- [`bruce_data()`](https://asgouveiaa.github.io/bruce/reference/bruce_data.md)
  : Assemble and validate input data
- [`bruce_covariates()`](https://asgouveiaa.github.io/bruce/reference/bruce_covariates.md)
  : Attach spatial covariates
- [`bruce_qc()`](https://asgouveiaa.github.io/bruce/reference/bruce_qc.md)
  : Quality control of stations against the background

## Model

- [`bruce_mesh()`](https://asgouveiaa.github.io/bruce/reference/bruce_mesh.md)
  : Build the SPDE mesh
- [`bruce_spec()`](https://asgouveiaa.github.io/bruce/reference/bruce_spec.md)
  : Model specification
- [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md)
  : Fit the model for one time window
- [`bruce_check()`](https://asgouveiaa.github.io/bruce/reference/bruce_check.md)
  : Diagnostics of a fit

## Use the fit

- [`bruce_grid()`](https://asgouveiaa.github.io/bruce/reference/bruce_grid.md)
  : Prediction grid at a chosen resolution
- [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)
  : Predict on a raster grid
- [`bruce_update()`](https://asgouveiaa.github.io/bruce/reference/bruce_update.md)
  : Map a new day quickly with the hyperparameters of an earlier fit
- [`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md)
  : Decompose the predicted surface into its components
- [`plot(`*`<bruce_decomposition>`*`)`](https://asgouveiaa.github.io/bruce/reference/plot.bruce_decomposition.md)
  : Plot a decomposition
- [`plot(`*`<bruce_prediction>`*`)`](https://asgouveiaa.github.io/bruce/reference/plot.bruce_prediction.md)
  : Plot a prediction: one map per day
- [`bruce_plot_series()`](https://asgouveiaa.github.io/bruce/reference/bruce_plot_series.md)
  : City series: the map against the reference stations
- [`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md)
  : Aggregate predictions over zones
- [`as_tibble(`*`<bruce_prediction>`*`)`](https://asgouveiaa.github.io/bruce/reference/as_tibble.bruce.md)
  [`as_tibble(`*`<bruce_decomposition>`*`)`](https://asgouveiaa.github.io/bruce/reference/as_tibble.bruce.md)
  : Prediction or decomposition as a long table

## Validation

- [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)
  : Cross-validation by station or by region
- [`plot(`*`<bruce_cv>`*`)`](https://asgouveiaa.github.io/bruce/reference/plot.bruce_cv.md)
  : Plot a cross-validation
- [`bruce_simulate()`](https://asgouveiaa.github.io/bruce/reference/bruce_simulate.md)
  : Simulate a synthetic city
