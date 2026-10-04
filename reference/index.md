# Package index

## Data

Declare what the city has: observations, stations, networks, background
and covariates.

- [`bruce_data()`](https://asgouveiaa.github.io/bruce/reference/bruce_data.md)
  : Assemble and validate input data
- [`bruce_covariates()`](https://asgouveiaa.github.io/bruce/reference/bruce_covariates.md)
  : Attach spatial covariates

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

- [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)
  : Predict on a raster grid
- [`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md)
  : Decompose the predicted surface into its components
- [`plot(`*`<bruce_decomposition>`*`)`](https://asgouveiaa.github.io/bruce/reference/plot.bruce_decomposition.md)
  : Plot a decomposition
- [`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md)
  : Aggregate predictions over zones
- [`as_tibble(`*`<bruce_prediction>`*`)`](https://asgouveiaa.github.io/bruce/reference/as_tibble.bruce.md)
  [`as_tibble(`*`<bruce_decomposition>`*`)`](https://asgouveiaa.github.io/bruce/reference/as_tibble.bruce.md)
  : Prediction or decomposition as a long table

## Validation

- [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)
  : Cross-validation by station or by region
- [`bruce_simulate()`](https://asgouveiaa.github.io/bruce/reference/bruce_simulate.md)
  : Simulate a synthetic city
