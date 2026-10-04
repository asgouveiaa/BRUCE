# bruce: Bayesian urban climate estimation from heterogeneous station networks

Daily urban temperature interpolation fusing networks of different
quality with an optional gridded background, fitted with R-INLA.

## Workflow

1.  [`bruce_data()`](https://asgouveiaa.github.io/bruce/reference/bruce_data.md)
    — declare observations, stations, networks and background

2.  [`bruce_covariates()`](https://asgouveiaa.github.io/bruce/reference/bruce_covariates.md)
    — attach covariate rasters

3.  [`bruce_mesh()`](https://asgouveiaa.github.io/bruce/reference/bruce_mesh.md)
    — build the SPDE mesh

4.  [`bruce_spec()`](https://asgouveiaa.github.io/bruce/reference/bruce_spec.md)
    — priors, optional terms, collapse handling

5.  [`bruce_fit()`](https://asgouveiaa.github.io/bruce/reference/bruce_fit.md)
    — fit one window;
    [`bruce_check()`](https://asgouveiaa.github.io/bruce/reference/bruce_check.md)
    to inspect attempts

6.  [`bruce_predict()`](https://asgouveiaa.github.io/bruce/reference/bruce_predict.md)
    — rasters of mean, sd and quantiles

7.  [`bruce_cv()`](https://asgouveiaa.github.io/bruce/reference/bruce_cv.md)
    — leave-station-out / leave-region-out validation

8.  [`bruce_zonal()`](https://asgouveiaa.github.io/bruce/reference/bruce_zonal.md)
    — weighted means by zone (e.g. population exposure)

9.  [`bruce_decompose()`](https://asgouveiaa.github.io/bruce/reference/bruce_decompose.md)
    — covariate vs spatial-field share of the surface, with
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

[`bruce_simulate()`](https://asgouveiaa.github.io/bruce/reference/bruce_simulate.md)
generates a synthetic city for testing.

## Author

**Maintainer**: Ayrton Gouveia <ayrtonsgouveia19@gmail.com>
