#' bruce: Bayesian urban climate estimation from heterogeneous station networks
#'
#' Daily urban temperature interpolation fusing networks of different quality
#' with an optional gridded background, fitted with R-INLA.
#'
#' @section Workflow:
#' 1. [bruce_data()] — declare observations, stations, networks and background
#' 2. [bruce_covariates()] — attach covariate rasters
#' 3. [bruce_mesh()] — build the SPDE mesh
#' 4. [bruce_spec()] — priors, optional terms, collapse handling
#' 5. [bruce_fit()] — fit one window; [bruce_check()] to inspect attempts
#' 6. [bruce_predict()] — rasters of mean, sd and quantiles
#' 7. [bruce_cv()] — leave-station-out / leave-region-out validation
#' 8. [bruce_zonal()] — weighted means by zone (e.g. population exposure)
#' 9. [bruce_decompose()] — covariate vs spatial-field share of the surface, with `plot()`
#'
#' [bruce_simulate()] generates a synthetic city for testing.
#'
#' @keywords internal
"_PACKAGE"
