#' Model specification
#'
#' Collects priors, optional terms and the fitting strategy. Defaults
#' reproduce the model of the original Rio de Janeiro application, except
#' where noted.
#'
#' @param range_prior PC prior on the Matern range: `"nn"` anchors the median
#'   at the median nearest-neighbour distance between stations (the
#'   resolution limit of the network), or give a number in metres. The Rio de
#'   Janeiro application uses `15000`; `"nn"` is the default here because it
#'   adapts to any network without a city-specific constant.
#' @param range_prior_prob Probability that the range is below `range_prior`.
#' @param sigma_prior PC prior on the field standard deviation, `c(u, p)`
#'   meaning P(sigma > u) = p.
#' @param group_rho_prior PC prior (`pc.cor1`) on the day-to-day AR(1)
#'   correlation of the spatial field.
#' @param group_rho_fixed Fix that correlation at a value instead of
#'   estimating it.
#' @param station_effect Include an iid station random effect (absorbs
#'   persistent station-level siting bias).
#' @param station_prior PC prior on its standard deviation.
#' @param day_effect Include a city-wide daily AR(1) term.
#' @param day_rho_prior,day_prec_prior Priors of the daily term.
#' @param day_sd_fixed Fix the daily term's standard deviation instead of
#'   estimating it.
#' @param biased_prec_prior PC prior on the observation standard deviation of
#'   networks with `obs_sd = NA`.
#' @param intercept_prior `c(mean, precision)` for the intercept; `NA` mean
#'   uses the mean of reference observations.
#' @param covariate_prec,bias_prec Prior precision of covariate and network
#'   offset coefficients (mean 0).
#' @param int_strategy INLA integration strategy.
#' @param on_collapse What to do when the fit degenerates (field day-to-day
#'   correlation near 1, or field standard deviation exploding):
#'   `"retry_then_fix"` restarts from a reset mode and, if it collapses again,
#'   refits with the correlation fixed at `fallback_group_rho`; `"fix"` goes
#'   straight to the fixed refit; `"warn"` only warns; `"none"` does nothing.
#' @param collapse_rho,collapse_sd_ratio Thresholds of the collapse check: a
#'   fit is flagged when the posterior mean correlation exceeds
#'   `collapse_rho` or the field standard deviation exceeds
#'   `collapse_sd_ratio` times the standard deviation of the observations.
#' @param fallback_group_rho Value used when the correlation must be fixed.
#' @param num_threads Passed to INLA (`NULL` uses INLA's default).
#'
#' @return An object of class `bruce_spec`.
#' @export
bruce_spec <- function(range_prior = "nn", range_prior_prob = 0.5,
                       sigma_prior = c(8, 0.01),
                       group_rho_prior = c(0.5, 0.75), group_rho_fixed = NULL,
                       station_effect = TRUE, station_prior = c(1, 0.01),
                       day_effect = TRUE, day_rho_prior = c(0.5, 0.75),
                       day_prec_prior = c(3, 0.01), day_sd_fixed = NULL,
                       biased_prec_prior = c(3, 0.01),
                       intercept_prior = c(NA, 0.001),
                       covariate_prec = 1 / 4, bias_prec = 1 / 9,
                       int_strategy = "ccd",
                       on_collapse = c("retry_then_fix", "fix", "warn", "none"),
                       collapse_rho = 0.98, collapse_sd_ratio = 10,
                       fallback_group_rho = 0.55,
                       num_threads = NULL) {
  on_collapse <- match.arg(on_collapse)
  if (!(identical(range_prior, "nn") || (is.numeric(range_prior) && range_prior > 0))) {
    stop("`range_prior` must be 'nn' or a positive number of metres", call. = FALSE)
  }
  if (!is.null(group_rho_fixed) && (group_rho_fixed <= -1 || group_rho_fixed >= 1)) {
    stop("`group_rho_fixed` must be in (-1, 1)", call. = FALSE)
  }
  structure(as.list(environment()), class = "bruce_spec")
}
