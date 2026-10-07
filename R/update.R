#' Map a new day quickly with the hyperparameters of an earlier fit
#'
#' Operational use: learn once, map every day. The full model (hyperparameters
#' estimated) is fitted once, for example on the previous month. Each new day
#' is then mapped by refitting only a short window ending on that day, with
#' the hyperparameters fixed at those of the full fit: the spatial field,
#' network offsets, covariate effects and daily level are still estimated
#' from the new data; only range, field sd, noise levels and correlations are
#' borrowed. One run takes about a minute instead of tens of minutes.
#'
#' In the Rio de Janeiro tests (analises/rio/10), hyperparameters borrowed
#' from another month changed the leave-station-out RMSE by up to 0.07 C and
#' the coverage by up to 0.05. They can still go stale (a heat wave, a change
#' of season, a new network): refit the full model regularly, for example
#' every month. A warning is given when the day is more than `max_age` days
#' from the full fit.
#'
#' @param fit A `bruce_fit` with estimated hyperparameters (the full fit).
#' @param data A `bruce_data` with the recent observations (and background
#'   and covariates built the same way as for `fit`: same networks, same
#'   covariate names).
#' @param date Day to map; default the last day with observations.
#' @param days_before Days before `date` included in the short window (the
#'   field's day-to-day correlation lets them inform `date`). `date` itself
#'   may have no data: the map is then a one-day-ahead forecast.
#' @param grid,domain,resample,cap Passed to [bruce_predict()].
#' @param n_samples Posterior samples for the uncertainty (`0`: mean only,
#'   fastest).
#' @param max_age Warn when `date` is more than this many days from the
#'   window of `fit`.
#' @param seed Passed to [bruce_predict()].
#'
#' @return A list (class `bruce_update`) with `prediction` (a
#'   `bruce_prediction` for `date`), `fit` (the short-window fit) and
#'   `minutes`.
#' @seealso [bruce_fit()] with `theta_init` and `fix_theta`, which this wraps.
#' @export
bruce_update <- function(fit, data, date = NULL, days_before = 7, grid = NULL, domain = NULL,
                         n_samples = 0, resample = "average", cap = TRUE, max_age = 45, seed = NULL) {
  stopifnot(inherits(fit, "bruce_fit"), inherits(data, "bruce_data"))
  check_inla()
  if (isTRUE(fit$fixed_theta)) {
    warning("`fit` itself had fixed hyperparameters; borrowing from a fit that estimated them is safer",
            call. = FALSE)
  }
  cov_fit <- sort(fit$cov_names %||% character(0))
  cov_new <- sort(data$covariates$names %||% character(0))
  if (!identical(cov_fit, cov_new)) {
    stop(sprintf("Covariates differ: %s in `fit`, %s in `data`. Build `data` with the same covariates.",
                 if (length(cov_fit)) paste(cov_fit, collapse = ", ") else "none",
                 if (length(cov_new)) paste(cov_new, collapse = ", ") else "none"), call. = FALSE)
  }
  date <- as.Date(date %||% max(data$obs$date))
  window <- c(date - days_before, date)
  if (!any(data$obs$date >= window[1] & data$obs$date <= window[2])) {
    stop(sprintf("No observations between %s and %s", window[1], window[2]), call. = FALSE)
  }
  age <- min(abs(as.numeric(date - range(fit$days))))
  if (age > max_age) {
    warning(sprintf(paste0("Hyperparameters were learned on %s to %s, %d days from %s. ",
                           "Refit the full model regularly (e.g. monthly)."),
                    min(fit$days), max(fit$days), round(age), date), call. = FALSE)
  }
  spec <- fit$spec
  if (!is.null(fit$group_rho_fixed)) spec$group_rho_fixed <- fit$group_rho_fixed

  t0 <- Sys.time()
  f2 <- bruce_fit(data, mesh = fit$mesh, spec = spec, window = window, theta_init = fit,
                  fix_theta = TRUE, keep_config = n_samples > 0)
  pred <- bruce_predict(f2, grid = grid, domain = domain, dates = date, n_samples = n_samples,
                        resample = resample, cap = cap, seed = seed)
  structure(list(prediction = pred, fit = f2, date = date,
                 minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                 learned_on = range(fit$days)),
            class = "bruce_update")
}

#' @method print bruce_update
#' @export
print.bruce_update <- function(x, ...) {
  cat(sprintf("<bruce_update> map of %s (window from %s), hyperparameters learned on %s to %s, %.1f min\n",
              x$date, min(x$fit$days), x$learned_on[1], x$learned_on[2], x$minutes))
  n_obs <- sum(x$fit$data$obs$date == x$date)
  if (n_obs == 0) cat("No observations on that day: the map is a one-day-ahead forecast\n")
  print(x$prediction)
  invisible(x)
}
