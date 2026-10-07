#' Decompose the predicted surface into its components
#'
#' Splits the posterior-mean prediction on a grid into the part explained by
#' covariates, the city-wide daily term and the spatial field, and measures
#' how much of the spatial variance of the surface comes from each.
#'
#' @inheritParams bruce_predict
#' @param fit A `bruce_fit`.
#'
#' @details For each day the prediction at cell \eqn{s} is
#'   \deqn{\hat y(s) = \beta_0 + \underbrace{x(s)^\top\beta}_{\text{covariates}} +
#'   \underbrace{d_t}_{\text{day}} + \underbrace{w(s,t)}_{\text{field}}}
#'   using posterior means. The intercept and the daily term are constant in
#'   space, so the spatial variance across cells is
#'   \deqn{\mathrm{Var}(\hat y) = \mathrm{Var}(x^\top\beta) + \mathrm{Var}(w) +
#'   2\,\mathrm{Cov}(x^\top\beta, w).}
#'   `share_covariates` and `share_field` are each component's variance over
#'   `var_covariates + var_field`; the covariance term is reported separately
#'   so that a large overlap between the two parts is visible rather than
#'   hidden in the shares.
#'
#' @return An object of class `bruce_decomposition`: `total`, `covariates`,
#'   `field` (SpatRasters, one layer per day), `by_covariate` (list of
#'   SpatRasters, one per covariate), `day_effect` (one value per day),
#'   `variance` (one row per day) and `stations` (locations used in the fit,
#'   for plotting).
#' @export
bruce_decompose <- function(fit, grid = NULL, domain = NULL, dates = NULL, cap = TRUE,
                            resample = "average") {
  stopifnot(inherits(fit, "bruce_fit"))
  check_inla()
  days <- fit$days
  dates <- if (is.null(dates)) days else as.Date(dates)
  if (!all(dates %in% days)) stop("Some `dates` are outside the fitted window", call. = FALSE)
  g_idx <- match(dates, days)

  gp <- prepare_grid(fit, grid, domain, cap, resample)
  comps <- latent_components(fit, 0, NULL)
  beta <- comps$fixed[, 1]
  names(beta) <- c("intercept", fit$cov_names)
  n <- fit$spde$n.spde
  nc <- terra::ncell(gp$template)

  per_cov <- lapply(stats::setNames(fit$cov_names, fit$cov_names),
                    function(cn) gp$X[, cn] * beta[[cn]])
  cov_part <- if (length(per_cov)) Reduce(`+`, per_cov) else rep(0, length(gp$cells))

  m_tot <- m_field <- matrix(NA_real_, nc, length(dates))
  var_rows <- vector("list", length(dates))
  for (k in seq_along(dates)) {
    g <- g_idx[k]
    field <- as.vector(gp$A %*% comps$field[((g - 1) * n + 1):(g * n), 1])
    day <- comps$day[g, 1]
    m_field[gp$cells, k] <- field
    m_tot[gp$cells, k] <- beta[["intercept"]] + cov_part + day + field
    vc <- stats::var(cov_part); vf <- stats::var(field); cv <- stats::cov(cov_part, field)
    var_rows[[k]] <- data.frame(
      date = dates[k], var_total = stats::var(m_tot[gp$cells, k]),
      var_covariates = vc, var_field = vf, cov_term = 2 * cv,
      share_covariates = vc / (vc + vf), share_field = vf / (vc + vf),
      field_min = min(field), field_max = max(field))
  }

  to_rast <- function(m) {
    r <- terra::rast(gp$template, nlyrs = ncol(m))
    terra::values(r) <- m
    names(r) <- as.character(dates)
    r
  }
  one_layer <- function(v) {
    m <- matrix(NA_real_, nc, 1); m[gp$cells, 1] <- v
    r <- terra::rast(gp$template); terra::values(r) <- m; r
  }

  st <- unique(fit$data$obs[fit$data$obs$date %in% days &
                              !(fit$data$obs$station_id %in% fit$drop_stations),
                            c("station_id", "source", "x", "y")])

  structure(list(
    total = to_rast(m_tot),
    covariates = one_layer(cov_part),
    field = to_rast(m_field),
    by_covariate = lapply(per_cov, one_layer),
    day_effect = stats::setNames(comps$day[g_idx, 1], as.character(dates)),
    intercept = beta[["intercept"]],
    variance = tibble::as_tibble(do.call(rbind, var_rows)),
    stations = st, dates = dates),
    class = "bruce_decomposition")
}

#' @method print bruce_decomposition
#' @export
print.bruce_decomposition <- function(x, ...) {
  cat(sprintf("<bruce_decomposition> %d day(s)\n", length(x$dates)))
  v <- x$variance
  cat(sprintf("Share of spatial variance (mean over days): covariates %.0f%%, field %.0f%%\n",
              100 * mean(v$share_covariates), 100 * mean(v$share_field)))
  print(data.frame(date = v$date, round(v[, -1], 3)), row.names = FALSE)
  invisible(x)
}

#' Plot a decomposition
#'
#' Three maps for one day: total prediction, covariate part (centred, so it
#' reads as a deviation from the city mean) and spatial field with the
#' stations used in the fit. The covariate and field panels share a
#' symmetric colour scale so their magnitudes can be compared directly.
#'
#' @param x A `bruce_decomposition`.
#' @param date Day to plot (default: the first).
#' @param stations Overlay station locations on the field panel.
#' @param ... Unused.
#' @method plot bruce_decomposition
#' @export
plot.bruce_decomposition <- function(x, date = NULL, stations = TRUE, ...) {
  k <- if (is.null(date)) 1 else match(as.Date(date), x$dates)
  if (is.na(k)) stop("`date` is not in the decomposition", call. = FALSE)
  tot <- x$total[[k]]
  cov <- x$covariates - mean(terra::values(x$covariates), na.rm = TRUE)
  fld <- x$field[[k]]
  lim <- 1.001 * max(abs(c(terra::values(cov), terra::values(fld))), na.rm = TRUE)
  div <- grDevices::hcl.colors(50, "Blue-Red 3")

  op <- graphics::par(mfrow = c(1, 3), mar = c(2, 2, 3, 4))
  on.exit(graphics::par(op))
  terra::plot(tot, main = sprintf("Total (%s)", format(x$dates[k])),
              col = grDevices::hcl.colors(50, "YlOrRd", rev = TRUE))
  terra::plot(cov, main = sprintf("Covariates (%.0f%%)", 100 * x$variance$share_covariates[k]),
              col = div, range = c(-lim, lim))
  terra::plot(fld, main = sprintf("Spatial field (%.0f%%)", 100 * x$variance$share_field[k]),
              col = div, range = c(-lim, lim))
  if (stations && nrow(x$stations)) {
    src <- factor(x$stations$source)
    pch <- c(21, 24, 22, 23, 25)[(as.integer(src) - 1) %% 5 + 1]
    graphics::points(x$stations$x, x$stations$y, pch = pch, cex = 0.8, bg = "white")
    terra::add_legend("topleft", legend = levels(src), pch = c(21, 24, 22, 23, 25)[seq_along(levels(src))],
                     pt.bg = "white", bg = "white", box.col = "grey60", cex = 0.85)
  }
  invisible(x)
}
