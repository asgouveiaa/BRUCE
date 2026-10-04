#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble

#' Prediction or decomposition as a long table
#'
#' One row per grid cell and day, ready for `dplyr` and `ggplot2`
#' (e.g. `geom_raster(aes(x, y, fill = mean))`).
#'
#' @param x A `bruce_prediction` (from [bruce_predict()]) or a
#'   `bruce_decomposition` (from [bruce_decompose()]).
#' @param ... Unused.
#'
#' @return A tibble with the cell centre coordinates `x`, `y` (in the grid's
#'   projection) and `date`, plus
#'   * for a prediction: `mean` and, when posterior samples were drawn, `sd`
#'     and one column per quantile (`q0.05`, `q0.95`, ...);
#'   * for a decomposition: `total`, `covariates`, `field` and one column per
#'     covariate's contribution (`cov_elev`, ...).
#'
#'   Cells outside the domain are dropped.
#'
#' @examples
#' \dontrun{
#' pred <- bruce_predict(fit, n_samples = 200)
#' as_tibble(pred) |>
#'   dplyr::filter(date == as.Date("2025-01-15")) |>
#'   ggplot2::ggplot(ggplot2::aes(x, y, fill = mean)) +
#'   ggplot2::geom_raster() + ggplot2::coord_equal()
#' }
#' @name as_tibble.bruce
NULL

#' @rdname as_tibble.bruce
#' @export
as_tibble.bruce_prediction <- function(x, ...) {
  long <- function(r) as.vector(terra::values(r, mat = TRUE)[x$cells, , drop = FALSE])
  xy <- terra::xyFromCell(x$template, x$cells)
  n <- length(x$cells); nd <- length(x$dates)
  dates <- rep(x$dates, each = n); m <- long(x$mean)
  out <- tibble::tibble(x = rep(xy[, 1], nd), y = rep(xy[, 2], nd), date = dates, mean = m)
  if (!is.null(x$sd)) out$sd <- long(x$sd)
  for (q in names(x$quantiles)) out[[q]] <- long(x$quantiles[[q]])
  out
}

#' @rdname as_tibble.bruce
#' @export
as_tibble.bruce_decomposition <- function(x, ...) {
  cells <- which(!is.na(terra::values(x$total[[1]], mat = FALSE)))
  xy <- terra::xyFromCell(x$total, cells)
  n <- length(cells); nd <- length(x$dates)
  long <- function(r) as.vector(terra::values(r, mat = TRUE)[cells, , drop = FALSE])
  dates <- rep(x$dates, each = n)
  tot <- long(x$total); cov <- rep(long(x$covariates), nd); fld <- long(x$field)
  out <- tibble::tibble(x = rep(xy[, 1], nd), y = rep(xy[, 2], nd), date = dates,
                        total = tot, covariates = cov, field = fld)
  for (cn in names(x$by_covariate)) out[[paste0("cov_", cn)]] <- rep(long(x$by_covariate[[cn]]), nd)
  out
}
