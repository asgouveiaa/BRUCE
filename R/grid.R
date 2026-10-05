#' Prediction grid at a chosen resolution
#'
#' Builds a raster covering a polygon at the requested cell size, with cells
#' outside the polygon set to `NA` (cell centre inside the polygon), ready
#' for the `grid` argument of [bruce_predict()], [bruce_decompose()] and
#' [bruce_zonal()].
#'
#' @param domain Polygon(s) (`sf`, `sfc` or `SpatVector`).
#' @param resolution Cell size in the units of the projection (metres).
#' @param crs Projection of the grid: anything [sf::st_crs()] accepts, or a
#'   `bruce_data` / `bruce_fit` to use the model's projection. `NULL` keeps
#'   the domain's projection, which must then be projected.
#'
#' @details The grid resolution sets where the surface is evaluated, not how
#'   much spatial detail the model has. The spatial field varies on the scale
#'   of the mesh and of the station network (kilometres), so a finer grid
#'   adds detail only through the covariates. Two consequences:
#'   * covariates need native data at least as fine as the grid, or the map
#'     shows the blocks of the coarser raster;
#'   * a coefficient refers to the covariate as it was measured at the
#'     stations. If urban fraction was the mean over a 500 m block, predicting
#'     at 100 m with urban fraction averaged over 100 m cells gives larger
#'     contrasts than the model learned. Keep the covariate's support (e.g. a
#'     500 m moving-window mean computed at 100 m) when refining the grid.
#'
#'   Prediction time and memory grow with the number of cells: halving the
#'   cell size multiplies them by four, and posterior samples kept with
#'   `keep_samples = TRUE` take cells x samples x days x 8 bytes.
#'
#' @return A single-layer [terra::SpatRaster] named `grid`, 1 inside the
#'   domain and `NA` outside.
#' @examples
#' sim <- bruce_simulate(n_days = 2, seed = 1)
#' g <- bruce_grid(sim$domain, resolution = 1000)
#' g
#' @export
bruce_grid <- function(domain, resolution, crs = NULL) {
  stopifnot(is.numeric(resolution), length(resolution) %in% 1:2, all(resolution > 0))
  dom <- if (inherits(domain, "SpatVector")) sf::st_as_sf(domain) else sf::st_as_sf(sf::st_geometry(domain))
  if (inherits(crs, "bruce_fit")) crs <- crs$data$crs
  if (inherits(crs, "bruce_data")) crs <- crs$crs
  if (!is.null(crs)) dom <- sf::st_transform(dom, crs)
  if (isTRUE(sf::st_is_longlat(dom))) {
    stop("`domain` is in longitude/latitude; give a projected `crs` (e.g. the model's)", call. = FALSE)
  }
  v <- terra::vect(dom)
  g <- terra::rast(v, resolution = resolution)
  terra::values(g) <- 1
  g <- terra::mask(g, v)
  names(g) <- "grid"
  g
}
