#' Build the SPDE mesh
#'
#' Triangulates around all observation locations (stations and background
#' cells) with a non-convex hull boundary.
#'
#' @param d A `bruce_data` object.
#' @param max_edge Length-2 vector: maximum triangle edge inside and outside
#'   the hull (metres).
#' @param cutoff Minimum distance between mesh vertices (metres).
#' @param offset Length-2 vector: inner and outer extension (metres).
#' @param convex Convexity parameter of the non-convex hull (metres if
#'   positive, fraction of the domain diameter if negative). The hull is built
#'   with the same algorithm as `INLA::inla.nonconvex.hull()`, so meshes match
#'   those built with `INLA::inla.mesh.2d()` from the same points.
#'
#' @details When left `NULL`, sizes are derived from the data: the median
#'   nearest-neighbour distance between stations (`nn`) and the diagonal of
#'   the bounding box (`L`). The inner edge is `clamp(2 * nn, L / 40, L / 10)`,
#'   the outer edge 8 times that, the cutoff 0.6 times the inner edge, and the
#'   offsets `inner / 5` and `L / 2`. These are heuristics; for a
#'   city-scale network they land in the same range as the values used in the
#'   original Rio de Janeiro application (`max_edge = c(5000, 40000)`,
#'   `cutoff = 3000`, `offset = c(1000, 50000)`), but check the number of
#'   vertices (`mesh$n`) before fitting a full window.
#'
#' @return An `fm_mesh_2d` object with attribute `bruce_nn` (the station
#'   nearest-neighbour median, used by the default range prior).
#' @export
bruce_mesh <- function(d, max_edge = NULL, cutoff = NULL, offset = NULL, convex = -0.1) {
  stopifnot(inherits(d, "bruce_data"))
  st_xy <- unique(as.matrix(d$obs[, c("x", "y")]))
  loc <- st_xy
  if (!is.null(d$background)) loc <- rbind(loc, unique(as.matrix(d$background[, c("x", "y")])))
  loc <- unique(loc)

  nn <- nn_median(st_xy)
  L <- sqrt(diff(range(loc[, 1]))^2 + diff(range(loc[, 2]))^2)
  inner <- min(max(2 * nn, L / 40, na.rm = TRUE), L / 10)
  max_edge <- max_edge %||% c(inner, 8 * inner)
  cutoff <- cutoff %||% (0.6 * max_edge[1])
  offset <- offset %||% c(max_edge[1] / 5, L / 2)

  bnd <- inla_style_hull(loc, convex)
  mesh <- fmesher::fm_mesh_2d(loc = loc, boundary = bnd, max.edge = max_edge,
                              cutoff = cutoff, offset = offset)
  attr(mesh, "bruce_nn") <- nn
  attr(mesh, "bruce_settings") <- list(max_edge = max_edge, cutoff = cutoff, offset = offset)
  message(sprintf("Mesh: %d vertices (max_edge = %s, cutoff = %.0f m; station NN median = %.0f m)",
                  mesh$n, paste(round(max_edge), collapse = "/"), cutoff, nn))
  mesh
}

# Non-convex hull with INLA's algorithm (what INLA::inla.nonconvex.hull()
# uses). Older fmesher versions expose it as fm_nonconvex_hull_inla(); newer
# ones as fm_nonconvex_hull(method = "fm"). Plain fm_nonconvex_hull() in the
# older versions builds a different (sf buffer) hull, which changes the mesh.
inla_style_hull <- function(loc, convex) {
  if ("method" %in% names(formals(fmesher::fm_nonconvex_hull))) {
    return(fmesher::fm_nonconvex_hull(loc, convex = convex, method = "fm", format = "fm"))
  }
  fmesher::fm_nonconvex_hull_inla(loc, convex = convex)
}
