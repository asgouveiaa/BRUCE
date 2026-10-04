#' Attach spatial covariates
#'
#' Extracts covariates at station locations and at background cells,
#' standardises them, and keeps the rasters for prediction.
#'
#' @param d A `bruce_data` object.
#' @param rasters Either a [terra::SpatRaster] whose layer names are the
#'   covariate names, or a **named list** of single-layer SpatRasters. With a
#'   list, each covariate keeps its own projection and resolution and is
#'   extracted directly from its native raster, with no resampling to a
#'   common grid. This matters for fine-scale covariates such as sky view
#'   factor, where resampling can pick a neighbouring cell. Names must be
#'   syntactic R names, e.g. `elev`, `urb`, `svf`.
#' @param background_method How each covariate is represented at background
#'   cells: `"buffer_mean"` (mean within `background_buffer` metres, matching
#'   the footprint of a coarse grid cell), `"point"` (value at the cell
#'   centre), or `"zero"` (set to the covariate mean; for covariates that have
#'   no meaning at the background's scale, such as sky view factor). A single
#'   value or a named vector per covariate.
#' @param background_buffer Buffer radius in metres for `"buffer_mean"`.
#' @param scale_on Which locations define the standardisation of each
#'   covariate: `"background"` or `"stations"`. Defaults to `"background"`
#'   when a background exists and the covariate is not `"zero"`, otherwise
#'   `"stations"`.
#' @param na_value Raw value used where extraction returns `NA` (outside the
#'   raster or over missing cells), as a named vector per covariate, e.g.
#'   `c(urb = 0)` for an urban fraction raster that only covers the city.
#'   Applied before standardisation, at stations, background cells and the
#'   prediction grid. Covariates not listed get the covariate mean (0 on the
#'   standardised scale).
#'
#' @details The observed range at stations and background cells is stored
#'   and used by [bruce_predict()] to cap covariates on the prediction grid,
#'   so the model never extrapolates a linear effect beyond the data.
#'
#' @return `d` with a `covariates` element.
#' @export
bruce_covariates <- function(d, rasters, background_method = "buffer_mean",
                             background_buffer = 5500, scale_on = NULL, na_value = NULL) {
  stopifnot(inherits(d, "bruce_data"))
  rlist <- as_raster_list(rasters)
  nms <- names(rlist)
  reserved <- c("intercept", "day_idx", "station_idx", "s", "y")
  if (any(nms %in% reserved) || any(startsWith(nms, "bias_"))) {
    stop("Covariate names cannot be any of: ", paste(reserved, collapse = ", "),
         " or start with 'bias_'", call. = FALSE)
  }

  expand <- function(arg, default) {
    if (is.null(arg)) return(stats::setNames(rep(default, length(nms)), nms))
    if (is.null(names(arg))) return(stats::setNames(rep(arg, length.out = length(nms)), nms))
    bad <- setdiff(names(arg), nms)
    if (length(bad)) stop("Unknown covariate name(s): ", paste(bad, collapse = ", "), call. = FALSE)
    out <- stats::setNames(rep(default, length(nms)), nms)
    out[names(arg)] <- arg
    out
  }
  bg_method <- expand(background_method, "buffer_mean")
  if (!all(bg_method %in% c("buffer_mean", "point", "zero"))) {
    stop("`background_method` must be 'buffer_mean', 'point' or 'zero'", call. = FALSE)
  }
  na_val <- expand(na_value, NA_real_)
  has_bg <- !is.null(d$background)
  default_scale <- ifelse(has_bg & bg_method != "zero", "background", "stations")
  sc_on <- if (is.null(scale_on)) default_scale else expand(scale_on, NA)
  sc_on[is.na(sc_on)] <- default_scale[is.na(sc_on)]

  fill <- function(v, cn) { if (!is.na(na_val[[cn]])) v[is.na(v)] <- na_val[[cn]]; v }

  # ---- extraction, each covariate from its own raster --------------------------
  st <- unique(d$obs[, c("station_id", "x", "y")])
  st <- st[!duplicated(st$station_id), ]
  raw_st <- as.data.frame(lapply(stats::setNames(nms, nms), function(cn)
    fill(extract_points(rlist[[cn]], st[, c("x", "y")], d$crs), cn)))

  raw_bg <- NULL
  if (has_bg) {
    bg_loc <- unique(d$background[, c("x", "y")])
    raw_bg <- as.data.frame(lapply(stats::setNames(nms, nms), function(cn) {
      v <- switch(bg_method[[cn]],
                  buffer_mean = extract_buffer(rlist[[cn]], bg_loc, d$crs, background_buffer),
                  point = extract_points(rlist[[cn]], bg_loc, d$crs),
                  zero = rep(NA_real_, nrow(bg_loc)))
      if (bg_method[[cn]] == "zero") v else fill(v, cn)
    }))
  }

  # ---- standardisation ------------------------------------------------------
  scaling <- data.frame(name = nms, mean = NA_real_, sd = NA_real_, min = NA_real_,
                        max = NA_real_, scale_on = unname(sc_on),
                        background_method = unname(bg_method),
                        na_value = unname(na_val), stringsAsFactors = FALSE)
  st_out <- data.frame(station_id = st$station_id)
  bg_out <- if (has_bg) data.frame(x = bg_loc$x, y = bg_loc$y) else NULL
  n_na <- c(stations = 0, background = 0)

  for (k in seq_along(nms)) {
    cn <- nms[k]
    ref <- if (sc_on[cn] == "background" && has_bg) raw_bg[[cn]] else raw_st[[cn]]
    mu <- mean(ref, na.rm = TRUE)
    s <- stats::sd(ref, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-6) s <- 1
    pool <- c(raw_st[[cn]], if (has_bg && bg_method[cn] != "zero") raw_bg[[cn]])
    scaling[k, c("mean", "sd", "min", "max")] <- c(mu, s, range(pool, na.rm = TRUE))

    z <- (raw_st[[cn]] - mu) / s
    n_na["stations"] <- n_na["stations"] + sum(is.na(z))
    z[is.na(z)] <- 0
    st_out[[cn]] <- z
    if (has_bg) {
      zb <- if (bg_method[cn] == "zero") rep(0, nrow(bg_loc)) else (raw_bg[[cn]] - mu) / s
      n_na["background"] <- n_na["background"] + sum(is.na(zb))
      zb[is.na(zb)] <- 0
      bg_out[[cn]] <- zb
    }
  }
  if (any(n_na > 0)) {
    message(sprintf("Missing covariate values set to the mean: %d at stations, %d at background cells",
                    n_na["stations"], n_na["background"]))
  }

  d$covariates <- list(names = nms, rasters = lapply(rlist, terra::wrap), scaling = scaling,
                       station = st_out, background = bg_out)
  d
}

# A SpatRaster or a named list of SpatRasters -> named list of single layers.
as_raster_list <- function(rasters) {
  if (inherits(rasters, "SpatRaster")) {
    rl <- lapply(names(rasters), function(n) rasters[[n]])
    names(rl) <- names(rasters)
  } else if (is.list(rasters) && length(rasters) && all(vapply(rasters, inherits, TRUE, "SpatRaster"))) {
    if (is.null(names(rasters)) || any(names(rasters) == "")) {
      stop("A list of rasters must be named (one name per covariate)", call. = FALSE)
    }
    rl <- lapply(rasters, function(r) r[[1]])
  } else {
    stop("`rasters` must be a SpatRaster or a named list of SpatRasters", call. = FALSE)
  }
  nms <- names(rl)
  if (!identical(make.names(nms), nms) || anyDuplicated(nms)) {
    stop("Covariate names must be unique syntactic R names", call. = FALSE)
  }
  for (n in nms) names(rl[[n]]) <- n
  rl
}

# Raster values at points given in the data CRS (projected to the raster's).
extract_points <- function(r, xy, crs) {
  pts <- terra::vect(as.matrix(xy[, c("x", "y")]), crs = crs$wkt)
  if (!terra::same.crs(r, pts)) pts <- terra::project(pts, terra::crs(r))
  terra::extract(r, pts, ID = FALSE)[[1]]
}

# Mean raster value within a buffer (metres, in the data CRS) around points.
extract_buffer <- function(r, xy, crs, width) {
  pts <- terra::vect(as.matrix(xy[, c("x", "y")]), crs = crs$wkt)
  buf <- terra::buffer(pts, width = width)
  if (!terra::same.crs(r, buf)) buf <- terra::project(buf, terra::crs(r))
  v <- terra::extract(r, buf, fun = mean, na.rm = TRUE, ID = FALSE)[[1]]
  v[!is.finite(v)] <- NA
  v
}

# Standardise (and optionally cap) raw covariate values using stored scaling.
scale_covariates <- function(raw, scaling, cap = TRUE) {
  out <- raw
  for (k in seq_len(nrow(scaling))) {
    cn <- scaling$name[k]
    v <- raw[[cn]]
    nv <- if ("na_value" %in% names(scaling)) scaling$na_value[k] else NA
    if (!is.na(nv)) v[is.na(v)] <- nv
    if (cap) v <- pmin(pmax(v, scaling$min[k]), scaling$max[k])
    z <- (v - scaling$mean[k]) / scaling$sd[k]
    z[is.na(z)] <- 0
    out[[cn]] <- z
  }
  out
}
