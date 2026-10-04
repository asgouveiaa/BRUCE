#' Attach spatial covariates
#'
#' Extracts covariates at station locations and at background cells,
#' standardises them, and keeps the rasters for prediction.
#'
#' @param d A `bruce_data` object.
#' @param rasters A [terra::SpatRaster] whose layer names are the covariate
#'   names (syntactic R names), e.g. `elev`, `urb`, `svf`.
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
#'
#' @details Covariate values that are missing after extraction are set to the
#'   covariate mean (0 on the standardised scale). The observed range at
#'   stations and background cells is stored and used by [bruce_predict()] to
#'   cap covariates on the prediction grid, so the model never extrapolates a
#'   linear effect beyond the data.
#'
#' @return `d` with a `covariates` element.
#' @export
bruce_covariates <- function(d, rasters, background_method = "buffer_mean",
                             background_buffer = 5500, scale_on = NULL) {
  stopifnot(inherits(d, "bruce_data"))
  if (!inherits(rasters, "SpatRaster")) stop("`rasters` must be a terra SpatRaster", call. = FALSE)
  nms <- names(rasters)
  if (!identical(make.names(nms), nms)) stop("Layer names must be syntactic R names", call. = FALSE)
  reserved <- c("intercept", "day_idx", "station_idx", "s", "y")
  if (any(nms %in% reserved) || any(startsWith(nms, "bias_"))) {
    stop("Layer names cannot be any of: ", paste(reserved, collapse = ", "),
         " or start with 'bias_'", call. = FALSE)
  }
  if (!terra::same.crs(rasters, d$crs$wkt)) rasters <- terra::project(rasters, d$crs$wkt)

  expand <- function(arg, default) {
    if (is.null(arg)) return(stats::setNames(rep(default, length(nms)), nms))
    if (is.null(names(arg))) return(stats::setNames(rep(arg, length.out = length(nms)), nms))
    out <- stats::setNames(rep(default, length(nms)), nms)
    out[names(arg)] <- arg
    out
  }
  bg_method <- expand(background_method, "buffer_mean")
  if (!all(bg_method %in% c("buffer_mean", "point", "zero"))) {
    stop("`background_method` must be 'buffer_mean', 'point' or 'zero'", call. = FALSE)
  }
  has_bg <- !is.null(d$background)
  default_scale <- ifelse(has_bg & bg_method != "zero", "background", "stations")
  sc_on <- if (is.null(scale_on)) default_scale else expand(scale_on, NA)
  sc_on[is.na(sc_on)] <- default_scale[is.na(sc_on)]

  vect_pts <- function(xy) terra::vect(as.matrix(xy), crs = d$crs$wkt)

  # ---- stations -------------------------------------------------------------
  st <- unique(d$obs[, c("station_id", "x", "y")])
  raw_st <- terra::extract(rasters, vect_pts(st[, c("x", "y")]))
  raw_st <- raw_st[, nms, drop = FALSE]

  # ---- background cells -----------------------------------------------------
  raw_bg <- NULL
  if (has_bg) {
    bg_loc <- unique(d$background[, c("x", "y")])
    raw_bg <- as.data.frame(matrix(NA_real_, nrow(bg_loc), length(nms), dimnames = list(NULL, nms)))
    pts <- vect_pts(bg_loc)
    if (any(bg_method == "buffer_mean")) {
      lyr <- nms[bg_method == "buffer_mean"]
      buf <- terra::buffer(pts, width = background_buffer)
      ex <- terra::extract(rasters[[lyr]], buf, fun = mean, na.rm = TRUE)
      raw_bg[, lyr] <- ex[, lyr]
    }
    if (any(bg_method == "point")) {
      lyr <- nms[bg_method == "point"]
      ex <- terra::extract(rasters[[lyr]], pts)
      raw_bg[, lyr] <- ex[, lyr]
    }
  }

  # ---- standardisation ------------------------------------------------------
  scaling <- data.frame(name = nms, mean = NA_real_, sd = NA_real_, min = NA_real_,
                        max = NA_real_, scale_on = unname(sc_on),
                        background_method = unname(bg_method), stringsAsFactors = FALSE)
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

  d$covariates <- list(names = nms, rasters = terra::wrap(rasters), scaling = scaling,
                       station = st_out, background = bg_out)
  d
}

# Standardise (and optionally cap) raw covariate values using stored scaling.
scale_covariates <- function(raw, scaling, cap = TRUE) {
  out <- raw
  for (k in seq_len(nrow(scaling))) {
    cn <- scaling$name[k]
    v <- raw[[cn]]
    if (cap) v <- pmin(pmax(v, scaling$min[k]), scaling$max[k])
    z <- (v - scaling$mean[k]) / scaling$sd[k]
    z[is.na(z)] <- 0
    out[[cn]] <- z
  }
  out
}
