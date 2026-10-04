# Internal helpers -------------------------------------------------------------

check_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required for this function. Install it with:\n",
      "  install.packages('INLA', repos = c(getOption('repos'), ",
      "INLA = 'https://inla.r-inla-download.org/R/stable'), dep = TRUE)",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

need_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop(sprintf("`%s` is missing column(s): %s", what, paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

# Adds projected x/y (metres) to a data frame that has either x/y already in
# `crs` or lon/lat in WGS84.
project_xy <- function(df, crs, what) {
  if (all(c("x", "y") %in% names(df))) return(df)
  if (!all(c("lon", "lat") %in% names(df))) {
    stop(sprintf("`%s` needs either x/y (in `crs`) or lon/lat (WGS84) columns", what),
         call. = FALSE)
  }
  pts <- sf::st_as_sf(df[, c("lon", "lat")], coords = c("lon", "lat"), crs = 4326)
  xy <- sf::st_coordinates(sf::st_transform(pts, crs))
  df$x <- xy[, 1]
  df$y <- xy[, 2]
  df
}

# Median nearest-neighbour distance between distinct locations.
nn_median <- function(xy) {
  xy <- unique(as.matrix(xy))
  if (nrow(xy) < 2) return(NA_real_)
  if (nrow(xy) > 3000) xy <- xy[sample(nrow(xy), 3000), ]
  dm <- as.matrix(stats::dist(xy))
  diag(dm) <- Inf
  stats::median(apply(dm, 1, min))
}

# CRPS of an ensemble (energy form, O(n log n)).
crps_ensemble <- function(y, x) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n == 0) return(NA_real_)
  term1 <- mean(abs(x - y))
  term2 <- 2 / n^2 * sum((2 * seq_len(n) - n - 1) * x)  # E|X - X'|
  term1 - 0.5 * term2
}

# Internal precision of an AR(1) correlation on INLA's scale.
rho_to_theta <- function(rho) log((1 + rho) / (1 - rho))

window_days <- function(d, window) {
  if (is.null(window)) {
    rng <- range(d$obs$date)
  } else {
    rng <- as.Date(window)
    if (length(rng) != 2 || anyNA(rng)) stop("`window` must be two dates", call. = FALSE)
  }
  days <- seq(rng[1], rng[2], by = "day")
  if (length(days) > 62) {
    message(sprintf(
      "Window spans %d days; the model was developed for monthly windows and cost grows with window length.",
      length(days)))
  }
  days
}

`%||%` <- function(a, b) if (is.null(a)) b else a
