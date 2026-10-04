#' Simulate a synthetic city
#'
#' Generates data from the BRUCE data-generating process on an invented city,
#' for testing the pipeline and checking parameter recovery: two covariates
#' (`elev`, a pair of hills; `urb`, an urban core), a reference network, a
#' biased low-cost network, an optional coarse background product, a
#' Matern field with AR(1) dependence between days and a city-wide daily
#' AR(1) term.
#'
#' @param n_reference,n_biased Number of reference and low-cost stations.
#' @param n_days Number of days.
#' @param extent Side of the square domain (metres).
#' @param range,sigma Matern range (metres) and standard deviation of the
#'   field (smoothness 1).
#' @param group_rho Day-to-day correlation of the field.
#' @param intercept,beta_elev,beta_urb Mean level and covariate effects
#'   (`beta_elev` per metre of elevation, `beta_urb` per unit urban fraction).
#' @param day_sd,day_rho Standard deviation and autocorrelation of the daily
#'   term.
#' @param bias Offset of the low-cost network.
#' @param sd_reference,sd_biased Observation noise of each network.
#' @param sd_station Standard deviation of persistent per-station offsets in
#'   the low-cost network.
#' @param background Simulate a coarse background product.
#' @param background_res,background_sd Its cell size (metres) and noise.
#' @param res Covariate raster resolution (metres).
#' @param crs EPSG code used for the coordinates.
#' @param origin Lower-left corner of the domain in `crs`.
#' @param clustered If `TRUE`, low-cost stations concentrate in the urban
#'   core, mimicking the deployment pattern that drove field collapse in the
#'   original application.
#' @param seed Random seed.
#'
#' @return A list with `obs`, `stations`, `sources`, `background`,
#'   `covariates` (SpatRaster), `domain` (sf polygon), `crs`, `truth`
#'   (SpatRaster, one layer per day, latent temperature) and `params`. The
#'   first four plug straight into [bruce_data()].
#' @export
bruce_simulate <- function(n_reference = 12, n_biased = 40, n_days = 10, extent = 30000,
                           range = 8000, sigma = 1, group_rho = 0.6,
                           intercept = 25, beta_elev = -0.0065, beta_urb = 1.5,
                           day_sd = 1.5, day_rho = 0.7, bias = 1.5,
                           sd_reference = 0.3, sd_biased = 1, sd_station = 0.5,
                           background = TRUE, background_res = 9000, background_sd = 0.7,
                           res = 250, crs = 31983, origin = c(650000, 7450000),
                           clustered = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  crs_str <- paste0("EPSG:", crs)

  # ---- covariates -----------------------------------------------------------
  r <- terra::rast(xmin = origin[1], xmax = origin[1] + extent,
                   ymin = origin[2], ymax = origin[2] + extent,
                   resolution = res, crs = crs_str)
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  u <- (xy[, 1] - origin[1]) / extent
  v <- (xy[, 2] - origin[2]) / extent
  elev <- 800 * exp(-((u - 0.25)^2 + (v - 0.7)^2) / 0.02) +
          400 * exp(-((u - 0.75)^2 + (v - 0.3)^2) / 0.03)
  urb <- stats::plogis(12 * (0.3 - sqrt((u - 0.55)^2 + (v - 0.45)^2)))
  cov <- c(r, r)
  terra::values(cov) <- cbind(elev, urb)
  names(cov) <- c("elev", "urb")

  # ---- spatial field on a node lattice, AR(1) in time -------------------------
  n_side <- 30
  gx <- seq(0, extent, length.out = n_side)
  nodes <- as.matrix(expand.grid(gx, gx))
  D <- as.matrix(stats::dist(nodes))
  kappa <- sqrt(8) / range
  C <- matern1(D, kappa) * sigma^2
  L <- t(chol(C + diag(1e-8, nrow(C))))
  W <- matrix(0, nrow(nodes), n_days)
  W[, 1] <- L %*% stats::rnorm(nrow(nodes))
  for (t in seq_len(n_days)[-1]) {
    W[, t] <- group_rho * W[, t - 1] + sqrt(1 - group_rho^2) * (L %*% stats::rnorm(nrow(nodes)))
  }
  day <- numeric(n_days)
  day[1] <- stats::rnorm(1, 0, day_sd)
  for (t in seq_len(n_days)[-1]) {
    day[t] <- day_rho * day[t - 1] + sqrt(1 - day_rho^2) * stats::rnorm(1, 0, day_sd)
  }
  node_of <- function(px, py) {
    ix <- pmin(pmax(round((px - origin[1]) / extent * (n_side - 1)) + 1, 1), n_side)
    iy <- pmin(pmax(round((py - origin[2]) / extent * (n_side - 1)) + 1, 1), n_side)
    (iy - 1) * n_side + ix
  }
  truth_fun <- function(px, py, e, ub, t) {
    intercept + beta_elev * e + beta_urb * ub + day[t] + W[node_of(px, py), t]
  }

  dates <- as.Date("2025-01-01") + seq_len(n_days) - 1
  node_cells <- node_of(xy[, 1], xy[, 2])
  truth_m <- sapply(seq_len(n_days), function(t)
    intercept + beta_elev * elev + beta_urb * urb + day[t] + W[node_cells, t])
  truth <- terra::rast(r, nlyrs = n_days)
  terra::values(truth) <- truth_m
  names(truth) <- as.character(dates)

  # ---- stations ---------------------------------------------------------------
  place <- function(n, cluster) {
    if (!cluster) {
      return(cbind(origin[1] + stats::runif(n, 0.05, 0.95) * extent,
                   origin[2] + stats::runif(n, 0.05, 0.95) * extent))
    }
    cbind(origin[1] + pmin(pmax(stats::rnorm(n, 0.55, 0.08), 0.02), 0.98) * extent,
          origin[2] + pmin(pmax(stats::rnorm(n, 0.45, 0.08), 0.02), 0.98) * extent)
  }
  p_ref <- place(n_reference, FALSE)
  p_bia <- place(n_biased, clustered)
  stations <- data.frame(
    station_id = c(sprintf("REF%02d", seq_len(n_reference)), sprintf("IOT%03d", seq_len(n_biased))),
    source = rep(c("reference", "lowcost"), c(n_reference, n_biased)),
    x = c(p_ref[, 1], p_bia[, 1]), y = c(p_ref[, 2], p_bia[, 2]),
    stringsAsFactors = FALSE)
  cv_st <- terra::extract(cov, terra::vect(as.matrix(stations[, c("x", "y")]), crs = crs_str))
  st_off <- ifelse(stations$source == "lowcost", stats::rnorm(nrow(stations), 0, sd_station), 0)

  obs <- do.call(rbind, lapply(seq_len(n_days), function(t) {
    mu <- truth_fun(stations$x, stations$y, cv_st$elev, cv_st$urb, t)
    noise <- ifelse(stations$source == "reference",
                    stats::rnorm(nrow(stations), 0, sd_reference),
                    bias + st_off + stats::rnorm(nrow(stations), 0, sd_biased))
    data.frame(station_id = stations$station_id, date = dates[t], value = mu + noise,
               truth = mu, stringsAsFactors = FALSE)
  }))

  sources <- data.frame(source = c("reference", "lowcost"), bias = c(FALSE, TRUE),
                        obs_sd = c(sd_reference, NA))

  # ---- background: block means of the truth plus noise ------------------------
  bg <- NULL
  if (background) {
    fact <- max(1, round(background_res / res))
    agg <- terra::aggregate(truth, fact = fact, fun = "mean")
    bxy <- terra::xyFromCell(agg, seq_len(terra::ncell(agg)))
    bv <- terra::values(agg, mat = TRUE)
    bg <- do.call(rbind, lapply(seq_len(n_days), function(t)
      data.frame(x = bxy[, 1], y = bxy[, 2], date = dates[t],
                 value = bv[, t] + stats::rnorm(nrow(bxy), 0, background_sd))))
    bg <- bg[!is.na(bg$value), ]
  }

  domain <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(c(
    xmin = origin[1], ymin = origin[2], xmax = origin[1] + extent, ymax = origin[2] + extent),
    crs = sf::st_crs(crs))))

  list(obs = obs, stations = stations, sources = sources, background = bg,
       covariates = cov, domain = domain, crs = crs, truth = truth,
       params = list(intercept = intercept, beta_elev = beta_elev, beta_urb = beta_urb,
                     bias = bias, range = range, sigma = sigma, group_rho = group_rho,
                     day_sd = day_sd, day_rho = day_rho, sd_reference = sd_reference,
                     sd_biased = sd_biased, sd_station = sd_station,
                     background_sd = background_sd))
}

# Matern covariance with smoothness 1 (alpha = 2 in SPDE terms), unit variance.
matern1 <- function(d, kappa) {
  kd <- kappa * d
  out <- kd * besselK(pmax(kd, 1e-10), 1)
  out[d == 0] <- 1
  out
}
