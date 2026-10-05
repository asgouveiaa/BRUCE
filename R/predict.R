#' Predict on a raster grid
#'
#' Predicts the variable at the reference-network scale (no network offset,
#' no station effect) for each day of the fitted window.
#'
#' @param fit A `bruce_fit`.
#' @param grid A [terra::SpatRaster] defining the prediction grid. Cells
#'   where its first layer is `NA` are skipped, so a masked raster works as a
#'   domain. `NULL` uses the covariate rasters' own grid.
#' @param domain Optional polygon (`sf` or `SpatVector`) to mask the grid.
#' @param dates Subset of the fitted days; `NULL` predicts all.
#' @param n_samples Number of joint posterior samples. `0` returns posterior
#'   means only (fast, from INLA's marginal summaries). `> 0` returns mean,
#'   standard deviation and quantiles from samples; requires the fit to have
#'   been made with `keep_config = TRUE`.
#' @param probs Quantiles to return when sampling.
#' @param keep_samples Keep the per-cell samples (needed by [bruce_zonal()]
#'   for uncertainty of zonal means). Memory grows as cells x samples x days.
#' @param cap Cap covariates to the range observed at stations and background
#'   cells (avoids extrapolating linear effects).
#' @param resample How each covariate is brought to the grid: `"average"`
#'   (block mean of the native raster over each grid cell, via
#'   [terra::project()]), `"bilinear"`, or `"near"` (the native raster's value
#'   at the grid cell centre, read directly from that raster in its own
#'   projection). A single value or a named vector, e.g.
#'   `c(elev = "average", svf = "near")`. Rasters that already share the
#'   grid's geometry are used as they are.
#' @param seed Seed for posterior sampling. A number makes samples
#'   reproducible (INLA then samples on one thread); `NULL` samples in
#'   parallel, faster but not reproducible.
#'
#' @return An object of class `bruce_prediction` with `mean`, `sd` and
#'   `quantiles` (SpatRasters with one layer per day) and, optionally,
#'   `samples`.
#' @export
bruce_predict <- function(fit, grid = NULL, domain = NULL, dates = NULL, n_samples = 0,
                          probs = c(0.05, 0.95), keep_samples = FALSE, cap = TRUE,
                          resample = "average", seed = 1) {
  stopifnot(inherits(fit, "bruce_fit"))
  check_inla()
  days <- fit$days
  if (!is.null(dates)) {
    dates <- as.Date(dates)
    if (!all(dates %in% days)) stop("Some `dates` are outside the fitted window", call. = FALSE)
  } else dates <- days
  g_idx <- match(dates, days)

  gp <- prepare_grid(fit, grid, domain, cap, resample)
  comps <- latent_components(fit, n_samples, seed)
  S <- ncol(comps$fixed)

  tmpl <- gp$template
  nc <- terra::ncell(tmpl)
  empty <- function() matrix(NA_real_, nc, length(dates))
  m_mean <- empty()
  m_sd <- if (S > 1) empty() else NULL
  m_q <- if (S > 1) lapply(probs, function(p) empty()) else NULL
  samples <- if (keep_samples && S > 1) vector("list", length(dates)) else NULL

  for (k in seq_along(dates)) {
    P <- eval_predictor(fit, comps, gp$X, gp$A, g_idx[k])
    m_mean[gp$cells, k] <- rowMeans(P)
    if (S > 1) {
      m_sd[gp$cells, k] <- row_sd(P)
      qs <- row_quantiles(P, probs)
      for (j in seq_along(probs)) m_q[[j]][gp$cells, k] <- qs[, j]
      if (!is.null(samples)) samples[[k]] <- P
    }
  }
  if (!is.null(samples)) names(samples) <- as.character(dates)

  to_rast <- function(m) {
    r <- terra::rast(tmpl, nlyrs = length(dates))
    terra::values(r) <- m
    names(r) <- as.character(dates)
    r
  }
  r_mean <- to_rast(m_mean)
  r_sd <- if (S > 1) to_rast(m_sd) else NULL
  r_q <- if (S > 1) stats::setNames(lapply(m_q, to_rast), paste0("q", probs)) else NULL

  structure(list(mean = r_mean, sd = r_sd, quantiles = r_q, samples = samples,
                 cells = gp$cells, template = tmpl, dates = dates, n_samples = S),
            class = "bruce_prediction")
}

prepare_grid <- function(fit, grid, domain, cap, resample = "average") {
  cov_names <- fit$cov_names
  rasters <- if (length(cov_names)) unwrap_covariates(fit$data$covariates$rasters) else NULL
  if (is.null(grid)) {
    if (is.null(rasters)) stop("`grid` is required when the model has no covariates", call. = FALSE)
    grid <- rasters[[1]]
  }
  if (!terra::same.crs(grid, fit$data$crs$wkt)) grid <- terra::project(grid, fit$data$crs$wkt)
  tmpl <- terra::rast(grid[[1]])
  ok <- if (terra::hasValues(grid)) !is.na(terra::values(grid[[1]], mat = FALSE)) else rep(TRUE, terra::ncell(tmpl))
  if (!is.null(domain)) {
    dv <- if (inherits(domain, "SpatVector")) domain else terra::vect(sf::st_transform(sf::st_as_sf(domain), fit$data$crs))
    inside <- terra::rasterize(dv, tmpl, field = 1, background = NA)
    ok <- ok & !is.na(terra::values(inside, mat = FALSE))
  }
  cells <- which(ok)
  if (!length(cells)) stop("No grid cells to predict", call. = FALSE)
  xy <- terra::xyFromCell(tmpl, cells)

  X <- matrix(1, length(cells), 1, dimnames = list(NULL, "intercept"))
  if (length(cov_names)) {
    meth <- stats::setNames(rep("average", length(cov_names)), cov_names)
    if (is.null(names(resample))) meth[] <- resample[1] else meth[names(resample)] <- resample
    if (!all(meth %in% c("average", "bilinear", "near"))) {
      stop("`resample` must be 'average', 'bilinear' or 'near'", call. = FALSE)
    }
    raw <- as.data.frame(lapply(stats::setNames(cov_names, cov_names), function(cn) {
      r1 <- rasters[[cn]]
      if (terra::same.crs(r1, tmpl) && terra::compareGeom(r1, tmpl, stopOnError = FALSE)) {
        return(terra::values(r1, mat = FALSE)[cells])
      }
      if (meth[[cn]] == "near") {
        return(extract_points(r1, data.frame(x = xy[, 1], y = xy[, 2]), fit$data$crs))
      }
      terra::values(terra::project(r1, tmpl, method = meth[[cn]]), mat = FALSE)[cells]
    }))
    n_na <- sum(!stats::complete.cases(raw))
    if (n_na) message(sprintf("%d grid cell(s) with missing covariates (filled with `na_value` or the covariate mean)", n_na))
    z <- scale_covariates(raw, fit$data$covariates$scaling, cap = cap)
    X <- cbind(X, as.matrix(z[, cov_names, drop = FALSE]))
  }
  A <- INLA::inla.spde.make.A(fit$mesh, loc = xy)
  list(template = tmpl, cells = cells, X = X, A = A)
}

# Stored covariate rasters -> named list of SpatRasters (also reads fits made
# with 0.1.0, which stored one packed multi-layer raster).
unwrap_covariates <- function(x) {
  if (inherits(x, c("PackedSpatRaster", "SpatRaster"))) return(as_raster_list(terra::unwrap(x)))
  lapply(x, terra::unwrap)
}

# Posterior means (S = 1) or joint samples (S = n_samples) of the latent
# components needed for prediction.
latent_components <- function(fit, n_samples, seed) {
  n_days <- length(fit$days)
  coef_names <- c("intercept", fit$cov_names)
  if (n_samples <= 0) {
    fx <- fit$inla$summary.fixed
    fixed <- matrix(fx[coef_names, "mean"], ncol = 1)
    day <- if (fit$spec$day_effect) matrix(fit$inla$summary.random$day_idx$mean, ncol = 1)
           else matrix(0, n_days, 1)
    field <- matrix(fit$inla$summary.random$s$mean, ncol = 1)
    return(list(fixed = fixed, day = day, field = field))
  }
  if (is.null(fit$inla$misc$configs)) {
    stop("Posterior sampling needs a fit made with `keep_config = TRUE`", call. = FALSE)
  }
  # A non-zero seed makes INLA sample serially (reproducible) and warn every
  # time; seed = NULL uses seed 0, parallel but not reproducible.
  smp <- withCallingHandlers(
    INLA::inla.posterior.sample(n_samples, fit$inla, seed = seed %||% 0L),
    warning = function(w) {
      if (grepl("parallel model is disabled", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  nm <- rownames(smp[[1]]$latent)
  i_fix <- match(paste0(coef_names, ":1"), nm)
  if (anyNA(i_fix)) stop("Could not find fixed effects in posterior samples", call. = FALSE)
  i_day <- grep("^day_idx:", nm)
  i_s <- grep("^s:", nm)
  take <- function(i) do.call(cbind, lapply(smp, function(s) s$latent[i, 1]))
  list(fixed = take(i_fix),
       day = if (fit$spec$day_effect) take(i_day) else matrix(0, n_days, n_samples),
       field = take(i_s))
}

# Linear predictor for rows of X / A on fitted day g: matrix rows x samples.
eval_predictor <- function(fit, comps, X, A, g) {
  n <- fit$spde$n.spde
  rows <- ((g - 1) * n + 1):(g * n)
  S <- ncol(comps$fixed)
  P <- X %*% comps$fixed
  P <- P + matrix(comps$day[g, ], nrow(X), S, byrow = TRUE)
  P <- P + as.matrix(A %*% comps$field[rows, , drop = FALSE])
  P
}
