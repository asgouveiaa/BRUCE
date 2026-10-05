#' Cross-validation by station or by region
#'
#' Refits the model with stations held out and scores predictions at the
#' held-out reference stations.
#'
#' @param d,mesh,spec,window As in [bruce_fit()].
#' @param scheme `"lso"` (leave-station-out: one fold per test station,
#'   removing only that station) or `"lro"` (leave-region-out: one fold per
#'   region, removing every station of every network in the region; the
#'   background stays).
#' @param regions For `"lro"`: polygons (`sf`) with an id column.
#' @param region_col Name of the id column in `regions`.
#' @param test_sources Networks scored at held-out locations; default all
#'   reference networks (`bias = FALSE`).
#' @param n_samples Posterior samples per fold. With `0`, scores use
#'   posterior means and CRPS/coverage are not computed.
#' @param level Nominal level of the predictive interval for coverage.
#' @param noise_sd Observation noise added to samples before scoring (`0`
#'   scores the latent field, as in the original application; use the
#'   reference network's `obs_sd` to score the predictive distribution of an
#'   observation).
#' @param checkpoint_dir Optional directory; each finished fold is saved
#'   there and skipped on rerun, so long runs can be resumed.
#' @param theta_init Passed to [bruce_fit()]. Starting each fold from the
#'   mode of the full-data fit (`full_fit$inla$mode$theta`) usually saves
#'   most of the mode search, since removing a few stations barely moves it.
#' @param fix_theta Keep the hyperparameters fixed at `theta_init` in every
#'   fold (see [bruce_fit()]): each fold becomes a single Gaussian
#'   approximation, much faster, at the cost of ignoring how hyperparameters
#'   would change without the held-out data. With `theta_init = NULL` the
#'   full data are fitted once first and their mode is used.
#' @param fold_ids Run only these folds: station ids for `"lso"`, region ids
#'   for `"lro"`. Useful for quick checks on a subset.
#' @param seed Seed for posterior sampling (`NULL` for parallel, non-reproducible
#'   sampling; see [bruce_predict()]).
#' @param verbose Print progress.
#'
#' @return An object of class `bruce_cv` with per-observation `predictions`,
#'   per-fold `folds` (scores and fit status) and `overall` scores, all as
#'   tibbles.
#' @export
bruce_cv <- function(d, mesh, spec = bruce_spec(), window = NULL,
                     scheme = c("lso", "lro"), regions = NULL, region_col = NULL,
                     test_sources = NULL, n_samples = 200, level = 0.9, noise_sd = 0,
                     checkpoint_dir = NULL, theta_init = NULL, fix_theta = FALSE,
                     fold_ids = NULL, seed = 1, verbose = TRUE) {
  stopifnot(inherits(d, "bruce_data"))
  check_inla()
  scheme <- match.arg(scheme)
  days <- window_days(d, window)
  test_sources <- test_sources %||% d$sources$source[!d$sources$bias]
  st <- unique(d$obs[d$obs$date %in% days, c("station_id", "source", "x", "y")])

  folds <- if (scheme == "lso") {
    ids <- st$station_id[st$source %in% test_sources]
    lapply(ids, function(i) list(id = i, test = i, drop = i))
  } else {
    if (is.null(regions) || is.null(region_col)) {
      stop("`regions` and `region_col` are needed for scheme = 'lro'", call. = FALSE)
    }
    pts <- sf::st_as_sf(st, coords = c("x", "y"), crs = d$crs)
    reg <- sf::st_transform(sf::st_as_sf(regions)[, region_col], d$crs)
    st$region <- sf::st_join(pts, reg)[[region_col]]
    rids <- sort(unique(st$region[st$source %in% test_sources & !is.na(st$region)]))
    lapply(rids, function(r) list(
      id = as.character(r),
      test = st$station_id[st$region %in% r & st$source %in% test_sources],
      drop = st$station_id[st$region %in% r]))
  }
  if (!is.null(fold_ids)) {
    miss <- setdiff(as.character(fold_ids), vapply(folds, function(f) as.character(f$id), ""))
    if (length(miss)) warning("`fold_ids` not found among the folds: ", paste(miss, collapse = ", "), call. = FALSE)
    folds <- Filter(function(f) as.character(f$id) %in% as.character(fold_ids), folds)
  }
  if (!length(folds)) stop("No folds: no test stations in this window", call. = FALSE)
  if (fix_theta && is.null(theta_init)) {
    if (verbose) message("fix_theta: fitting the full data once to get the hyperparameter mode")
    theta_init <- bruce_fit(d, mesh, spec, window = range(days))$inla$mode$theta
  }
  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

  alpha <- (1 - level) / 2
  run_fold <- function(fd) {
    ck <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, sprintf("fold_%s_%s.rds", scheme, make.names(fd$id)))
    if (!is.null(ck) && file.exists(ck)) {
      if (verbose) message(sprintf("[%s] fold %s: loaded from checkpoint", scheme, fd$id))
      return(readRDS(ck))
    }
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      bruce_fit(d, mesh, spec, window = range(days), theta_init = theta_init,
                fix_theta = fix_theta, keep_config = n_samples > 0, drop_stations = fd$drop),
      error = function(e) e)
    if (inherits(fit, "error")) {
      if (verbose) message(sprintf("[%s] fold %s: fit failed (%s)", scheme, fd$id, conditionMessage(fit)))
      return(list(pred = NULL, info = data.frame(fold = fd$id, status = "error",
                  group_rho = NA, field_sd = NA, minutes = NA)))
    }
    te <- d$obs[d$obs$date %in% days & d$obs$station_id %in% fd$test, ]
    te$day_idx <- match(te$date, days)
    X <- matrix(1, nrow(te), 1, dimnames = list(NULL, "intercept"))
    if (length(fit$cov_names)) {
      cz <- d$covariates$station[match(te$station_id, d$covariates$station$station_id), fit$cov_names, drop = FALSE]
      X <- cbind(X, as.matrix(cz))
    }
    A <- INLA::inla.spde.make.A(mesh, loc = cbind(te$x, te$y))
    comps <- latent_components(fit, n_samples, seed)
    out <- vector("list", length(unique(te$day_idx)))
    for (k in seq_along(out)) {
      g <- sort(unique(te$day_idx))[k]
      r <- which(te$day_idx == g)
      P <- eval_predictor(fit, comps, X[r, , drop = FALSE], A[r, , drop = FALSE], g)
      if (noise_sd > 0 && ncol(P) > 1) P <- P + stats::rnorm(length(P), 0, noise_sd)
      res <- data.frame(fold = fd$id, station_id = te$station_id[r], date = te$date[r],
                        obs = te$value[r], pred = rowMeans(P))
      if (ncol(P) > 1) {
        res$pred_sd <- row_sd(P)
        qq <- row_quantiles(P, c(alpha, 1 - alpha))
        res$lower <- qq[, 1]
        res$upper <- qq[, 2]
        res$crps <- vapply(seq_along(r), function(i) crps_ensemble(te$value[r[i]], P[i, ]), 0)
      }
      out[[k]] <- res
    }
    pred <- do.call(rbind, out)
    last <- fit$attempts[nrow(fit$attempts), ]
    info <- data.frame(fold = fd$id, status = fit$status, group_rho = last$group_rho,
                       field_sd = last$field_sd,
                       minutes = (proc.time()[["elapsed"]] - t0) / 60)
    res <- list(pred = pred, info = info)
    if (!is.null(ck)) saveRDS(res, ck)
    if (verbose) message(sprintf("[%s] fold %s: %d obs, RMSE %.3f, status %s (%.1f min)",
                                 scheme, fd$id, nrow(pred),
                                 sqrt(mean((pred$pred - pred$obs)^2)), fit$status, info$minutes))
    res
  }

  res <- lapply(folds, run_fold)
  pred <- do.call(rbind, lapply(res, `[[`, "pred"))
  info <- do.call(rbind, lapply(res, `[[`, "info"))
  if (is.null(pred)) stop("All folds failed", call. = FALSE)

  score <- function(p) {
    e <- p$pred - p$obs
    out <- data.frame(n = nrow(p), rmse = sqrt(mean(e^2)), mae = mean(abs(e)), bias = mean(e))
    if ("crps" %in% names(p)) {
      out$crps <- mean(p$crps)
      out$coverage <- mean(p$obs >= p$lower & p$obs <= p$upper)
    }
    out
  }
  by_fold <- do.call(rbind, lapply(split(pred, pred$fold), function(p) cbind(fold = p$fold[1], score(p))))
  folds_df <- merge(info, by_fold, by = "fold", all.x = TRUE)

  structure(list(predictions = tibble::as_tibble(pred), folds = tibble::as_tibble(folds_df),
                 overall = tibble::as_tibble(score(pred)),
                 scheme = scheme, level = level, days = days),
            class = "bruce_cv")
}
