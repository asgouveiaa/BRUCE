#' Plot a prediction: one map per day
#'
#' One small map per day. With `scale = "common"` all days share one colour
#' scale, so days can be compared directly (from the 0.5% to the 99.5%
#' quantile of all days, so a few extreme cells do not wash out the rest).
#' With `scale = "free"` each day gets its own scale and legend, which shows
#' the pattern within the city on every day, even on days much cooler or
#' warmer than the rest. With `scale = "anomaly"` each day is shown minus
#' its own city mean, on one symmetric scale: how much hotter or cooler each
#' place is than the city that day, comparable across days.
#'
#' @param x A `bruce_prediction`.
#' @param what `"mean"` (default) or `"sd"` (uncertainty; needs a prediction
#'   made with `n_samples > 0`).
#' @param scale `"common"` (default), `"free"` or `"anomaly"` (see Details;
#'   `"anomaly"` only for `what = "mean"`).
#' @param dates Days to plot (default: all).
#' @param ncol Maps per row (default: up to 6).
#' @param main Overall title.
#' @param ... Unused.
#' @method plot bruce_prediction
#' @export
plot.bruce_prediction <- function(x, what = c("mean", "sd"), scale = c("common", "free", "anomaly"),
                                  dates = NULL, ncol = NULL, main = NULL, ...) {
  what <- match.arg(what)
  scale <- match.arg(scale)
  if (scale == "anomaly" && what == "sd") stop("`scale = \"anomaly\"` is for `what = \"mean\"`", call. = FALSE)
  r <- x[[what]]
  if (is.null(r)) stop("This prediction has no sd: predict with n_samples > 0", call. = FALSE)
  if (!is.null(dates)) {
    k <- match(as.character(as.Date(dates)), names(r))
    if (anyNA(k)) stop("Some `dates` are not in the prediction", call. = FALSE)
    r <- r[[k]]
  }
  n <- terra::nlyr(r)
  nc <- ncol %||% min(6, n)
  nl <- ceiling(n / nc)
  faixa <- function(v) {
    l <- stats::quantile(v, c(0.005, 0.995), na.rm = TRUE, names = FALSE)
    if (diff(l) <= 0) l <- l + c(-0.5, 0.5)
    l
  }
  if (scale == "anomaly") {
    r <- r - terra::global(r, "mean", na.rm = TRUE)[, 1]
    a <- max(abs(faixa(terra::values(r))))
    lim <- c(-a, a)
    cores <- grDevices::hcl.colors(50, "Blue-Red 3")
  } else {
    lim <- faixa(terra::values(r))
    cores <- if (what == "mean") grDevices::hcl.colors(50, "YlOrRd", rev = TRUE) else
      grDevices::hcl.colors(50, "Purples", rev = TRUE)
  }
  op <- graphics::par(mfrow = c(nl, nc), mar = c(0.5, 0.5, 1.6, 3.5),
                      oma = c(0, 0, if (!is.null(main)) 2 else 0, 0))
  on.exit(graphics::par(op))
  for (i in seq_len(n)) {
    if (scale == "free") {
      terra::plot(r[[i]], range = faixa(terra::values(r[[i]])), col = cores, axes = FALSE,
                  main = names(r)[i], legend = TRUE)
    } else {
      terra::plot(r[[i]], range = lim, col = cores, axes = FALSE, main = names(r)[i],
                  legend = (i %% nc == 0) || (i == n))
    }
  }
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, cex = 1.1)
  invisible(x)
}

#' Plot a cross-validation
#'
#' Left: each held-out reading against its prediction, with the predictive
#' interval as a vertical segment and the 1:1 line. Right: RMSE of each
#' held-out station (LSO) or region (LRO).
#'
#' @param x A `bruce_cv`.
#' @param ... Unused.
#' @method plot bruce_cv
#' @export
plot.bruce_cv <- function(x, ...) {
  p <- as.data.frame(x$predictions)
  if (!nrow(p)) stop("No predictions to plot", call. = FALSE)
  grp <- if (x$scheme == "lro" && "fold" %in% names(p)) as.character(p$fold) else p$station_id
  e <- p$pred - p$obs
  has_int <- all(c("lower", "upper") %in% names(p))
  lim <- range(c(p$obs, p$pred, if (has_int) c(p$lower, p$upper)), na.rm = TRUE)
  g <- factor(grp)
  cores <- grDevices::hcl.colors(max(3, nlevels(g)), "Dark 3")[as.integer(g)]

  op <- graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  on.exit(graphics::par(op))
  graphics::plot(p$obs, p$pred, pch = 19, cex = 0.7, col = cores, xlim = lim, ylim = lim,
                 xlab = "Held-out reading", ylab = "Prediction",
                 main = sprintf("%s: RMSE %.2f, bias %+.2f", toupper(x$scheme), sqrt(mean(e^2)), mean(e)))
  if (has_int) graphics::segments(p$obs, p$lower, p$obs, p$upper, col = grDevices::adjustcolor(cores, 0.25))
  graphics::abline(0, 1, lty = 2)
  rmse_g <- sort(tapply(e, g, function(v) sqrt(mean(v^2))))
  op2 <- graphics::par(mar = c(4.5, max(4, 0.45 * max(nchar(names(rmse_g)))), 3, 1))
  graphics::barplot(rmse_g, horiz = TRUE, las = 1, cex.names = 0.7, xlab = "RMSE",
                    main = if (x$scheme == "lso") "By station" else "By region")
  graphics::par(op2)
  invisible(x)
}

#' City series: the map against the reference stations
#'
#' For each day, the median of the mapped surface over the grid and its 5-95%
#' spread across cells, with the daily mean of the reference stations'
#' readings. A quick check that the map follows the stations through the
#' period (level and day-to-day changes).
#'
#' @param prediction A `bruce_prediction`.
#' @param data The `bruce_data` (or `bruce_fit`) with the observations.
#' @param sources Networks averaged as points (default: the reference
#'   networks).
#' @param ylab Axis label.
#' @return A data frame with the daily values, invisibly.
#' @export
bruce_plot_series <- function(prediction, data, sources = NULL, ylab = "Value") {
  stopifnot(inherits(prediction, "bruce_prediction"))
  if (inherits(data, "bruce_fit")) data <- data$data
  stopifnot(inherits(data, "bruce_data"))
  sources <- sources %||% data$sources$source[!data$sources$bias]
  r <- prediction$mean
  q <- t(vapply(seq_len(terra::nlyr(r)), function(i)
    stats::quantile(terra::values(r[[i]]), c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE), numeric(3)))
  dias <- as.Date(names(r))
  ob <- data$obs[data$obs$source %in% sources & data$obs$date %in% dias, ]
  est <- tapply(ob$value, as.character(ob$date), mean)
  out <- data.frame(date = dias, map_p05 = q[, 1], map_median = q[, 2], map_p95 = q[, 3],
                    stations_mean = unname(est[as.character(dias)]))
  graphics::plot(dias, q[, 2], type = "n", ylim = range(c(q, out$stations_mean), na.rm = TRUE),
                 xlab = "", ylab = ylab,
                 main = "Map (median and 5-95% over the grid) and reference stations")
  graphics::polygon(c(dias, rev(dias)), c(q[, 1], rev(q[, 3])),
                    col = grDevices::adjustcolor("orange", 0.3), border = NA)
  graphics::lines(dias, q[, 2], lwd = 2, col = "darkorange3")
  graphics::points(dias, out$stations_mean, pch = 19, col = "steelblue")
  graphics::legend("bottomleft", c("map: median", "map: 5-95%", "reference stations: mean"),
                   col = c("darkorange3", grDevices::adjustcolor("orange", 0.5), "steelblue"),
                   lwd = c(2, 8, NA), pch = c(NA, NA, 19), bty = "n")
  invisible(out)
}
