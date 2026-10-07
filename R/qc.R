#' Quality control of stations against the background
#'
#' Compares every station with the background product (the nearest
#' background cell, same day) and flags stations that behave unlike the rest
#' of their own network. A network's typical difference to the background is
#' expected (that is its bias, which the model estimates); what is suspect is
#' a station far from its own network, or one that swings much more than the
#' others from day to day. A broken humidity sensor, a station with wrong
#' coordinates or a sensor in direct sun show up this way.
#'
#' Without a background, each station is compared with the daily mean of the
#' reference stations instead.
#'
#' Stations of reference networks (`bias = FALSE`) are reported but never
#' dropped: a reference station far from the background is more likely real
#' local climate (a valley, the waterfront) than a fault.
#'
#' @param d A `bruce_data`.
#' @param threshold Flag a station whose mean difference is more than
#'   `threshold` (in the variable's units) away from its network's median
#'   difference, or whose day-to-day standard deviation of the difference is
#'   above `threshold`.
#' @param min_days Flag stations with fewer days of data than this.
#' @param drop `FALSE` (default) returns the report. `TRUE` returns `d` with
#'   the flagged stations of biased networks removed, and the report in
#'   `d$qc`.
#'
#' @return With `drop = FALSE`, a tibble (class `bruce_qc`) with one row per
#'   station: `n_days`, `n_compared` (days with a background value),
#'   `mean_diff` and `sd_diff` (station minus background), `deviation` (mean
#'   difference minus the network's median), `flag` (reason, empty if none)
#'   and `drop`. With `drop = TRUE`, the filtered `bruce_data`.
#' @examples
#' sim <- bruce_simulate(n_days = 5, seed = 1)
#' d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
#' bruce_qc(d)
#' @export
bruce_qc <- function(d, threshold = 3, min_days = 5, drop = FALSE) {
  stopifnot(inherits(d, "bruce_data"))
  if (!is.numeric(threshold) || threshold <= 0) stop("`threshold` must be positive", call. = FALSE)
  obs <- as.data.frame(d$obs)
  ref_src <- d$sources$source[!d$sources$bias]

  if (!is.null(d[["background"]])) {
    bg <- as.data.frame(d[["background"]])
    cel <- unique(bg[, c("x", "y")])
    st <- unique(obs[, c("station_id", "x", "y")])
    st$cel <- vapply(seq_len(nrow(st)), function(i)
      which.min((cel$x - st$x[i])^2 + (cel$y - st$y[i])^2), 1L)
    bg$cel <- match(paste(bg$x, bg$y), paste(cel$x, cel$y))
    obs$cel <- st$cel[match(obs$station_id, st$station_id)]
    obs$against <- bg$value[match(paste(obs$cel, obs$date), paste(bg$cel, bg$date))]
    against <- "background"
  } else {
    ref <- obs[obs$source %in% ref_src, ]
    if (!nrow(ref)) stop("No background and no reference network to compare with", call. = FALSE)
    m <- tapply(ref$value, as.character(ref$date), mean)
    obs$against <- unname(m[as.character(obs$date)])
    against <- "reference stations' daily mean"
  }
  obs$dif <- obs$value - obs$against

  rows <- lapply(split(obs, obs$station_id), function(z) {
    ok <- !is.na(z$dif)
    data.frame(station_id = z$station_id[1], source = z$source[1],
               n_days = length(unique(z$date)), n_compared = sum(ok),
               mean_diff = if (any(ok)) mean(z$dif[ok]) else NA_real_,
               sd_diff = if (sum(ok) > 1) stats::sd(z$dif[ok]) else NA_real_)
  })
  rep <- do.call(rbind, rows)
  med <- tapply(rep$mean_diff, rep$source, stats::median, na.rm = TRUE)
  rep$network_median <- unname(med[rep$source])
  rep$deviation <- rep$mean_diff - rep$network_median
  rep$role <- ifelse(rep$source %in% ref_src, "reference", "biased")

  why <- function(i) {
    r <- rep[i, ]
    w <- character(0)
    if (r$n_days < min_days) w <- c(w, sprintf("only %d day(s)", r$n_days))
    if (r$n_compared >= min_days) {
      if (abs(r$deviation) > threshold) w <- c(w, sprintf("%+.1f from its network", r$deviation))
      if (!is.na(r$sd_diff) && r$sd_diff > threshold) w <- c(w, sprintf("day-to-day sd %.1f", r$sd_diff))
    }
    paste(w, collapse = "; ")
  }
  rep$flag <- vapply(seq_len(nrow(rep)), why, "")
  rep$drop <- nzchar(rep$flag) & rep$role == "biased"
  rep <- rep[order(-nzchar(rep$flag), -abs(rep$deviation)), c(
    "station_id", "source", "role", "n_days", "n_compared", "mean_diff", "sd_diff",
    "network_median", "deviation", "flag", "drop")]
  rownames(rep) <- NULL
  rep <- tibble::as_tibble(rep)
  attr(rep, "against") <- against
  attr(rep, "threshold") <- threshold
  class(rep) <- c("bruce_qc", class(rep))

  if (!drop) return(rep)
  out <- rep$station_id[rep$drop]
  if (length(out)) {
    message(sprintf("bruce_qc: dropped %d station(s): %s", length(out), paste(out, collapse = ", ")))
    d$obs <- d$obs[!d$obs$station_id %in% out, ]
    if (!is.null(d$stations)) d$stations <- d$stations[!d$stations$station_id %in% out, ]
    if (!is.null(d$covariates$station)) {
      d$covariates$station <- d$covariates$station[!d$covariates$station$station_id %in% out, ]
    }
  } else {
    message("bruce_qc: no station dropped")
  }
  d$qc <- rep
  d
}

#' @method print bruce_qc
#' @export
print.bruce_qc <- function(x, ...) {
  cat(sprintf("<bruce_qc> %d stations compared with the %s (threshold %g)\n",
              nrow(x), attr(x, "against") %||% "background", attr(x, "threshold") %||% NA))
  df <- as.data.frame(x)
  net <- do.call(rbind, lapply(split(df, df$source), function(z)
    data.frame(source = z$source[1], role = z$role[1], stations = nrow(z),
               median_diff = round(z$network_median[1], 2),
               flagged = sum(nzchar(z$flag)), dropped = sum(z$drop))))
  cat("Station minus background, by network:\n")
  print(net, row.names = FALSE)
  fl <- df[nzchar(df$flag), ]
  if (nrow(fl)) {
    cat("\nFlagged stations:\n")
    fl$mean_diff <- round(fl$mean_diff, 2); fl$sd_diff <- round(fl$sd_diff, 2)
    fl$action <- ifelse(fl$drop, "drop", "keep (reference network: check by hand)")
    print(fl[, c("station_id", "source", "n_days", "mean_diff", "sd_diff", "flag", "action")], row.names = FALSE)
  } else {
    cat("\nNo station flagged.\n")
  }
  invisible(x)
}
