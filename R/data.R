#' Assemble and validate input data
#'
#' The data contract of the package. Everything specific to a city enters here:
#' which networks exist, which ones are trusted, and where the stations are.
#'
#' @param obs Data frame of station observations with columns `station_id`,
#'   `date` (coercible to `Date`) and `value` (e.g. daily median or maximum
#'   temperature).
#' @param stations Data frame with one row per station: `station_id`, `source`,
#'   and either `x`/`y` in `crs` or `lon`/`lat` in WGS84.
#' @param sources Data frame describing each network, with columns
#'   * `source`: network name, matching `stations$source`;
#'   * `bias`: `FALSE` for reference networks (e.g. official stations),
#'     `TRUE` for networks that get an estimated offset (low-cost sensors,
#'     citizen stations, airports with a known systematic difference);
#'   * `obs_sd`: fixed observation standard deviation, or `NA` to estimate it
#'     with a PC prior.
#'   At least one reference network (`bias = FALSE`) is needed to identify the
#'   intercept.
#' @param crs Projected coordinate reference system in metres (EPSG code or
#'   anything accepted by [sf::st_crs()]).
#' @param background Optional gridded product on the same variable (e.g.
#'   ERA5-Land daily statistic), as a data frame with `date`, `value` and
#'   `x`/`y` or `lon`/`lat` of cell centres. Enters as an additional unbiased
#'   likelihood with fixed standard deviation `background_sd`.
#' @param background_sd Fixed observation standard deviation of the background.
#'
#' @return An object of class `bruce_data`.
#' @export
bruce_data <- function(obs, stations, sources, crs, background = NULL,
                       background_sd = 0.7) {
  obs <- as.data.frame(obs)
  stations <- as.data.frame(stations)
  sources <- as.data.frame(sources)
  need_cols(obs, c("station_id", "date", "value"), "obs")
  need_cols(stations, c("station_id", "source"), "stations")
  need_cols(sources, c("source", "bias", "obs_sd"), "sources")

  crs <- sf::st_crs(crs)
  if (is.na(crs)) stop("`crs` is not a valid coordinate reference system", call. = FALSE)
  if (isTRUE(sf::st_is_longlat(crs))) {
    stop("`crs` must be a projected system in metres (e.g. a UTM zone)", call. = FALSE)
  }

  # ---- sources --------------------------------------------------------------
  if (anyDuplicated(sources$source)) stop("`sources$source` has duplicates", call. = FALSE)
  if (!is.logical(sources$bias) || anyNA(sources$bias)) {
    stop("`sources$bias` must be TRUE/FALSE", call. = FALSE)
  }
  sources$obs_sd <- as.numeric(sources$obs_sd)
  if (any(sources$obs_sd <= 0, na.rm = TRUE)) stop("`sources$obs_sd` must be > 0 or NA", call. = FALSE)

  # ---- stations -------------------------------------------------------------
  stations$station_id <- as.character(stations$station_id)
  if (anyDuplicated(stations$station_id)) stop("`stations$station_id` has duplicates", call. = FALSE)
  bad_src <- setdiff(unique(stations$source), sources$source)
  if (length(bad_src)) {
    stop("stations use sources not declared in `sources`: ", paste(bad_src, collapse = ", "),
         call. = FALSE)
  }
  stations <- project_xy(stations, crs, "stations")

  # ---- observations ---------------------------------------------------------
  obs$station_id <- as.character(obs$station_id)
  obs$date <- as.Date(obs$date)
  if (anyNA(obs$date)) stop("`obs$date` has values that are not dates", call. = FALSE)
  n_na <- sum(is.na(obs$value))
  if (n_na) {
    message(sprintf("Dropping %d observation(s) with missing value", n_na))
    obs <- obs[!is.na(obs$value), ]
  }
  unknown <- setdiff(unique(obs$station_id), stations$station_id)
  if (length(unknown)) {
    stop(sprintf("%d station_id(s) in `obs` are not in `stations` (e.g. %s)",
                 length(unknown), paste(utils::head(unknown, 5), collapse = ", ")), call. = FALSE)
  }
  if (anyDuplicated(obs[, c("station_id", "date")])) {
    stop("`obs` has more than one value per station and date; aggregate to daily first",
         call. = FALSE)
  }
  obs <- merge(obs[, c("station_id", "date", "value")],
               stations[, c("station_id", "source", "x", "y")],
               by = "station_id", sort = FALSE)
  obs <- obs[order(obs$date, obs$station_id), ]
  rownames(obs) <- NULL

  used <- sources[sources$source %in% unique(obs$source), ]
  if (!any(!used$bias)) {
    msg <- "No reference network (bias = FALSE) has observations"
    if (is.null(background)) stop(msg, "; the intercept is not identifiable", call. = FALSE)
    warning(msg, "; the intercept rests on the background product alone", call. = FALSE)
  }

  # ---- background -----------------------------------------------------------
  if (!is.null(background)) {
    background <- as.data.frame(background)
    need_cols(background, c("date", "value"), "background")
    background$date <- as.Date(background$date)
    background <- project_xy(background, crs, "background")
    background <- background[!is.na(background$value), c("x", "y", "date", "value")]
  }

  structure(
    list(obs = obs, stations = stations, sources = sources,
         background = background, background_sd = background_sd,
         crs = crs, covariates = NULL),
    class = "bruce_data"
  )
}
