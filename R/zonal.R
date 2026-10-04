#' Aggregate predictions over zones
#'
#' Weighted zonal means (e.g. population-weighted exposure by neighbourhood
#' or health territory), propagating posterior uncertainty when samples were
#' kept.
#'
#' @param pred A `bruce_prediction`.
#' @param zones Polygons (`sf`).
#' @param id_col Name of the zone id column.
#' @param weights Optional [terra::SpatRaster] of weights (e.g. population).
#'   Resampled to the prediction grid by averaging, so only relative
#'   weights within a zone matter.
#' @param probs Quantiles of the zonal mean (only with samples).
#'
#' @return A data frame with one row per zone and date: `mean` and, with
#'   samples, `sd` and the requested quantiles of the zonal mean.
#' @export
bruce_zonal <- function(pred, zones, id_col, weights = NULL, probs = c(0.05, 0.95)) {
  stopifnot(inherits(pred, "bruce_prediction"))
  tmpl <- pred$template
  zones <- sf::st_transform(sf::st_as_sf(zones), terra::crs(tmpl))
  zones$.zid <- seq_len(nrow(zones))
  zr <- terra::rasterize(terra::vect(zones), tmpl, field = ".zid")
  zid <- terra::values(zr, mat = FALSE)[pred$cells]

  w <- rep(1, length(pred$cells))
  if (!is.null(weights)) {
    wr <- weights[[1]]
    if (!terra::same.crs(wr, tmpl)) wr <- terra::project(wr, tmpl, method = "average")
    if (!terra::compareGeom(wr, tmpl, stopOnError = FALSE)) wr <- terra::resample(wr, tmpl, method = "average")
    w <- terra::values(wr, mat = FALSE)[pred$cells]
    w[is.na(w)] <- 0
  }
  keep <- !is.na(zid) & w > 0
  zid <- zid[keep]; w <- w[keep]
  wsum <- rowsum(w, zid)

  out <- lapply(seq_along(pred$dates), function(k) {
    if (!is.null(pred$samples)) {
      M <- pred$samples[[k]][keep, , drop = FALSE]
      zm <- rowsum(M * w, zid) / as.vector(wsum)
      res <- data.frame(zid = as.integer(rownames(zm)), date = pred$dates[k],
                        mean = rowMeans(zm), sd = apply(zm, 1, stats::sd))
      qs <- t(apply(zm, 1, stats::quantile, probs = probs, names = FALSE))
      if (length(probs) == 1) qs <- t(qs)
      for (j in seq_along(probs)) res[[paste0("q", probs[j])]] <- qs[, j]
    } else {
      v <- terra::values(pred$mean[[k]], mat = FALSE)[pred$cells][keep]
      zm <- rowsum(v * w, zid) / as.vector(wsum)
      res <- data.frame(zid = as.integer(rownames(zm)), date = pred$dates[k], mean = zm[, 1])
    }
    res
  })
  out <- do.call(rbind, out)
  out[[id_col]] <- sf::st_drop_geometry(zones)[[id_col]][out$zid]
  out$zid <- NULL
  out[, c(id_col, setdiff(names(out), id_col))]
}
