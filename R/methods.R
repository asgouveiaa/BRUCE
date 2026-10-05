#' Diagnostics of a fit
#'
#' Reports every fitting attempt (default, restart, fixed correlation), the
#' key hyperparameters and whether the field collapsed.
#'
#' @param fit A `bruce_fit`.
#' @return A data frame of attempts, invisibly; printed.
#' @export
bruce_check <- function(fit) {
  stopifnot(inherits(fit, "bruce_fit"))
  cat(sprintf("Status: %s\n", fit$status))
  if (!is.null(fit$group_rho_fixed)) {
    cat(sprintf("Day-to-day field correlation fixed at %.2f\n", fit$group_rho_fixed))
  }
  cat(sprintf("Range prior median: %.0f m (posterior mean range: %.0f m)\n", fit$range_prior,
              fit$attempts$range[nrow(fit$attempts)]))
  mr <- fit$mesh_resolution
  if (!is.null(mr) && is.finite(mr$ratio)) {
    cat(sprintf("Mesh resolution: typical edge near the data %.0f m; range / edge = %.1f (%s)\n",
                mr$edge, mr$ratio,
                if (mr$ok) "ok" else sprintf("below %g: range estimate may depend on the mesh", mr$min_ratio)))
  }
  if (!is.null(fit$minutes)) cat(sprintf("Fitting time: %.1f min\n", fit$minutes))
  print(fit$attempts, row.names = FALSE)
  invisible(fit$attempts)
}

#' @export
print.bruce_data <- function(x, ...) {
  cat("<bruce_data>\n")
  cat(sprintf("  %d observations, %d stations, %s to %s\n", nrow(x$obs),
              length(unique(x$obs$station_id)), min(x$obs$date), max(x$obs$date)))
  tab <- table(x$obs$source)
  src <- merge(x$sources, data.frame(source = names(tab), n_obs = as.integer(tab)), all.x = TRUE)
  src$role <- ifelse(src$bias, "biased", "reference")
  print(src[, c("source", "role", "obs_sd", "n_obs")], row.names = FALSE)
  if (!is.null(x$background)) {
    cat(sprintf("  background: %d cells x days (sd = %.2f)\n", nrow(x$background), x$background_sd))
  }
  if (!is.null(x$covariates)) cat("  covariates:", paste(x$covariates$names, collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.bruce_fit <- function(x, ...) {
  cat("<bruce_fit>", format(min(x$days)), "to", format(max(x$days)),
      sprintf("(%d days), status: %s\n", length(x$days), x$status))
  print(x$blocks, row.names = FALSE)
  fx <- x$inla$summary.fixed[, c("mean", "sd", "0.025quant", "0.975quant")]
  print(round(fx, 3))
  invisible(x)
}

#' @export
summary.bruce_fit <- function(object, ...) {
  print(object)
  hp <- object$inla$summary.hyperpar
  if (isTRUE(object$fixed_theta) || is.null(hp) || !nrow(hp)) {
    cat("\nHyperparameters fixed at `theta_init` (not estimated in this fit)\n")
  } else {
    cat("\nHyperparameters:\n")
    print(round(hp[, c("mean", "sd", "0.025quant", "0.975quant")], 3))
  }
  cat("\n")
  bruce_check(object)
  invisible(object)
}

#' @export
print.bruce_cv <- function(x, ...) {
  cat(sprintf("<bruce_cv> %s, %d folds, %d held-out observations\n",
              toupper(x$scheme), nrow(x$folds), nrow(x$predictions)))
  print(round(x$overall, 3), row.names = FALSE)
  bad <- x$folds$status[x$folds$status != "ok"]
  if (length(bad)) cat(sprintf("Folds not 'ok': %d (%s)\n", length(bad),
                               paste(unique(bad), collapse = ", ")))
  invisible(x)
}

#' @export
print.bruce_prediction <- function(x, ...) {
  cat(sprintf("<bruce_prediction> %d day(s), %d cells, %s\n", length(x$dates), length(x$cells),
              if (x$n_samples > 1) sprintf("%d posterior samples", x$n_samples) else "posterior means"))
  invisible(x)
}
