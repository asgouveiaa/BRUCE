#' Fit the model for one time window
#'
#' @param d A `bruce_data` object (usually after [bruce_covariates()]).
#' @param mesh Mesh from [bruce_mesh()].
#' @param spec Model specification from [bruce_spec()].
#' @param window Two dates (first and last day). `NULL` uses the full range
#'   of `d$obs`. Windows are fitted independently; one month is the tested
#'   scale.
#' @param theta_init Optional starting point for the hyperparameters: a
#'   `bruce_fit` (recommended: its structure is checked against this model
#'   and a clear error explains any mismatch) or a vector on INLA's internal
#'   scale (`other_fit$inla$mode$theta`). Without `fix_theta` it is only a
#'   starting point (warm start): the hyperparameters are still estimated.
#' @param fix_theta Keep the hyperparameters fixed at `theta_init` instead
#'   of estimating them: INLA then computes a single Gaussian approximation
#'   of the latent field, which is much faster. Meant for operational maps
#'   with hyperparameters learned from a previous fit, and for fast
#'   cross-validation. Hyperparameter uncertainty is not propagated, and
#'   `theta_init` must come from a model with the same structure (same
#'   networks with estimated noise, same optional terms).
#' @param keep_config Keep the latent configuration so that posterior samples
#'   can be drawn (needed by `bruce_predict(n_samples > 0)`).
#' @param cpo Compute CPO/PIT.
#' @param drop_stations Station ids excluded from fitting (used by
#'   [bruce_cv()]).
#' @param verbose Passed to INLA.
#'
#' @details The linear predictor for an observation of network `k` at station
#'   `i`, location `s` and day `t` is
#'   \deqn{\eta = \beta_0 + \delta_k + x(s)^\top\beta + u_i + d_t + w(s, t)}
#'   where \eqn{\delta_k} is the network offset (only for `bias = TRUE`),
#'   \eqn{u_i} the iid station effect, \eqn{d_t} the city-wide daily AR(1)
#'   term and \eqn{w(s,t)} a Matern SPDE field with AR(1) dependence between
#'   days. Each network, and the background, has its own Gaussian
#'   likelihood. Background cells share \eqn{\beta_0}, covariates, \eqn{d_t}
#'   and \eqn{w} but no offset or station effect.
#'
#'   After fitting, the field is checked for collapse (see [bruce_spec()]);
#'   every attempt is recorded in `fit$attempts` and the outcome in
#'   `fit$status`.
#'
#'   The mesh resolution is also checked against the estimated range: the
#'   typical edge of the triangles that contain the data is compared with the
#'   posterior mean range (`fit$mesh_resolution`). When the range is less than
#'   five times that edge, a warning is issued: the field is poorly resolved
#'   and the range estimate can change with the triangulation (in the Rio de
#'   Janeiro application, two meshes of about 1,000 vertices built with the
#'   same settings by different fmesher versions gave ranges of 7.7 and
#'   4.1 km).
#'
#' @return An object of class `bruce_fit`.
#' @export
bruce_fit <- function(d, mesh, spec = bruce_spec(), window = NULL, theta_init = NULL,
                      fix_theta = FALSE, keep_config = FALSE, cpo = FALSE,
                      drop_stations = NULL, verbose = FALSE) {
  stopifnot(inherits(d, "bruce_data"), inherits(spec, "bruce_spec"))
  check_inla()
  if (fix_theta && is.null(theta_init)) stop("`fix_theta = TRUE` needs `theta_init`", call. = FALSE)
  t0 <- proc.time()[["elapsed"]]
  b <- build_model(d, mesh, spec, window, drop_stations)
  y_sd <- stats::sd(b$obs$value)

  run <- function(theta = NULL) {
    tryCatch(
      run_inla(b, spec, theta = theta, fixed = fix_theta, keep_config = keep_config, cpo = cpo,
               verbose = verbose),
      error = function(e) {
        if (is.null(theta)) stop(e)
        stop("INLA failed when started from `theta_init` (", conditionMessage(e), "). ",
             "`theta_init` must come from a model with the same structure; pass the source ",
             "bruce_fit itself as `theta_init` to have the structure checked.", call. = FALSE)
      })
  }
  record <- function(fit, label) {
    chk <- check_collapse(fit, spec, y_sd)
    data.frame(attempt = length(attempts) + 1, strategy = label,
               group_rho = chk$group_rho, field_sd = chk$field_sd, range = chk$range,
               collapsed = chk$collapsed, stringsAsFactors = FALSE)
  }

  attempts <- list()
  fixed_rho <- spec$group_rho_fixed
  if (!is.null(fixed_rho)) set_group_rho(b, fixed_rho, spec)

  theta_source <- NULL
  if (inherits(theta_init, "bruce_fit")) {
    src <- theta_init
    theta_source <- sprintf("the fit of %s to %s", format(min(src$days)), format(max(src$days)))
    if (!is.null(src$theta_signature)) {
      dif <- compare_signature(src$theta_signature, theta_signature(b, spec, fixed_rho))
      if (length(dif)) {
        stop("`theta_init` comes from a model with a different structure:\n  - ",
             paste(dif, collapse = "\n  - "),
             "\nFit without `theta_init`, or borrow it from a fit with the same networks and options.",
             call. = FALSE)
      }
    }
    theta_init <- src$inla$mode$theta
  } else if (!is.null(theta_init)) {
    theta_source <- "user-supplied values"
  }

  fit <- run(theta_init)
  if (!is.null(theta_init) && length(fit$mode$theta) != length(theta_init)) {
    stop(sprintf(paste0("`theta_init` has %d values but this model has %d hyperparameters; ",
                        "it must come from a model with the same structure"),
                 length(theta_init), length(fit$mode$theta)), call. = FALSE)
  }
  attempts[[1]] <- record(fit, if (fix_theta) "hyperparameters fixed at theta_init"
                               else if (is.null(theta_init)) "default" else "user theta_init")
  status <- if (attempts[[1]]$collapsed) "collapsed" else "ok"

  if (status == "collapsed" && fix_theta) {
    warning("Spatial field collapsed at the fixed `theta_init`; results are unreliable", call. = FALSE)
  } else if (status == "collapsed" && is.null(fixed_rho) && spec$on_collapse != "none") {
    if (spec$on_collapse == "warn") {
      warning("Spatial field collapsed (see fit$attempts); results are unreliable", call. = FALSE)
    }
    if (spec$on_collapse == "retry_then_fix") {
      theta2 <- reset_theta(fit, y_sd)
      if (!is.null(theta2)) {
        fit2 <- tryCatch(run(theta2), error = function(e) NULL)
        if (!is.null(fit2)) {
          attempts[[length(attempts) + 1]] <- record(fit2, "restart from reset mode")
          if (!attempts[[length(attempts)]]$collapsed) {
            fit <- fit2
            status <- "ok_after_restart"
          }
        }
      }
    }
    if (status == "collapsed" && spec$on_collapse %in% c("retry_then_fix", "fix")) {
      fixed_rho <- spec$fallback_group_rho
      set_group_rho(b, fixed_rho, spec)
      fit <- run(NULL)
      attempts[[length(attempts) + 1]] <- record(
        fit, sprintf("group correlation fixed at %.2f", fixed_rho))
      status <- "fixed_group_rho"
      message(sprintf("Field collapsed; refitted with day-to-day correlation fixed at %.2f",
                      fixed_rho))
    }
  }

  obs_xy <- rbind(as.matrix(b$obs[, c("x", "y")]),
                  if (!is.null(d$background)) as.matrix(d$background[, c("x", "y")]))
  mres <- tryCatch(mesh_resolution(mesh, obs_xy, attempts[[length(attempts)]]$range),
                   error = function(e) list(edge = NA_real_, range = NA_real_, ratio = NA_real_,
                                            min_ratio = 5, ok = NA))
  if (isTRUE(!mres$ok)) {
    warning(sprintf(paste0(
      "Estimated range (%.0f m) is only %.1f times the typical mesh edge near the data (%.0f m). ",
      "The SPDE field is poorly resolved at this mesh, and the range estimate can change with ",
      "the triangulation; consider a finer mesh (max_edge and cutoff of about range / 5). ",
      "See bruce_check()."), mres$range, mres$ratio, mres$edge), call. = FALSE)
  }

  structure(
    list(inla = fit, data = d, mesh = mesh, spde = b$spde, spec = spec, mesh_resolution = mres,
         days = b$days, station_ids = b$station_ids, cov_names = b$cov_names,
         bias_cols = b$bias_cols, blocks = b$blocks, range_prior = b$range_prior,
         group_rho_fixed = fixed_rho, drop_stations = drop_stations,
         attempts = do.call(rbind, attempts), status = status, fixed_theta = fix_theta,
         theta_init_source = theta_source,
         theta_signature = theta_signature(b, spec, fixed_rho),
         group_rho_reason = if (!is.null(spec$group_rho_fixed)) "set in the specification"
                            else if (status == "fixed_group_rho") "the field collapsed when it was estimated",
         minutes = (proc.time()[["elapsed"]] - t0) / 60,
         inla_version = as.character(utils::packageVersion("INLA"))),
    class = "bruce_fit"
  )
}

# Builds stacks, formula and controls ------------------------------------------
build_model <- function(d, mesh, spec, window, drop_stations = NULL) {
  cov_names <- if (is.null(d$covariates)) character(0) else d$covariates$names
  days <- window_days(d, window)
  n_days <- length(days)

  obs <- d$obs[d$obs$date %in% days & !(d$obs$station_id %in% drop_stations), ]
  if (nrow(obs) == 0) stop("No observations in this window", call. = FALSE)
  obs$day_idx <- match(obs$date, days)
  if (length(cov_names)) {
    obs <- merge(obs, d$covariates$station, by = "station_id", all.x = TRUE, sort = FALSE)
  }
  station_ids <- sort(unique(obs$station_id))
  obs$station_idx <- match(obs$station_id, station_ids)

  bg <- NULL
  if (!is.null(d$background)) {
    bg <- d$background[d$background$date %in% days, ]
    bg$day_idx <- match(bg$date, days)
    if (length(cov_names)) bg <- merge(bg, d$covariates$background, by = c("x", "y"), sort = FALSE)
    if (nrow(bg) == 0) bg <- NULL
  }

  # ---- SPDE -----------------------------------------------------------------
  r0 <- if (identical(spec$range_prior, "nn")) {
    attr(mesh, "bruce_nn") %||% nn_median(obs[, c("x", "y")])
  } else spec$range_prior
  spde <- INLA::inla.spde2.pcmatern(mesh, alpha = 2, constr = TRUE,
                                    prior.range = c(r0, spec$range_prior_prob),
                                    prior.sigma = spec$sigma_prior)
  s_index <- INLA::inla.spde.make.index("s", n.spde = spde$n.spde, n.group = n_days)

  # ---- one likelihood block per network (+ background) ---------------------
  src_used <- d$sources[d$sources$source %in% unique(obs$source), ]
  bias_sources <- src_used$source[src_used$bias]
  bias_cols <- if (length(bias_sources)) make.unique(paste0("bias_", make.names(bias_sources))) else character(0)

  blocks <- list()
  if (!is.null(bg)) {
    blocks[[1]] <- list(tag = "background", source = NA_character_, df = bg, sd = d$background_sd)
  }
  for (k in seq_len(nrow(src_used))) {
    blocks[[length(blocks) + 1]] <- list(
      tag = src_used$source[k], source = src_used$source[k],
      df = obs[obs$source == src_used$source[k], ], sd = src_used$obs_sd[k])
  }
  n_lik <- length(blocks)

  make_effects <- function(df, src) {
    n <- nrow(df)
    e <- data.frame(intercept = rep(1, n))
    for (k in seq_along(bias_sources)) {
      e[[bias_cols[k]]] <- rep(as.numeric(!is.na(src) && src == bias_sources[k]), n)
    }
    for (cn in cov_names) e[[cn]] <- df[[cn]]
    e$day_idx <- df$day_idx
    e$station_idx <- if (is.na(src)) rep(NA_integer_, n) else df$station_idx
    e
  }

  stks <- vector("list", n_lik)
  control_family <- vector("list", n_lik)
  for (j in seq_len(n_lik)) {
    bl <- blocks[[j]]
    y <- matrix(NA_real_, nrow(bl$df), n_lik)
    y[, j] <- bl$df$value
    A <- INLA::inla.spde.make.A(mesh, loc = cbind(bl$df$x, bl$df$y),
                                group = bl$df$day_idx, n.group = n_days)
    stks[[j]] <- INLA::inla.stack(tag = bl$tag, data = list(y = y), A = list(1, A),
                                  effects = list(make_effects(bl$df, bl$source), s_index))
    control_family[[j]] <- if (!is.na(bl$sd)) {
      list(hyper = list(prec = list(initial = log(1 / bl$sd^2), fixed = TRUE)))
    } else {
      list(hyper = list(prec = list(prior = "pc.prec", param = spec$biased_prec_prior)))
    }
  }
  stk <- do.call(INLA::inla.stack, unname(stks))

  # ---- formula (objects referenced by name live in `env`) -------------------
  env <- new.env(parent = globalenv())
  env$f <- INLA::f
  env$spde_obj <- spde
  env$day_values <- seq_len(n_days)
  env$n_groups <- n_days      # all days of the window, even those without any data
  env$hyp_station <- list(prec = list(prior = "pc.prec", param = spec$station_prior))
  env$hyp_day <- list(
    rho = list(prior = "pc.cor1", param = spec$day_rho_prior),
    prec = if (is.null(spec$day_sd_fixed)) list(prior = "pc.prec", param = spec$day_prec_prior)
           else list(initial = log(1 / spec$day_sd_fixed^2), fixed = TRUE))
  env$ctrl_group <- list(model = "ar1", hyper = list(
    rho = list(prior = "pc.cor1", param = spec$group_rho_prior)))

  terms <- c("intercept", bias_cols, cov_names)
  if (spec$station_effect) terms <- c(terms, "f(station_idx, model = 'iid', hyper = hyp_station)")
  if (spec$day_effect) terms <- c(terms, "f(day_idx, model = 'ar1', values = day_values, hyper = hyp_day)")
  terms <- c(terms, "f(s, model = spde_obj, group = s.group, ngroup = n_groups, control.group = ctrl_group)")
  formula <- stats::as.formula(paste("y ~ 0 +", paste(terms, collapse = " + ")), env = env)

  ref_vals <- obs$value[obs$source %in% src_used$source[!src_used$bias]]
  int_mean <- if (is.na(spec$intercept_prior[1])) mean(if (length(ref_vals)) ref_vals else obs$value)
              else spec$intercept_prior[1]
  mean_l <- list(intercept = int_mean, default = 0)
  prec_l <- list(intercept = spec$intercept_prior[2], default = 0.001)
  for (cn in cov_names) { mean_l[[cn]] <- 0; prec_l[[cn]] <- spec$covariate_prec }
  for (bc in bias_cols) { mean_l[[bc]] <- 0; prec_l[[bc]] <- spec$bias_prec }

  list(stack = stk, formula = formula, env = env, n_lik = n_lik,
       control_family = control_family,
       control_fixed = list(mean = mean_l, prec = prec_l),
       spde = spde, days = days, station_ids = station_ids, cov_names = cov_names,
       bias_cols = stats::setNames(bias_cols, bias_sources),
       blocks = data.frame(likelihood = seq_len(n_lik),
                           tag = vapply(blocks, `[[`, "", "tag"),
                           n = vapply(blocks, function(bl) nrow(bl$df), 0L),
                           obs_sd = vapply(blocks, `[[`, 0, "sd")),
       range_prior = r0, obs = obs)
}

run_inla <- function(b, spec, theta = NULL, fixed = FALSE, keep_config = FALSE, cpo = FALSE,
                     verbose = FALSE) {
  e <- b$env
  data <- INLA::inla.stack.data(b$stack, spde_obj = e$spde_obj, hyp_station = e$hyp_station,
                                hyp_day = e$hyp_day, ctrl_group = e$ctrl_group,
                                day_values = e$day_values, n_groups = e$n_groups)
  args <- list(
    formula = b$formula, family = rep("gaussian", b$n_lik), data = data,
    control.predictor = list(A = INLA::inla.stack.A(b$stack), compute = FALSE),
    control.family = b$control_family, control.fixed = b$control_fixed,
    control.compute = list(config = keep_config, cpo = cpo),
    control.inla = c(list(int.strategy = spec$int_strategy),
                     if (!is.null(spec$strategy)) list(strategy = spec$strategy)),
    verbose = verbose
  )
  if (!is.null(theta)) {
    args$control.mode <- if (fixed) list(theta = theta, fixed = TRUE)
                         else list(theta = theta, restart = TRUE)
  }
  if (!is.null(spec$num_threads)) args$num.threads <- spec$num_threads
  do.call(INLA::inla, args)
}

set_group_rho <- function(b, rho, spec) {
  b$env$ctrl_group <- list(model = "ar1", hyper = list(
    rho = list(initial = rho_to_theta(rho), fixed = TRUE)))
  invisible(b)
}

hyper_row <- function(hp, pattern) {
  i <- grep(pattern, rownames(hp))
  if (length(i)) hp[i[1], "mean"] else NA_real_
}

check_collapse <- function(fit, spec, y_sd) {
  hp <- fit$summary.hyperpar
  rho <- hyper_row(hp, "^GroupRho for s")
  sdv <- hyper_row(hp, "^Stdev for s")
  rng <- hyper_row(hp, "^Range for s")
  collapsed <- isTRUE(rho > spec$collapse_rho) || isTRUE(sdv > spec$collapse_sd_ratio * y_sd)
  list(group_rho = rho, field_sd = sdv, range = rng, collapsed = collapsed)
}

# Moves the collapsed hyperparameters back to a plausible region so INLA can
# restart the mode search from there.
reset_theta <- function(fit, y_sd) {
  th <- fit$mode$theta
  nm <- names(th)
  if (is.null(nm)) return(NULL)
  k_rho <- which(grepl("group", nm, ignore.case = TRUE) & grepl("for s$", nm))
  k_sd <- which(grepl("stdev", nm, ignore.case = TRUE) & grepl("for s$", nm))
  if (length(k_rho) != 1) return(NULL)
  th[k_rho] <- rho_to_theta(0.8)
  if (length(k_sd) == 1) th[k_sd] <- log(y_sd)
  unname(th)
}

# Typical edge length of the mesh where the data are (median edge of the
# triangles that contain observation locations), compared with the estimated
# range. Rule of thumb for SPDE models: edges of about range / 5 or finer.
mesh_resolution <- function(mesh, xy, range, min_ratio = 5) {
  out <- list(edge = NA_real_, range = range, ratio = NA_real_, min_ratio = min_ratio, ok = NA)
  tv <- tryCatch(mesh$graph$tv, error = function(e) NULL)
  if (is.null(tv) || !is.finite(range)) return(out)
  loc <- mesh$loc[, 1:2, drop = FALSE]
  xy <- unique(xy)
  tri <- tryCatch({
    A <- fmesher::fm_basis(mesh, loc = xy)
    unique(apply(as.matrix(A) > 0, 1, function(r) {
      v <- which(r); hit <- which(apply(tv, 1, function(t) all(v %in% t)))
      if (length(hit)) hit[1] else NA_integer_
    }))
  }, error = function(e) NULL)
  tri <- tri[!is.na(tri)]
  if (!length(tri)) tri <- seq_len(nrow(tv))
  e <- unlist(lapply(tri, function(k) {
    p <- loc[tv[k, ], , drop = FALSE]
    c(sqrt(sum((p[1, ] - p[2, ])^2)), sqrt(sum((p[2, ] - p[3, ])^2)), sqrt(sum((p[3, ] - p[1, ])^2)))
  }))
  out$edge <- stats::median(e)
  out$ratio <- range / out$edge
  out$ok <- out$ratio >= min_ratio
  out
}

# What the hyperparameter vector of a model contains: used to check that a
# theta borrowed from another fit fits this one.
theta_signature <- function(b, spec, fixed_rho) {
  list(estimated_noise = b$blocks$tag[is.na(b$blocks$obs_sd)],
       station_effect = isTRUE(spec$station_effect),
       day_effect = isTRUE(spec$day_effect),
       day_sd_fixed = !is.null(spec$day_sd_fixed),
       group_rho_fixed = !is.null(fixed_rho))
}

compare_signature <- function(from, to) {
  fmt <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
  yn <- function(x) if (isTRUE(x)) "yes" else "no"
  out <- character(0)
  if (!identical(as.character(from$estimated_noise), as.character(to$estimated_noise))) {
    out <- c(out, sprintf("networks with estimated noise: %s in theta's model, %s here",
                          fmt(from$estimated_noise), fmt(to$estimated_noise)))
  }
  lab <- c(station_effect = "station effect", day_effect = "daily term",
           day_sd_fixed = "daily term sd fixed",
           group_rho_fixed = "day-to-day field correlation fixed")
  for (k in names(lab)) {
    if (!identical(isTRUE(from[[k]]), isTRUE(to[[k]]))) {
      out <- c(out, sprintf("%s: %s in theta's model, %s here", lab[[k]], yn(from[[k]]), yn(to[[k]])))
    }
  }
  out
}
