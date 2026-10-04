#' Fit the model for one time window
#'
#' @param d A `bruce_data` object (usually after [bruce_covariates()]).
#' @param mesh Mesh from [bruce_mesh()].
#' @param spec Model specification from [bruce_spec()].
#' @param window Two dates (first and last day). `NULL` uses the full range
#'   of `d$obs`. Windows are fitted independently; one month is the tested
#'   scale.
#' @param theta_init Optional starting mode for the hyperparameters (on
#'   INLA's internal scale, e.g. `other_fit$inla$mode$theta` from a stable
#'   fit with the same structure).
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
#' @return An object of class `bruce_fit`.
#' @export
bruce_fit <- function(d, mesh, spec = bruce_spec(), window = NULL, theta_init = NULL,
                      keep_config = FALSE, cpo = FALSE, drop_stations = NULL,
                      verbose = FALSE) {
  stopifnot(inherits(d, "bruce_data"), inherits(spec, "bruce_spec"))
  check_inla()
  b <- build_model(d, mesh, spec, window, drop_stations)
  y_sd <- stats::sd(b$obs$value)

  run <- function(theta = NULL) {
    run_inla(b, spec, theta = theta, keep_config = keep_config, cpo = cpo, verbose = verbose)
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

  fit <- run(theta_init)
  attempts[[1]] <- record(fit, if (is.null(theta_init)) "default" else "user theta_init")
  status <- if (attempts[[1]]$collapsed) "collapsed" else "ok"

  if (status == "collapsed" && is.null(fixed_rho) && spec$on_collapse != "none") {
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

  structure(
    list(inla = fit, data = d, mesh = mesh, spde = b$spde, spec = spec,
         days = b$days, station_ids = b$station_ids, cov_names = b$cov_names,
         bias_cols = b$bias_cols, blocks = b$blocks, range_prior = b$range_prior,
         group_rho_fixed = fixed_rho, drop_stations = drop_stations,
         attempts = do.call(rbind, attempts), status = status,
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
  terms <- c(terms, "f(s, model = spde_obj, group = s.group, control.group = ctrl_group)")
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

run_inla <- function(b, spec, theta = NULL, keep_config = FALSE, cpo = FALSE, verbose = FALSE) {
  e <- b$env
  data <- INLA::inla.stack.data(b$stack, spde_obj = e$spde_obj, hyp_station = e$hyp_station,
                                hyp_day = e$hyp_day, ctrl_group = e$ctrl_group,
                                day_values = e$day_values)
  args <- list(
    formula = b$formula, family = rep("gaussian", b$n_lik), data = data,
    control.predictor = list(A = INLA::inla.stack.A(b$stack), compute = FALSE),
    control.family = b$control_family, control.fixed = b$control_fixed,
    control.compute = list(config = keep_config, cpo = cpo),
    control.inla = list(int.strategy = spec$int_strategy),
    verbose = verbose
  )
  if (!is.null(theta)) args$control.mode <- list(theta = theta, restart = TRUE)
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
