skip_if_not_installed("INLA")
skip_if_not_installed("fmesher")
skip_on_cran()

sim <- bruce_simulate(n_reference = 12, n_biased = 30, n_days = 5, seed = 3)
d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
d <- bruce_covariates(d, sim$covariates, background_buffer = 4500)
mesh <- bruce_mesh(d)

test_that("fit runs and recovers the main effects", {
  fit <- bruce_fit(d, mesh, keep_config = TRUE)
  expect_s3_class(fit, "bruce_fit")
  fx <- fit$inla$summary.fixed
  ci <- function(nm) fx[nm, c("0.025quant", "0.975quant")]
  # elevation effect on the standardised scale
  b_elev <- sim$params$beta_elev * d$covariates$scaling$sd[d$covariates$scaling$name == "elev"]
  expect_true(ci("elev")[[1]] < b_elev && b_elev < ci("elev")[[2]])
  expect_true(ci("bias_lowcost")[[1]] < sim$params$bias && sim$params$bias < ci("bias_lowcost")[[2]])
  expect_true(fit$status %in% c("ok", "ok_after_restart", "fixed_group_rho"))
})

test_that("prediction returns rasters close to the truth", {
  fit <- bruce_fit(d, mesh, keep_config = TRUE)
  p <- bruce_predict(fit, n_samples = 50, keep_samples = TRUE)
  expect_equal(terra::nlyr(p$mean), 5)
  err <- terra::values(p$mean[[3]] - terra::resample(sim$truth[[3]], p$template))
  expect_lt(sqrt(mean(err^2, na.rm = TRUE)), 1.5)
  z <- bruce_zonal(p, sim$domain |> transform(id = "city"), "id")
  expect_equal(nrow(z), 5)
})

test_that("leave-station-out runs and checkpoints", {
  ck <- tempfile()
  cv <- bruce_cv(d, mesh, scheme = "lso", n_samples = 50, checkpoint_dir = ck, verbose = FALSE)
  expect_s3_class(cv, "bruce_cv")
  expect_equal(length(list.files(ck)), 12)
  expect_true(all(c("rmse", "crps", "coverage") %in% names(cv$overall)))
})

test_that("decomposition adds up to the prediction", {
  fit <- bruce_fit(d, mesh)
  dec <- bruce_decompose(fit, dates = fit$days[2])
  p <- bruce_predict(fit, dates = fit$days[2])
  diff <- terra::values(dec$total[[1]] - p$mean[[1]])
  expect_lt(max(abs(diff), na.rm = TRUE), 1e-8)
  expect_equal(dec$variance$share_covariates + dec$variance$share_field, 1)
  f <- tempfile(fileext = ".png"); grDevices::png(f, 1200, 400); plot(dec); grDevices::dev.off()
  expect_true(file.exists(f))
})

test_that("fixed hyperparameters reproduce the full fit and speed up CV", {
  fit <- bruce_fit(d, mesh)
  th <- fit$inla$mode$theta
  fx <- bruce_fit(d, mesh, theta_init = th, fix_theta = TRUE)
  expect_true(fx$fixed_theta)
  expect_equal(fx$inla$summary.fixed$mean, fit$inla$summary.fixed$mean, tolerance = 0.05)
  expect_error(bruce_fit(d, mesh, fix_theta = TRUE), "theta_init")
  ids <- unique(d$obs$station_id[d$obs$source == "reference"])[1:2]
  cv <- bruce_cv(d, mesh, n_samples = 20, theta_init = th, fix_theta = TRUE,
                 fold_ids = ids, verbose = FALSE)
  expect_equal(nrow(cv$folds), 2)
})

test_that("a window can end on a day without data (forecast of the next day)", {
  fit <- bruce_fit(d, mesh)
  fim <- max(d$obs$date) + 1
  fx <- bruce_fit(d, mesh, window = c(min(d$obs$date), fim),
                  theta_init = fit$inla$mode$theta, fix_theta = TRUE)
  expect_equal(length(fx$days), length(fit$days) + 1)
  expect_equal(nrow(fx$inla$summary.random$s), fx$spde$n.spde * length(fx$days))
  p <- bruce_predict(fx, dates = fim)
  expect_true(all(is.finite(terra::values(p$mean[[1]], na.rm = TRUE))))
})

test_that("bruce_cv takes a fit and scores the station-reading scale", {
  fit <- bruce_fit(d, mesh, keep_config = TRUE)
  ids <- unique(d$obs$station_id[d$obs$source == "reference"])[1:3]
  cv <- bruce_cv(fit, n_samples = 50, fold_ids = ids, verbose = FALSE)
  expect_equal(cv$scale, "observation")
  expect_true(all(c("coverage", "coverage_field", "crps", "crps_field") %in% names(cv$overall)))
  expect_gte(cv$overall$coverage, cv$overall$coverage_field)   # wider intervals
  cvf <- bruce_cv(fit, n_samples = 50, fold_ids = ids, scale = "field", verbose = FALSE)
  expect_false("coverage_field" %in% names(cvf$overall))
  expect_output(print(cv), "station-reading")
})

test_that("borrowed theta is checked for structure and reported", {
  fit <- bruce_fit(d, mesh)
  fx <- bruce_fit(d, mesh, theta_init = fit, fix_theta = TRUE)
  expect_match(fx$theta_init_source, "the fit of")
  expect_output(bruce_check(fx), "FIXED at the fit of")
  expect_error(bruce_fit(d, mesh, bruce_spec(group_rho_fixed = 0.5), theta_init = fit, fix_theta = TRUE),
               "different structure")
  expect_error(bruce_fit(d, mesh, theta_init = fit$inla$mode$theta[-1], fix_theta = TRUE))
})

test_that("bruce_update maps a new day with the full fit's hyperparameters", {
  fit <- bruce_fit(d, mesh, window = range(d$obs$date)[1] + c(0, 2))
  u <- bruce_update(fit, d, date = max(d$obs$date), days_before = 2, n_samples = 0)
  expect_s3_class(u, "bruce_update")
  expect_true(isTRUE(u$fit$fixed_theta))
  expect_equal(terra::nlyr(u$prediction$mean), 1)
  expect_equal(names(u$prediction$mean), as.character(max(d$obs$date)))
  expect_output(print(u), "bruce_update")
  # a day after the data: one-day-ahead forecast
  u2 <- bruce_update(fit, d, date = max(d$obs$date) + 1, days_before = 2)
  expect_output(print(u2), "forecast")
  expect_warning(bruce_update(fit, d, date = max(d$obs$date), days_before = 2, max_age = 0), "Refit")
  d0 <- d; d0$covariates <- NULL
  expect_error(bruce_update(fit, d0), "Covariates differ")
})
