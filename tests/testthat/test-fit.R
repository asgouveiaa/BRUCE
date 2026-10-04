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
