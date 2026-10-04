sim <- bruce_simulate(n_reference = 8, n_biased = 15, n_days = 4, seed = 1)

test_that("simulated data passes the data contract", {
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
  expect_s3_class(d, "bruce_data")
  expect_equal(nrow(d$obs), (8 + 15) * 4)
  expect_true(all(c("x", "y", "source") %in% names(d$obs)))
})

test_that("contract violations are caught", {
  st_bad <- sim$stations
  st_bad$source[1] <- "unknown"
  expect_error(bruce_data(sim$obs, st_bad, sim$sources, crs = sim$crs), "not declared")

  obs_dup <- rbind(sim$obs, sim$obs[1, ])
  expect_error(bruce_data(obs_dup, sim$stations, sim$sources, crs = sim$crs), "more than one value")

  src_all_biased <- sim$sources
  src_all_biased$bias <- TRUE
  expect_error(bruce_data(sim$obs, sim$stations, src_all_biased, crs = sim$crs), "not identifiable")

  expect_error(bruce_data(sim$obs, sim$stations, sim$sources, crs = 4326), "projected")
})

test_that("lon/lat input is projected", {
  pts <- sf::st_as_sf(sim$stations, coords = c("x", "y"), crs = sim$crs)
  ll <- sf::st_coordinates(sf::st_transform(pts, 4326))
  st_ll <- data.frame(station_id = sim$stations$station_id, source = sim$stations$source,
                      lon = ll[, 1], lat = ll[, 2])
  d <- bruce_data(sim$obs, st_ll, sim$sources, crs = sim$crs)
  m <- match(d$stations$station_id, sim$stations$station_id)
  expect_lt(max(abs(d$stations$x - sim$stations$x[m])), 0.01)
})

test_that("covariates are extracted, standardised and stored", {
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
  d <- bruce_covariates(d, sim$covariates, background_buffer = 4500)
  sc <- d$covariates$scaling
  expect_equal(sc$name, c("elev", "urb"))
  expect_true(all(sc$scale_on == "background"))
  expect_equal(mean(d$covariates$background$elev), 0, tolerance = 1e-8)
  expect_false(anyNA(d$covariates$station$elev))

  d2 <- bruce_covariates(d, sim$covariates, background_method = c(urb = "zero"))
  expect_true(all(d2$covariates$background$urb == 0))
  expect_equal(d2$covariates$scaling$scale_on, c("background", "stations"))
})

test_that("reserved covariate names are rejected", {
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs)
  r <- sim$covariates
  names(r) <- c("intercept", "urb")
  expect_error(bruce_covariates(d, r), "cannot be")
})
