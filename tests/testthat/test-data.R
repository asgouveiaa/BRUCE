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

test_that("a list of rasters is read in each raster's own projection", {
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
  d1 <- bruce_covariates(d, sim$covariates, background_buffer = 4500)
  # same urb layer, but stored in lon/lat: values at stations must match closely
  urb_ll <- terra::project(sim$covariates$urb, "EPSG:4326", method = "near")
  d2 <- bruce_covariates(d, list(elev = sim$covariates$elev, urb = urb_ll),
                         background_buffer = 4500)
  expect_equal(d2$covariates$names, c("elev", "urb"))
  expect_true(is.list(d2$covariates$rasters))
  expect_equal(d1$covariates$station$elev, d2$covariates$station$elev)
  expect_gt(cor(d1$covariates$station$urb, d2$covariates$station$urb), 0.99)
})

test_that("na_value fills extraction gaps before standardisation", {
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
  half <- terra::crop(sim$covariates$urb, terra::ext(650000, 665000, 7450000, 7480000))
  st <- unique(d$obs[, c("station_id", "x")])
  outside <- st$station_id[st$x > 665000]
  dn <- bruce_covariates(d, list(urb = half), na_value = c(urb = 0), background_method = "point")
  sc <- dn$covariates$scaling
  z_out <- dn$covariates$station$urb[dn$covariates$station$station_id %in% outside]
  expect_true(all(abs(z_out - (0 - sc$mean) / sc$sd) < 1e-10))
  expect_equal(sc$na_value, 0)
  expect_error(bruce_covariates(d, list(urb = half), na_value = c(svf = 0)), "Unknown covariate")
  expect_error(bruce_covariates(d, list(half)), "named")
})
