test_that("ensemble CRPS matches the closed form for a normal", {
  set.seed(2)
  mu <- 1; s <- 2; y <- 0.3
  z <- (y - mu) / s
  exact <- s * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi))
  approx <- bruce:::crps_ensemble(y, rnorm(2e5, mu, s))
  expect_equal(approx, exact, tolerance = 0.01)
})

test_that("AR(1) correlation maps to INLA's internal scale", {
  expect_equal(bruce:::rho_to_theta(0.55), 2 * atanh(0.55))
})

test_that("Matern(1) covariance is 1 at zero and decays", {
  v <- bruce:::matern1(c(0, 1000, 10000), kappa = sqrt(8) / 5000)
  expect_equal(v[1], 1)
  expect_true(all(diff(v) < 0))
  expect_equal(bruce:::matern1(5000, sqrt(8) / 5000), 0.139, tolerance = 0.01)
})

test_that("spec validates its inputs", {
  expect_s3_class(bruce_spec(), "bruce_spec")
  expect_error(bruce_spec(range_prior = "auto"), "range_prior")
  expect_error(bruce_spec(group_rho_fixed = 1), "group_rho_fixed")
})

test_that("mesh resolution compares range with edge length", {
  # two right triangles of side 1000 m; median edge = 1000
  m <- list(loc = cbind(c(0, 1000, 0, 1000), c(0, 0, 1000, 1000), 0),
            graph = list(tv = rbind(c(1, 2, 3), c(2, 4, 3))))
  r <- bruce:::mesh_resolution(m, cbind(500, 200), range = 3000)
  expect_equal(r$edge, 1000)
  expect_equal(r$ratio, 3)
  expect_false(r$ok)
  expect_true(bruce:::mesh_resolution(m, cbind(500, 200), range = 8000)$ok)
})

test_that("row_sd and row_quantiles match base R", {
  set.seed(1)
  P <- matrix(rnorm(300), 30)
  expect_equal(bruce:::row_sd(P), apply(P, 1, sd))
  q <- bruce:::row_quantiles(P, c(0.05, 0.5, 0.95))
  expect_equal(dim(q), c(30, 3))
  expect_equal(q, t(apply(P, 1, quantile, probs = c(0.05, 0.5, 0.95), names = FALSE)))
  expect_equal(bruce:::row_quantiles(P, 0.1)[, 1], apply(P, 1, quantile, probs = 0.1, names = FALSE))
})

test_that("bruce_grid builds a masked grid at the requested resolution", {
  sim <- bruce_simulate(n_days = 2, seed = 1)
  g1 <- bruce_grid(sim$domain, 1000)
  g2 <- bruce_grid(sim$domain, 500)
  expect_s4_class(g1, "SpatRaster")
  expect_equal(terra::res(g1), c(1000, 1000))
  n1 <- sum(!is.na(terra::values(g1))); n2 <- sum(!is.na(terra::values(g2)))
  expect_gt(n1, 0)
  expect_equal(n2 / n1, 4, tolerance = 0.1)
  ll <- sf::st_transform(sf::st_as_sf(sim$domain), 4326)
  expect_error(bruce_grid(ll, 1000), "longitude")
  expect_equal(terra::res(bruce_grid(ll, 1000, crs = sim$crs)), c(1000, 1000))
})

test_that("spec accepts the cheaper INLA strategies", {
  s <- bruce_spec(int_strategy = "eb", strategy = "gaussian")
  expect_equal(s$strategy, "gaussian")
  expect_null(bruce_spec()$strategy)
})

test_that("theta signatures are compared with readable differences", {
  a <- list(estimated_noise = c("aurassure", "wunder"), station_effect = TRUE,
            day_effect = TRUE, day_sd_fixed = FALSE, group_rho_fixed = FALSE)
  expect_length(bruce:::compare_signature(a, a), 0)
  b <- a; b$estimated_noise <- "wunder"; b$group_rho_fixed <- TRUE
  dif <- bruce:::compare_signature(a, b)
  expect_length(dif, 2)
  expect_match(dif[1], "aurassure, wunder in theta's model, wunder here")
  expect_match(dif[2], "correlation fixed: no in theta's model, yes here")
})

test_that("observation-scale offsets use network bias, station sd and noise", {
  fake <- list(
    inla = list(summary.fixed = data.frame(mean = 0.5, row.names = "bias_airport"),
                summary.hyperpar = data.frame(mean = c(4, 2), row.names = c(
                  "Precision for station_idx", "Precision for the Gaussian observations[3]")),
                mode = list(theta = c(`Log precision for station_idx` = log(4)))),
    bias_cols = c(airport = "bias_airport"), spec = list(station_effect = TRUE),
    blocks = data.frame(likelihood = 1:4, tag = c("background", "official", "wunder", "airport"),
                        n = 1, obs_sd = c(0.7, 0.3, NA, 0.3)))
  expect_equal(bruce:::station_sd(fake), 0.5)
  expect_equal(bruce:::noise_sd_source(fake, "official"), 0.3)
  expect_equal(bruce:::noise_sd_source(fake, "wunder"), 1 / sqrt(2))
  set.seed(1)
  off <- bruce:::reading_offsets(fake, c("official", "airport"), 20000)
  expect_equal(mean(off[1, ]), 0, tolerance = 0.02)
  expect_equal(mean(off[2, ]), 0.5, tolerance = 0.02)
  expect_equal(sd(off[1, ]), sqrt(0.5^2 + 0.3^2), tolerance = 0.02)
  # fixed hyperparameters: no posterior summary, read from the mode
  fake$inla$summary.hyperpar <- NULL
  expect_equal(bruce:::station_sd(fake), 0.5)
})
