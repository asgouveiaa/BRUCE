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
