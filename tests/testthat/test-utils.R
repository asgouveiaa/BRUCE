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
