test_that("bruce_qc flags a broken station of a biased network and keeps references", {
  sim <- bruce_simulate(n_days = 8, seed = 3)
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
  biased <- d$sources$source[d$sources$bias][1]
  ref <- d$sources$source[!d$sources$bias][1]
  bad <- d$obs$station_id[d$obs$source == biased][1]
  refst <- d$obs$station_id[d$obs$source == ref][1]
  d$obs$value[d$obs$station_id == bad] <- d$obs$value[d$obs$station_id == bad] - 8
  d$obs$value[d$obs$station_id == refst] <- d$obs$value[d$obs$station_id == refst] + 8

  q <- bruce_qc(d)
  expect_s3_class(q, "bruce_qc")
  expect_true(q$drop[q$station_id == bad])
  expect_match(q$flag[q$station_id == bad], "from its network")
  expect_true(nzchar(q$flag[q$station_id == refst]))
  expect_false(q$drop[q$station_id == refst])
  expect_output(print(q), "Flagged stations")

  expect_message(bruce_qc(d, drop = TRUE), "dropped 1")
  d2 <- suppressMessages(bruce_qc(d, drop = TRUE))
  expect_false(bad %in% d2$obs$station_id)
  expect_true(refst %in% d2$obs$station_id)
  expect_s3_class(d2$qc, "bruce_qc")
})

test_that("bruce_qc flags stations with few days and works without background", {
  sim <- bruce_simulate(n_days = 8, seed = 4)
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
  biased <- d$sources$source[d$sources$bias][1]
  st <- d$obs$station_id[d$obs$source == biased][2]
  keep <- d$obs$station_id != st | d$obs$date <= min(d$obs$date) + 1
  d$obs <- d$obs[keep, ]
  q <- bruce_qc(d)
  expect_match(q$flag[q$station_id == st], "only 2 day")
  expect_true(q$drop[q$station_id == st])
  d["background"] <- list(NULL)
  q2 <- bruce_qc(d)
  expect_equal(attr(q2, "against"), "reference stations' daily mean")
  expect_error(bruce_qc(d, threshold = 0), "positive")
})

test_that("plot methods draw without error", {
  sim <- bruce_simulate(n_days = 3, seed = 1)
  d <- bruce_data(sim$obs, sim$stations, sim$sources, crs = sim$crs, background = sim$background)
  g <- bruce_grid(sim$domain, 2000)
  r <- terra::rast(g, nlyrs = 3)
  terra::values(r) <- rnorm(terra::ncell(r) * 3, 25)
  names(r) <- as.character(seq(min(d$obs$date), by = 1, length.out = 3))
  pr <- structure(list(mean = r, sd = abs(r - 25), dates = as.Date(names(r)), cells = seq_len(terra::ncell(r)),
                       n_samples = 10), class = "bruce_prediction")
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  expect_invisible(plot(pr))
  expect_invisible(plot(pr, what = "sd", dates = names(r)[1:2], main = "x"))
  expect_invisible(plot(pr, scale = "free"))
  expect_invisible(plot(pr, scale = "anomaly"))
  expect_error(plot(pr, what = "sd", scale = "anomaly"), "anomaly")
  out <- bruce_plot_series(pr, d)
  expect_equal(nrow(out), 3)
  cvp <- tibble::tibble(fold = "a", station_id = rep(c("s1", "s2"), each = 5), date = Sys.Date(),
                        obs = rnorm(10), pred = rnorm(10), lower = -2, upper = 2)
  cv <- structure(list(predictions = cvp, scheme = "lso"), class = "bruce_cv")
  expect_invisible(plot(cv))
  grDevices::dev.off()
  expect_error(plot(structure(list(mean = r, sd = NULL), class = "bruce_prediction"), what = "sd"), "n_samples")
})
