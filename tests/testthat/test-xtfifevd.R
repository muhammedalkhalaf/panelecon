test_that("xtfifevd FEVD returns correct structure", {
  set.seed(123)
  n <- 10; tt <- 15
  uid  <- rep(1:n, each = tt)
  tval <- rep(1:tt, times = n)
  z_i  <- rep(rnorm(n), each = tt)
  x_it <- rnorm(n * tt)
  y    <- 0.5 * x_it + 0.8 * z_i + rnorm(n * tt, sd = 0.5)
  dat  <- data.frame(id = uid, time = tval, y = y, x = x_it, z = z_i)

  res <- xtfifevd(y ~ x, data = dat, index = c("id", "time"),
                  zinvariants = "z", method = "fevd")
  expect_s3_class(res, "xtfifevd")
  expect_equal(res$method, "FEVD")
  expect_length(res$beta_fe, 1)
  expect_length(res$gamma, 1)
  expect_true(!is.null(res$delta))
  # delta should be close to 1
  expect_true(abs(res$delta - 1) < 0.5)
})

test_that("xtfifevd FEF returns correct structure", {
  set.seed(42)
  n <- 8; tt <- 12
  uid  <- rep(1:n, each = tt)
  tval <- rep(1:tt, times = n)
  z_i  <- rep(rnorm(n), each = tt)
  x_it <- rnorm(n * tt)
  y    <- 0.3 * x_it + 1.2 * z_i + rnorm(n * tt)
  dat  <- data.frame(id = uid, time = tval, y = y, x = x_it, z = z_i)

  res <- xtfifevd(y ~ x, data = dat, index = c("id", "time"),
                  zinvariants = "z", method = "fef")
  expect_s3_class(res, "xtfifevd")
  expect_equal(res$method, "FEF")
  expect_true(all(res$se_gamma > 0))
})

test_that("xtfifevd beta_fe is close to true value", {
  set.seed(7)
  n <- 20; tt <- 30
  uid  <- rep(1:n, each = tt)
  z_i  <- rep(rnorm(n), each = tt)
  x_it <- rnorm(n * tt)
  y    <- 2.0 * x_it + 1.0 * z_i + rnorm(n * tt, sd = 0.3)
  dat  <- data.frame(id = uid, time = rep(1:tt, times = n),
                     y = y, x = x_it, z = z_i)
  res  <- xtfifevd(y ~ x, data = dat, index = c("id", "time"),
                   zinvariants = "z", method = "fef")
  expect_equal(res$beta_fe, 2.0, tolerance = 0.2)
})

test_that("print and summary do not error", {
  set.seed(1)
  n <- 6; tt <- 10
  uid  <- rep(1:n, each = tt)
  z_i  <- rep(rnorm(n), each = tt)
  x_it <- rnorm(n * tt)
  y    <- x_it + z_i + rnorm(n * tt)
  dat  <- data.frame(id = uid, time = rep(1:tt, n), y = y, x = x_it, z = z_i)
  res  <- xtfifevd(y ~ x, data = dat, index = c("id", "time"),
                   zinvariants = "z")
  expect_output(print(res))
  expect_output(summary(res))
})
