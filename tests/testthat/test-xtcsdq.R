test_that("xtcsdq returns correct structure", {
  set.seed(123)
  n <- 8; tt <- 20
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = rnorm(n * tt),
    x1   = rnorm(n * tt)
  )
  res <- xtcsdq(y ~ x1, data = dat, index = c("id", "time"),
                quantiles = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "xtcsdq")
  expect_equal(res$N, n)
  expect_equal(res$TT, tt)
  expect_length(res$T_tau, 3)
  expect_length(res$pval_T, 3)
})

test_that("xtcsdq portmanteau stat is mean of T_tau", {
  set.seed(42)
  n <- 6; tt <- 15
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = rnorm(n * tt),
    x1   = rnorm(n * tt)
  )
  res <- xtcsdq(y ~ x1, data = dat, index = c("id", "time"),
                quantiles = c(0.25, 0.75))
  expect_equal(res$M_K, mean(res$T_tau), tolerance = 1e-10)
})

test_that("xtcsdq p-values are in [0,1]", {
  set.seed(7)
  n <- 5; tt <- 12
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = rnorm(n * tt),
    x1   = rnorm(n * tt)
  )
  res <- xtcsdq(y ~ x1, data = dat, index = c("id", "time"),
                quantiles = 0.5)
  expect_true(all(res$pval_T >= 0 & res$pval_T <= 1))
  expect_true(all(res$pval_Ttilde >= 0 & res$pval_Ttilde <= 1))
})

test_that("print and summary do not error", {
  set.seed(1)
  n <- 5; tt <- 12
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = rnorm(n * tt),
    x1   = rnorm(n * tt)
  )
  res <- xtcsdq(y ~ x1, data = dat, index = c("id", "time"),
                quantiles = c(0.25, 0.5))
  expect_output(print(res))
  expect_output(summary(res))
})
