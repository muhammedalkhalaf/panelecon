test_that("xtqsh returns correct structure", {
  set.seed(123)
  n <- 8; tt <- 15
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = rnorm(n * tt),
    x1   = rnorm(n * tt)
  )
  res <- xtqsh(y ~ x1, data = dat, index = c("id", "time"), tau = c(0.25, 0.5, 0.75))
  expect_s3_class(res, "xtqsh")
  expect_length(res$tau, 3)
  expect_length(res$S, 3)
  expect_length(res$D, 3)
  expect_true(all(!is.na(res$pval_D)))
  expect_true(all(res$pval_D >= 0 & res$pval_D <= 1))
})

test_that("xtqsh beta_md has correct dimensions", {
  set.seed(42)
  n <- 6; tt <- 20
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = rnorm(n * tt),
    x1   = rnorm(n * tt),
    x2   = rnorm(n * tt)
  )
  res <- xtqsh(y ~ x1 + x2, data = dat, index = c("id", "time"), tau = 0.5)
  expect_equal(dim(res$beta_md), c(1L, 2L))
  expect_equal(colnames(res$beta_md), c("x1", "x2"))
})

test_that("xtqsh marginal tests work", {
  set.seed(7)
  n <- 6; tt <- 15
  dat <- data.frame(
    id   = rep(1:n, each = tt),
    time = rep(1:tt, times = n),
    y    = rnorm(n * tt),
    x1   = rnorm(n * tt)
  )
  res <- xtqsh(y ~ x1, data = dat, index = c("id", "time"),
               tau = c(0.25, 0.75), marginal = TRUE)
  expect_false(is.null(res$marginal))
  expect_true("x1" %in% names(res$marginal))
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
  res <- xtqsh(y ~ x1, data = dat, index = c("id", "time"), tau = 0.5)
  expect_output(print(res))
  expect_output(summary(res))
})
