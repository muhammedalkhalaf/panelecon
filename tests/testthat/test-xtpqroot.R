test_that("grunfeld_pqroot() returns correct structure", {
  dat <- grunfeld_pqroot()
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 100L)
  expect_named(dat, c("firm", "year", "invest", "mvalue"))
})

test_that("xtpqroot CIPS(tau) runs without error", {
  dat <- grunfeld_pqroot()
  res <- xtpqroot(dat, var = "invest", panel_id = "firm",
                  time_id = "year", test = "cipstau",
                  model = "intercept", maxlag = 1L,
                  quantiles = c(0.25, 0.50, 0.75), reps = 100L)
  expect_s3_class(res, "xtpqroot")
  expect_equal(res$test, "cipstau")
  expect_true(is.finite(res$cips))
  expect_equal(length(res$cipstau), 3L)
  expect_equal(length(res$cipstau_pv), 3L)
  expect_true(all(res$cipstau_pv >= 0) && all(res$cipstau_pv <= 1))
})

test_that("xtpqroot errors with invalid quantiles", {
  dat <- grunfeld_pqroot()
  expect_error(
    xtpqroot(dat, var = "invest", panel_id = "firm",
             time_id = "year", test = "cipstau",
             quantiles = c(0, 0.5))
  )
})

test_that("xtpqroot errors with trendshift in cipstau", {
  dat <- grunfeld_pqroot()
  expect_error(
    xtpqroot(dat, var = "invest", panel_id = "firm",
             time_id = "year", test = "cipstau",
             model = "trendshift")
  )
})

test_that("xtpqroot errors when variable not found", {
  dat <- grunfeld_pqroot()
  expect_error(
    xtpqroot(dat, var = "nonexistent", panel_id = "firm",
             time_id = "year", test = "cipstau")
  )
})

test_that("print.xtpqroot CIPS(tau) produces output without error", {
  dat <- grunfeld_pqroot()
  res <- xtpqroot(dat, var = "invest", panel_id = "firm",
                  time_id = "year", test = "cipstau",
                  model = "intercept", maxlag = 1L,
                  quantiles = c(0.25, 0.75), reps = 50L)
  expect_invisible(print(res))
})
