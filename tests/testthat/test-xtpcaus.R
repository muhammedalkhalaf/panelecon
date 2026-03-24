test_that("grunfeld_panel() returns correct structure", {
  dat <- grunfeld_panel()
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 100L)
  expect_named(dat, c("firm", "year", "invest", "mvalue"))
})

test_that("xtpcaus PFTY runs without error", {
  dat <- grunfeld_panel()
  res <- xtpcaus(dat, y = "invest", x = "mvalue",
                 panel_id = "firm", time_id = "year",
                 test = "pfty", pmax = 2L, dmax = 1L,
                 nboot = 99L, kmax = 2L, seed = 42L)
  expect_s3_class(res, "xtpcaus")
  expect_equal(res$test, "pfty")
  expect_true(is.finite(res$fisher))
  expect_true(res$fisher_pv >= 0 && res$fisher_pv <= 1)
  expect_equal(length(res$ind_wald), 5L)
  expect_equal(length(res$ind_pval_b), 5L)
})

test_that("xtpcaus PQC runs without error", {
  dat <- grunfeld_panel()
  res <- xtpcaus(dat, y = "invest", x = "mvalue",
                 panel_id = "firm", time_id = "year",
                 test = "pqc", pmax = 2L, nboot = 99L,
                 quantiles = c(0.25, 0.50, 0.75), seed = 42L)
  expect_s3_class(res, "xtpcaus")
  expect_equal(res$test, "pqc")
  expect_equal(length(res$wald_xy), 3L)
  expect_equal(length(res$pval_xy), 3L)
  expect_true(all(res$pval_xy >= 0) && all(res$pval_xy <= 1))
  expect_true(is.finite(res$supwald_xy))
})

test_that("xtpcaus errors on invalid quantiles", {
  dat <- grunfeld_panel()
  expect_error(
    xtpcaus(dat, y = "invest", x = "mvalue",
            panel_id = "firm", time_id = "year",
            test = "pqc", quantiles = c(0, 0.5))
  )
  expect_error(
    xtpcaus(dat, y = "invest", x = "mvalue",
            panel_id = "firm", time_id = "year",
            test = "pqc", quantiles = c(0.5, 1.0))
  )
})

test_that("xtpcaus errors on insufficient nboot", {
  dat <- grunfeld_panel()
  expect_error(
    xtpcaus(dat, y = "invest", x = "mvalue",
            panel_id = "firm", time_id = "year",
            test = "pfty", nboot = 10L)
  )
})

test_that("print.xtpcaus produces output for PFTY", {
  dat <- grunfeld_panel()
  res <- xtpcaus(dat, y = "invest", x = "mvalue",
                 panel_id = "firm", time_id = "year",
                 test = "pfty", pmax = 1L, dmax = 1L,
                 nboot = 99L, kmax = 1L, seed = 1L)
  expect_invisible(print(res))
})

test_that("print.xtpcaus produces output for PQC", {
  dat <- grunfeld_panel()
  res <- xtpcaus(dat, y = "invest", x = "mvalue",
                 panel_id = "firm", time_id = "year",
                 test = "pqc", pmax = 1L, nboot = 99L,
                 quantiles = c(0.25, 0.75), seed = 1L)
  expect_invisible(print(res))
})
