test_that("grunfeld_cmg() returns correct structure", {
  dat <- grunfeld_cmg()
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 100L)
  expect_named(dat, c("firm", "year", "invest", "mvalue"))
})

test_that("xtpcmg Group-Mean FM-OLS (q=2) runs without error", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "mg", q = 2L)
  expect_s3_class(res, "xtpcmg")
  expect_equal(res$model, "mg")
  expect_equal(res$q, 2L)
  expect_equal(length(res$coefficients), 2L)
  expect_true(all(is.finite(res$coefficients)))
  expect_true(all(res$se > 0))
  expect_true(all(res$pvalue >= 0) && all(res$pvalue <= 1))
})

test_that("xtpcmg Pooled FM-OLS (q=2) runs without error", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "pmg", q = 2L)
  expect_s3_class(res, "xtpcmg")
  expect_equal(res$model, "pmg")
  expect_equal(length(res$coefficients), 2L)
})

test_that("xtpcmg cubic (q=3) runs without error", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "mg", q = 3L)
  expect_s3_class(res, "xtpcmg")
  expect_equal(res$q, 3L)
  expect_equal(length(res$coefficients), 3L)
})

test_that("xtpcmg turning point computed for q=2", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "mg", q = 2L)
  if (!is.na(res$tp)) {
    expect_true(is.finite(res$tp))
    expect_true(is.finite(res$tp_se))
    expect_true(res$tp_lo <= res$tp_hi)
  }
})

test_that("xtpcmg MG model stores individual coefficients", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "mg", q = 2L)
  expect_false(is.null(res$ind_coef))
  expect_equal(nrow(res$ind_coef), 5L)
  expect_equal(ncol(res$ind_coef), 2L)
})

test_that("xtpcmg Swamy test computed for MG model", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "mg", q = 2L)
  if (!is.na(res$swamy_s)) {
    expect_true(res$swamy_s >= 0)
    expect_true(res$swamy_p >= 0 && res$swamy_p <= 1)
  }
})

test_that("print.xtpcmg produces output without error", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "mg", q = 2L)
  expect_invisible(print(res))
})

test_that("summary.xtpcmg produces output without error", {
  dat <- grunfeld_cmg()
  res <- xtpcmg(dat, y = "invest", x = "mvalue",
                panel_id = "firm", time_id = "year",
                model = "mg", q = 2L)
  expect_invisible(summary(res))
})

test_that("xtpcmg errors on invalid q", {
  dat <- grunfeld_cmg()
  expect_error(
    xtpcmg(dat, y = "invest", x = "mvalue",
           panel_id = "firm", time_id = "year",
           model = "mg", q = 4L)
  )
})

test_that("xtpcmg errors on missing variable", {
  dat <- grunfeld_cmg()
  expect_error(
    xtpcmg(dat, y = "nonexistent", x = "mvalue",
           panel_id = "firm", time_id = "year",
           model = "mg", q = 2L)
  )
})
