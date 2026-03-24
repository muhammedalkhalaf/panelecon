make_panel <- function(n = 5, t = 10, seed = 1) {
  set.seed(seed)
  data.frame(
    id   = rep(seq_len(n), each = t),
    time = rep(seq_len(t), times = n),
    y    = rnorm(n * t),
    x1   = rnorm(n * t),
    x2   = rnorm(n * t)
  )
}

test_that("xtpretest runs Hsiao tests", {
  df  <- make_panel()
  res <- xtpretest(df, y ~ x1, index = c("id", "time"), tests = "HSIAO")
  expect_type(res, "list")
  expect_named(res$hsiao, c("F1", "F1p", "F1_nd", "F1_dd",
                             "F2", "F2p", "F2_nd", "F2_dd",
                             "F3", "F3p", "F3_nd", "F3_dd",
                             "S1", "S2", "S3"), ignore.order = TRUE)
  expect_true(is.numeric(res$hsiao$F1))
  expect_true(res$hsiao$F1p >= 0 && res$hsiao$F1p <= 1)
})

test_that("xtpretest runs CSD test", {
  df  <- make_panel(seed = 5)
  res <- xtpretest(df, y ~ x1 + x2, index = c("id", "time"), tests = "CSD")
  expect_named(res$csd, c("cd_stat", "pval", "avg_rho"), ignore.order = TRUE)
  expect_true(is.numeric(res$csd$cd_stat))
})

test_that("xtpretest returns a recommendation", {
  df  <- make_panel()
  res <- xtpretest(df, y ~ x1, index = c("id", "time"),
                   tests = c("HSIAO", "CSD"))
  expect_true(is.character(res$recommendation))
  expect_true(nchar(res$recommendation) > 0L)
})

test_that("xtpretest errors on bad index", {
  df <- make_panel()
  expect_error(xtpretest(df, y ~ x1, index = c("bad_id", "time")))
})

test_that("xtpretest summary module works", {
  df  <- make_panel()
  res <- xtpretest(df, y ~ x1, index = c("id", "time"), tests = "SUMMARY")
  expect_type(res$summary, "list")
})
