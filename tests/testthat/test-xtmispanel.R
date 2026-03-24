make_panel_na <- function(seed = 1) {
  set.seed(seed)
  df <- data.frame(
    id   = rep(1:4, each = 8),
    time = rep(1:8, times = 4),
    y    = rnorm(32),
    x1   = rnorm(32)
  )
  df$y[c(3, 11, 20, 27)] <- NA
  df
}

test_that("xtmispanel detect module runs without error", {
  df  <- make_panel_na()
  res <- xtmispanel(df, vars = c("y", "x1"), index = c("id", "time"), detect = TRUE)
  expect_type(res, "list")
  expect_named(res, "detect", ignore.order = TRUE)
  expect_equal(res$detect$total_missing, 4L)
})

test_that("xtmispanel imputation with mean works", {
  df  <- make_panel_na()
  res <- xtmispanel(df, vars = "y", index = c("id", "time"),
                    detect = FALSE, impute = "mean", target = "y")
  expect_false(is.null(res$imputed))
  expect_false(anyNA(res$imputed$y_imp))
})

test_that("xtmispanel imputation with locf works", {
  df  <- make_panel_na()
  res <- xtmispanel(df, vars = "y", index = c("id", "time"),
                    detect = FALSE, impute = "locf", target = "y")
  expect_named(res$imputed, c("id", "time", "y", "x1", "y_imp"), ignore.order = TRUE)
})

test_that("xtmispanel imputation with linear interpolation works", {
  df  <- make_panel_na()
  res <- xtmispanel(df, vars = "y", index = c("id", "time"),
                    detect = FALSE, impute = "linear", target = "y")
  expect_false(anyNA(res$imputed$y_imp))
})

test_that("xtmispanel impute_stats reports correlation", {
  df  <- make_panel_na()
  res <- xtmispanel(df, vars = "y", index = c("id", "time"),
                    detect = FALSE, impute = "mean", target = "y")
  expect_true(!is.null(res$impute_stats$correlation))
})

test_that("xtmispanel errors on missing index columns", {
  df <- make_panel_na()
  expect_error(xtmispanel(df, vars = "y", index = c("bad", "time")))
})

test_that("xtmispanel sensitivity runs and returns best_method", {
  df  <- make_panel_na()
  res <- xtmispanel(df, vars = "y", index = c("id", "time"),
                    detect = FALSE, sensitivity = TRUE, target = "y")
  expect_true(!is.null(res$sensitivity$best_method))
  expect_true(res$sensitivity$best_method %in%
                c("mean","median","locf","nocb","linear","spline",
                  "pmm","hotdeck","knn","rf","em"))
})

test_that(".xmp_locf fills forward correctly", {
  x   <- c(1, NA, NA, 4)
  res <- xtmispanel:::.xmp_locf(x)
  expect_equal(res, c(1, 1, 1, 4))
})
