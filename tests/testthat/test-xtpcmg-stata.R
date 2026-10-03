## Reference values from the Stata module xtpcmg 1.0.2 (validated by its
## author against Wagner's MATLAB code) on the data generated below.
.pcmg_data <- function() {
  set.seed(20260928)
  do.call(rbind, lapply(1:10, function(i) {
    x1 <- cumsum(rnorm(50)); x2 <- cumsum(rnorm(50))
    u <- as.numeric(arima.sim(list(ar = 0.5), 50))
    y <- 1 + i / 10 + (0.5 + i / 50) * x1 - 0.3 * x2 + u
    data.frame(id = i, t = 1:50, y = y, x1 = x1, x2 = x2)
  }))
}

test_that("xtpcmg reproduces Stata xtpcmg", {
  d <- .pcmg_data()
  f <- function(...) xtpcmg(d, y = "y", x = "x1", panel_id = "id", time_id = "t", ...)
  chk <- function(r, b, se) {
    expect_equal(round(unname(r$coefficients), 6), b)
    expect_equal(round(unname(sqrt(diag(r$vcov))), 6), se)
  }
  chk(f(model = "mg"), c(0.717803, -0.042017), c(0.103348, 0.014335))
  chk(f(model = "pmg"), c(0.540448, 0.000337), c(0.036554, 0.003309))
  chk(f(model = "pmg", effects = "twoway"), c(0.535173, 0.002419), c(0.033325, 0.003327))
  chk(f(model = "mg", q = 3), c(0.068321, 0.100415, -0.015023), c(0.247626, 0.073499, 0.007351))
  chk(f(model = "pmg", controls = "x2"), c(0.576552, 0.001667, -0.267848), c(0.026270, 0.002428, 0.018113))
  chk(f(model = "pmg", kernel = "qs"), c(0.542464, 0.000562), c(0.038311, 0.003411))
  chk(f(model = "mg", corr_rob = TRUE), c(0.717803, -0.042017), c(0.100751, 0.014441))
})
