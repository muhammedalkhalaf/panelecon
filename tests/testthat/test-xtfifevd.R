# Tests for xtfifevd package

test_that("xtfifevd estimates time-varying coefficients correctly", {
  # Simulate data with known coefficients
  set.seed(123)
  N <- 100
  T <- 10
  n <- N * T
  
  id <- rep(1:N, each = T)
  time <- rep(1:T, N)
  alpha_i <- rep(rnorm(N, sd = 1), each = T)
  z <- rep(rnorm(N), each = T)
  x <- rnorm(n)
  
  # True model: y = 1 + 2*x + 0.5*z + alpha_i + eps
  y <- 1 + 2 * x + 0.5 * z + alpha_i + rnorm(n, sd = 0.5)
  
  data <- data.frame(id = id, time = time, y = y, x = x, z = z)
  
  # FEVD estimation
  fit <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time")
  
  # Check that time-varying coefficient is close to true value

  expect_equal(fit$beta["x"], c(x = 2), tolerance = 0.2)
  
  # Check that time-invariant coefficient is estimated
  expect_equal(fit$gamma["z"], c(z = 0.5), tolerance = 0.3)
  
  # Check structure

  expect_s3_class(fit, "xtfifevd")
  expect_equal(fit$method, "FEVD")
  expect_equal(fit$N, n)
  expect_equal(fit$N_g, N)
})


test_that("FEF produces same point estimates as FEVD", {
  set.seed(456)
  N <- 50
  T <- 8
  n <- N * T
  
  id <- rep(1:N, each = T)
  time <- rep(1:T, N)
  alpha_i <- rep(rnorm(N), each = T)
  z <- rep(runif(N, -1, 1), each = T)
  x <- rnorm(n)
  
  y <- 0.5 + 1.5 * x + 0.8 * z + alpha_i + rnorm(n, sd = 0.3)
  data <- data.frame(id = id, time = time, y = y, x = x, z = z)
  
  fit_fevd <- fevd(y ~ x | z, data = data, id = "id", time = "time")
  fit_fef <- fef(y ~ x | z, data = data, id = "id", time = "time")
  
  # Point estimates should be identical (Pesaran-Zhou Proposition 3)
  expect_equal(fit_fevd$beta, fit_fef$beta, tolerance = 1e-10)
  expect_equal(fit_fevd$gamma, fit_fef$gamma, tolerance = 1e-10)
})


test_that("FEF-IV requires instruments", {
  data <- data.frame(
    id = rep(1:10, each = 5),
    time = rep(1:5, 10),
    y = rnorm(50),
    x = rnorm(50),
    z = rep(rnorm(10), each = 5)
  )
  
  expect_error(
    fef_iv(y ~ x | z, data = data, id = "id", time = "time"),
    "instruments"
  )
})


test_that("FEF-IV works with valid instruments", {
  set.seed(789)
  N <- 60
  T <- 6
  n <- N * T
  
  id <- rep(1:N, each = T)
  time <- rep(1:T, N)
  
  # Generate instrument correlated with z but not with error
  iv <- rep(rnorm(N), each = T)
  z <- iv + rep(rnorm(N, sd = 0.5), each = T)  # z correlated with iv
  z <- rep(tapply(z, id, mean), each = T)      # Make time-invariant
  
  alpha_i <- rep(rnorm(N), each = T)
  x <- rnorm(n)
  y <- 1 + x + 0.5 * z + alpha_i + rnorm(n, sd = 0.4)
  
  data <- data.frame(id = id, time = time, y = y, x = x, z = z, iv = iv)
  
  fit <- fef_iv(y ~ x | z, data = data, id = "id", time = "time",
                instruments = ~ iv)
  
  expect_s3_class(fit, "xtfifevd")
  expect_equal(fit$method, "FEF-IV")
  expect_true(!is.null(fit$gamma))
})


test_that("summary.xtfifevd produces correct output", {
  set.seed(101)
  N <- 40
  T <- 5
  
  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = rnorm(N * T),
    x = rnorm(N * T),
    z = rep(rnorm(N), each = T)
  )
  
  fit <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time")
  s <- summary(fit)
  
  expect_s3_class(s, "summary.xtfifevd")
  expect_true("coefficients" %in% names(s))
  expect_equal(nrow(s$coefficients), 3)  # x, z, _cons
})


test_that("vcov returns correct structure", {
  set.seed(202)
  N <- 30
  T <- 4
  
  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = rnorm(N * T),
    x1 = rnorm(N * T),
    x2 = rnorm(N * T),
    z = rep(rnorm(N), each = T)
  )
  
  fit <- xtfifevd(y ~ x1 + x2 | z, data = data, id = "id", time = "time")
  V <- vcov(fit)
  
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(V))
  expect_equal(nrow(V), 4)  # x1, x2, z, _cons
  
  # Check symmetry
  expect_equal(V, t(V))
  
  # Check positive definiteness (all eigenvalues positive)
  expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
})


test_that("bw_ratio computes correctly", {
  set.seed(303)
  N <- 20
  T <- 6
  
  data <- data.frame(
    id = rep(1:N, each = T),
    z_const = rep(rnorm(N), each = T),   # Time-invariant
    z_vary = rnorm(N * T)                 # Time-varying
  )
  
  result <- bw_ratio(data, c("z_const", "z_vary"), id = "id")
  
  # Time-invariant variable should have within SD near 0
  expect_lt(result$sd_within[result$variable == "z_const"], 1e-10)
  
  # Time-varying variable should have finite ratio
  expect_true(is.finite(result$bw_ratio[result$variable == "z_vary"]))
})


test_that("confint produces valid intervals", {
  set.seed(404)
  N <- 50
  T <- 5
  
  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = rnorm(N * T),
    x = rnorm(N * T),
    z = rep(rnorm(N), each = T)
  )
  
  fit <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time")
  ci <- confint(fit)
  
  expect_equal(nrow(ci), 3)
  expect_equal(ncol(ci), 2)
  
  # Lower bound should be less than upper bound
  expect_true(all(ci[, 1] < ci[, 2]))
  
  # Check that coefficients are within intervals
  coefs <- coef(fit)
  for (i in seq_along(coefs)) {
    expect_gte(coefs[i], ci[i, 1])
    expect_lte(coefs[i], ci[i, 2])
  }
})


test_that("formula parsing handles multiple variables", {
  set.seed(505)
  N <- 40
  T <- 5
  n <- N * T
  
  data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    z1 = rep(rnorm(N), each = T),
    z2 = rep(rnorm(N), each = T)
  )
  
  fit <- xtfifevd(y ~ x1 + x2 + x3 | z1 + z2,
                  data = data, id = "id", time = "time")
  
  expect_equal(fit$k_x, 3)
  expect_equal(fit$k_z, 2)
  expect_equal(length(fit$coefficients), 6)  # 3 x + 2 z + 1 cons
})
