# Tests pinning the corrected computations to hand calculations
# (Pesaran and Zhou 2018, working paper equation numbers).

sim_panel <- function(seed = 1, N = 60, T = 5, unbal = FALSE) {
  set.seed(seed)
  eta <- rnorm(N)
  z1 <- rnorm(N) + 0.5
  z2 <- rbinom(N, 1, 0.4)
  r1 <- rnorm(N)
  r2 <- rnorm(N)
  Ti <- if (unbal) sample(2:T, N, TRUE) else rep(T, N)
  id <- rep(seq_len(N), Ti)
  tt <- unlist(lapply(Ti, seq_len))
  n <- length(id)
  x1 <- 1 + 0.7 * eta[id] + 0.5 * z1[id] + rnorm(n)
  x2 <- rnorm(n) + 0.3 * z2[id] + 2
  e <- rnorm(n, sd = 0.5 + z2[id])
  y <- 1 + 2 * x1 - x2 + z1[id] - 0.5 * z2[id] + eta[id] + e
  data.frame(id = id, time = tt, y = y, x1 = x1, x2 = x2,
             z1 = z1[id], z2 = z2[id], r1 = r1[id], r2 = r2[id])
}

# Hand computation of FEF / FEF-IV with the full Pesaran-Zhou covariance
byhand <- function(d, xn = c("x1", "x2"), zn = c("z1", "z2"), rn = NULL) {
  X <- as.matrix(d[, xn]); ids <- d$id
  xb <- apply(X, 2, function(v) ave(v, ids)); yb <- ave(d$y, ids)
  Xw <- X - xb; yw <- d$y - yb
  beta <- solve(crossprod(Xw), crossprod(Xw, yw))
  e <- as.vector(yw - Xw %*% beta)
  A <- solve(crossprod(Xw)); M <- 0
  for (i in unique(ids)) {
    s <- crossprod(Xw[ids == i, , drop = FALSE], e[ids == i])
    M <- M + s %*% t(s)
  }
  Vb <- A %*% M %*% A                                     # eq. (18)
  first <- !duplicated(ids); N <- sum(first)
  ybi <- yb[first]; xbi <- xb[first, , drop = FALSE]; Z <- as.matrix(d[first, zn])
  ui <- as.vector(ybi - xbi %*% beta)
  Zc <- scale(Z, scale = FALSE); Xc <- scale(xbi, scale = FALSE)
  zbar <- colMeans(Z); xbar <- colMeans(xbi)
  Qzz <- crossprod(Zc) / N; Qzx <- crossprod(Zc, Xc) / N
  if (is.null(rn)) {
    gam <- solve(crossprod(Zc), crossprod(Zc, ui - mean(ui)))   # eq. (4)
    Gz <- solve(Qzz) %*% t(Zc) / N
    Gb <- solve(Qzz) %*% Qzx
  } else {
    R <- as.matrix(d[first, rn]); Rc <- scale(R, scale = FALSE)
    Qzr <- crossprod(Zc, Rc) / N; Qrr <- crossprod(Rc) / N
    Qru <- crossprod(Rc, ui - mean(ui)) / N
    H <- solve(Qzr %*% solve(Qrr) %*% t(Qzr)) %*% Qzr %*% solve(Qrr)
    gam <- H %*% Qru                                            # eq. (48)
    Gz <- H %*% t(Rc) / N
    Gb <- H %*% (crossprod(Rc, Xc) / N)
  }
  alpha <- mean(ui) - sum(zbar * gam)                           # eq. (5)
  vs <- as.vector(ui - alpha - Z %*% gam)                       # eq. (20)
  Vgamma <- Gz %*% (vs^2 * t(Gz)) + Gb %*% Vb %*% t(Gb)       # eq. (17)/(51)
  # cross-covariance from eq. (A.11) and delta method for eq. (5)
  Cgb <- -Gb %*% Vb
  w <- as.vector(1 - N * crossprod(zbar, Gz))
  cvec <- as.vector(xbar - crossprod(Gb, zbar))
  Valpha <- sum(w^2 * vs^2) / N^2 + as.numeric(t(cvec) %*% Vb %*% cvec)
  Cab <- -t(cvec) %*% Vb
  Cag <- Gz %*% (w * vs^2) / N + Gb %*% Vb %*% cvec
  list(beta = drop(beta), gamma = drop(gam), alpha = alpha, Vb = Vb,
       Vgamma = Vgamma, Cgb = Cgb, Valpha = Valpha, Cab = drop(Cab),
       Cag = drop(Cag))
}

test_that("FEF matches the Pesaran-Zhou formulas (eq. 4, 5, 17, 18, A.11)", {
  for (ub in c(FALSE, TRUE)) {
    d <- sim_panel(seed = 10 + ub, unbal = ub)
    f <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time")
    h <- byhand(d)
    expect_equal(unname(f$beta), unname(h$beta), tolerance = 1e-12)
    expect_equal(unname(f$gamma), unname(h$gamma), tolerance = 1e-12)
    expect_equal(unname(f$intercept), h$alpha, tolerance = 1e-12)
    V <- vcov(f)
    expect_equal(unname(V[1:2, 1:2]), unname(h$Vb), tolerance = 1e-12)
    expect_equal(unname(V[3:4, 3:4]), unname(h$Vgamma), tolerance = 1e-12)
    expect_equal(unname(f$V_gamma_pz), unname(h$Vgamma), tolerance = 1e-12)
    # cross-covariances and intercept variance
    expect_equal(unname(V[3:4, 1:2]), unname(h$Cgb), tolerance = 1e-12)
    expect_equal(unname(V[5, 5]), h$Valpha, tolerance = 1e-12)
    expect_equal(unname(V[5, 1:2]), unname(h$Cab), tolerance = 1e-12)
    expect_equal(unname(V[5, 3:4]), unname(h$Cag), tolerance = 1e-12)
    expect_equal(V, t(V))
  }
})

test_that("FEF-IV matches eq. 48 and 51 with the full covariance", {
  d <- sim_panel(seed = 3, N = 80)
  d$z1 <- d$z1 + 0.6 * ave(d$y, d$id) / 10   # make z1 endogenous-looking
  f <- fef_iv(y ~ x1 + x2 | z1 + z2, d, "id", "time",
              instruments = ~ r1 + r2 + z2)
  h <- byhand(d, rn = c("r1", "r2", "z2"))
  expect_equal(unname(f$gamma), unname(h$gamma), tolerance = 1e-12)
  expect_equal(unname(f$intercept), h$alpha, tolerance = 1e-12)
  V <- vcov(f)
  expect_equal(unname(V[3:4, 3:4]), unname(h$Vgamma), tolerance = 1e-12)
  expect_equal(unname(V[3:4, 1:2]), unname(h$Cgb), tolerance = 1e-12)
  expect_equal(unname(V[5, 5]), h$Valpha, tolerance = 1e-12)
  expect_equal(unname(V[5, 3:4]), unname(h$Cag), tolerance = 1e-12)
})

test_that("classical beta vcov equals sigma2_e (X'MX)^-1 from the dummy regression", {
  d <- sim_panel(seed = 4)
  f <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time", vcov_beta = "classical")
  lsdv <- lm(y ~ x1 + x2 + factor(id), data = d)
  expect_equal(unname(vcov(f)[1:2, 1:2]), unname(vcov(lsdv)[2:3, 2:3]),
               tolerance = 1e-10)
  expect_equal(f$sigma2_e, summary(lsdv)$sigma^2, tolerance = 1e-10)
  expect_equal(f$vcov_beta, "classical")
  # the robust matrix is still stored
  f2 <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time")
  expect_equal(f$V_beta_robust, f2$V_beta_robust)
  expect_equal(unname(vcov(f2)[1:2, 1:2]), unname(f2$V_beta_robust))
})

test_that("robust beta vcov equals plm::vcovHC(method = 'arellano', type = 'HC0')", {
  skip_if_not_installed("plm")
  d <- sim_panel(seed = 5, unbal = TRUE)
  f <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time")
  pw <- plm::plm(y ~ x1 + x2, data = plm::pdata.frame(d, index = c("id", "time")),
                 model = "within")
  Vplm <- plm::vcovHC(pw, method = "arellano", type = "HC0")
  expect_equal(unname(vcov(f)[1:2, 1:2]), unname(Vplm[1:2, 1:2]),
               tolerance = 1e-10)
})

test_that("FEVD stage 3 gives delta = 1 and coefficients equal to FEF", {
  for (ub in c(FALSE, TRUE)) {
    d <- sim_panel(seed = 20 + ub, unbal = ub)
    fv <- fevd(y ~ x1 + x2 | z1 + z2, d, "id", "time")
    ff <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time")
    expect_equal(fv$delta, 1, tolerance = 1e-10)
    expect_equal(coef(fv), coef(ff), tolerance = 1e-10)
    expect_equal(vcov(fv), vcov(ff))
    expect_equal(fv$fef$coefficients, coef(ff))
    # literal three stages with lm()
    X <- as.matrix(d[, c("x1", "x2")])
    uh <- ave(d$y, d$id) - as.vector(apply(X, 2, ave, d$id) %*% ff$beta)
    first <- !duplicated(d$id)
    s2 <- lm(uh[first] ~ z1 + z2, data = d[first, ])
    d$h <- uh - predict(s2, d)
    s3 <- lm(y ~ x1 + x2 + z1 + z2 + h, data = d)
    expect_equal(unname(fv$stage3$coefficients), unname(coef(s3)), tolerance = 1e-10)
    expect_equal(unname(fv$stage3$se_naive), unname(sqrt(diag(vcov(s3)))),
                 tolerance = 1e-8)
    # naive stage 3 SEs are smaller than the Pesaran-Zhou ones
    expect_true(all(fv$stage3$se_naive[c("z1", "z2")] <
                      sqrt(diag(vcov(fv)))[c("z1", "z2")]))
  }
})

test_that("formula parser handles transformations, factors and log(y)", {
  d <- sim_panel(seed = 30)
  d$x1 <- exp(d$x1 / 4)
  d$y <- d$y - min(d$y) + 1
  d$g <- factor(rep(sample(c("a", "b", "c"), 60, TRUE), each = 5))
  f <- fef(log(y) ~ log(x1) + I(x2^2) + x1:x2 | z1 + g, d, "id", "time")
  expect_equal(f$x_names, c("log(x1)", "I(x2^2)", "x1:x2"))
  expect_equal(f$z_names, c("z1", "gb", "gc"))
  expect_equal(f$y_name, "log(y)")
  # within regression on the transformed data
  lsdv <- lm(log(y) ~ log(x1) + I(x2^2) + x1:x2 + factor(id), data = d)
  expect_equal(unname(f$beta), unname(coef(lsdv)[c("log(x1)", "I(x2^2)", "x1:x2")]),
               tolerance = 1e-10)
  # stage 2 with factor contrasts, one row per unit
  first <- !duplicated(d$id)
  uh <- ave(log(d$y), d$id) -
    as.vector(apply(cbind(log(d$x1), d$x2^2, d$x1 * d$x2), 2, ave, d$id) %*% f$beta)
  s2 <- lm(uh[first] ~ z1 + g, data = d[first, ])
  expect_equal(unname(f$gamma), unname(coef(s2)[2:4]), tolerance = 1e-10)
  expect_equal(unname(f$intercept), unname(coef(s2)[1]), tolerance = 1e-10)
})

test_that("factor id with unused levels and character id work", {
  d <- sim_panel(seed = 31)
  d$id <- factor(d$id)
  d$y[d$id == "7"] <- NA
  f <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time")
  expect_equal(f$N_g, 59)
  expect_equal(f$N, 59 * 5)
  d2 <- sim_panel(seed = 31)
  d2 <- d2[d2$id != 7, ]
  d2$id <- paste0("p", d2$id)
  f2 <- fef(y ~ x1 + x2 | z1 + z2, d2, "id", "time")
  expect_equal(f$gamma, f2$gamma, tolerance = 1e-12)
})

test_that("duplicated (id, time) rows are an error", {
  d <- sim_panel(seed = 32)
  d2 <- rbind(d, d[1:3, ])
  expect_error(fef(y ~ x1 + x2 | z1 + z2, d2, "id", "time"), "Duplicated")
})

test_that("rarely changing z uses the unit mean with a warning", {
  d <- sim_panel(seed = 33)
  d$z1 <- d$z1 + rnorm(nrow(d), sd = 0.2)
  expect_warning(f <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time"),
                 "not strictly time-invariant")
  d$z1m <- ave(d$z1, d$id)
  f2 <- fef(y ~ x1 + x2 | z1m + z2, d, "id", "time")
  expect_equal(unname(f$gamma), unname(f2$gamma), tolerance = 1e-12)
  expect_equal(unname(vcov(f)), unname(vcov(f2)), tolerance = 1e-12)
})

test_that("sigma2_u is the stage 2 residual variance net of sigma2_e/T", {
  d <- sim_panel(seed = 34)
  f <- fef(y ~ x1 + x2 | z1 + z2, d, "id", "time")
  s2 <- sum(f$stage2_residuals^2) / (f$N_g - f$k_z - 1)
  expect_equal(f$sigma2_u, max(s2 - f$sigma2_e / 5, 0), tolerance = 1e-12)
})

test_that("a formula with two top-level bars is an error", {
  d <- sim_panel(seed = 35)
  expect_error(fef(y ~ x1 + x2 | z1 | z2, d, "id", "time"), "exactly one")
  expect_error(fef(y ~ x1 | x2 | z1 + z2, d, "id", "time"), "exactly one")
})
