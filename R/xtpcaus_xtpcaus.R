#' Panel Granger Causality Tests
#'
#' Tests whether \code{x} Granger-causes \code{y} in a balanced panel using
#' either the Panel Fourier Toda-Yamamoto (PFTY) test or the Panel Quantile
#' Causality (PQC) test.
#'
#' @param data A data frame in long format.
#' @param y Character. Name of the dependent (caused) variable.
#' @param x Character. Name of the independent (causing) variable.
#' @param panel_id Character. Name of the panel identifier variable.
#' @param time_id Character. Name of the time variable.
#' @param test Character. Test type: \code{"pfty"} for Panel Fourier
#'   Toda-Yamamoto or \code{"pqc"} for Panel Quantile Causality.
#' @param pmax Integer. Maximum lag order for selection. Default is \code{4}.
#' @param dmax Integer. Maximum integration order for Toda-Yamamoto
#'   augmentation. Default is \code{1}.
#' @param nboot Integer. Number of bootstrap replications. Minimum 99.
#'   Default is \code{499}.
#' @param kmax Integer. Maximum Fourier frequency (PFTY only). Default is
#'   \code{3}.
#' @param ic Character. Information criterion: \code{"aic"} or \code{"bic"}.
#'   Default is \code{"aic"}.
#' @param quantiles Numeric vector. Quantile grid for PQC test (values strictly
#'   between 0 and 1). Default is
#'   \code{c(0.1, 0.25, 0.50, 0.75, 0.90)}.
#' @param seed Integer. Random seed for bootstrap. \code{-1} means no seed.
#'   Default is \code{-1}.
#'
#' @return An object of class \code{"xtpcaus"} containing:
#' \describe{
#'   \item{test}{Character. \code{"pfty"} or \code{"pqc"}.}
#'   \item{N}{Integer. Number of panel units.}
#'   \item{TT}{Integer. Number of time periods.}
#'   \item{nboot}{Integer. Number of bootstrap replications.}
#'   \item{y}{Character. Name of the y variable.}
#'   \item{x}{Character. Name of the x variable.}
#'   For PFTY:
#'   \item{fisher}{Numeric. Fisher panel statistic.}
#'   \item{fisher_df}{Integer. Degrees of freedom (2*N).}
#'   \item{fisher_pv}{Numeric. Fisher p-value.}
#'   \item{wbar}{Numeric. Average individual Wald statistic.}
#'   \item{zbar}{Numeric. Dumitrescu-Hurlin Z-bar statistic.}
#'   \item{zbar_pv}{Numeric. Z-bar p-value.}
#'   \item{ind_wald}{Numeric vector. Individual Wald statistics (length N).}
#'   \item{ind_freq}{Integer vector. Optimal Fourier frequencies (length N).}
#'   \item{ind_pval_b}{Numeric vector. Bootstrap p-values (length N).}
#'   \item{ind_lags}{Integer vector. Selected lag orders (length N).}
#'   For PQC:
#'   \item{quantiles}{Numeric vector. Quantiles tested.}
#'   \item{wald_xy}{Numeric vector. Wald statistics per quantile (x => y).}
#'   \item{pval_xy}{Numeric vector. Bootstrap p-values per quantile (x => y).}
#'   \item{wald_yx}{Numeric vector. Wald statistics per quantile (y => x).}
#'   \item{pval_yx}{Numeric vector. Bootstrap p-values per quantile (y => x).}
#'   \item{supwald_xy}{Numeric. Sup-Wald statistic for x => y.}
#'   \item{supwald_yx}{Numeric. Sup-Wald statistic for y => x.}
#'   \item{p_opt}{Integer. Selected optimal lag.}
#' }
#'
#' @references
#' Chuang, C.C., Kuan, C.M. and Lin, H.Y. (2009).
#' Causality in quantiles and dynamic stock return-volume relations.
#' \emph{Journal of Banking and Finance}, 33(7), 1351--1360.
#' \doi{10.1016/j.jbankfin.2009.02.013}
#'
#' Emirmahmutoglu, F. and Kose, N. (2011).
#' Testing for Granger causality in heterogeneous mixed panels.
#' \emph{Economic Modelling}, 28(3), 870--876.
#' \doi{10.1016/j.econmod.2010.10.018}
#'
#' Toda, H.Y. and Yamamoto, T. (1995).
#' Statistical inference in vector autoregressions with possibly integrated
#' processes. \emph{Journal of Econometrics}, 66(1--2), 225--250.
#' \doi{10.1016/0304-4076(94)01616-8}
#'
#' Wang, K.M. and Nguyen, T.B. (2022).
#' A quantile panel-type analysis of income inequality and healthcare
#' expenditure. \emph{Economic Research}, 35(1), 873--893.
#' \doi{10.1080/1331677X.2021.1952089}
#'
#' Yilanci, V. and Gorus, M.S. (2020).
#' Does economic globalization have predictive power for ecological footprint.
#' \emph{Environmental Science and Pollution Research}, 27, 40552--40562.
#' \doi{10.1007/s11356-020-09895-x}
#'
#' @examples
#' dat <- grunfeld_panel()
#' # PFTY test (quick with few bootstrap reps)
#' \donttest{
#' res <- xtpcaus(dat, y = "invest", x = "mvalue",
#'                panel_id = "firm", time_id = "year",
#'                test = "pfty", pmax = 2L, dmax = 1L,
#'                nboot = 99L, kmax = 2L, seed = 42L)
#' print(res)
#' }
#'
#' # PQC test
#' \donttest{
#' res2 <- xtpcaus(dat, y = "invest", x = "mvalue",
#'                 panel_id = "firm", time_id = "year",
#'                 test = "pqc", pmax = 2L, nboot = 99L,
#'                 quantiles = c(0.25, 0.50, 0.75), seed = 42L)
#' print(res2)
#' }
#'
#' @export
xtpcaus <- function(data, y, x, panel_id, time_id,
                    test = c("pfty", "pqc"),
                    pmax = 4L, dmax = 1L, nboot = 499L,
                    kmax = 3L, ic = c("aic", "bic"),
                    quantiles = c(0.1, 0.25, 0.50, 0.75, 0.90),
                    seed = -1L) {

  test <- match.arg(test)
  ic   <- match.arg(ic)

  # Validate inputs
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (v in c(y, x, panel_id, time_id)) {
    if (!v %in% names(data)) stop(sprintf("Variable '%s' not found in data.", v))
  }
  if (nboot < 99L) stop("'nboot' must be at least 99.")
  if (any(quantiles <= 0) || any(quantiles >= 1)) {
    stop("All 'quantiles' values must be strictly between 0 and 1.")
  }

  # Sort and check balanced panel
  data <- data[order(data[[panel_id]], data[[time_id]]), ]
  panels <- sort(unique(data[[panel_id]]))
  N <- length(panels)
  TT <- length(sort(unique(data[[time_id]])))

  if (N < 2L) stop("At least 2 panel units are required.")
  if (TT < 10L) stop("At least 10 time periods are required.")
  if (nrow(data) != N * TT) stop("Panel must be strongly balanced.")

  if (!is.na(seed) && seed >= 0L) set.seed(seed)

  Y_mat <- matrix(data[[y]], nrow = TT, ncol = N)
  X_mat <- matrix(data[[x]], nrow = TT, ncol = N)

  pmax <- as.integer(pmax)
  dmax <- as.integer(dmax)
  nboot <- as.integer(nboot)
  kmax <- as.integer(kmax)

  if (test == "pfty") {
    res <- .xtpcaus_pfty(Y_mat, X_mat, N, TT, pmax, dmax, kmax, ic, nboot)
  } else {
    res <- .xtpcaus_pqc(Y_mat, X_mat, N, TT, pmax, dmax, quantiles, nboot)
  }

  res$test      <- test
  res$N         <- N
  res$TT        <- TT
  res$nboot     <- nboot
  res$y         <- y
  res$x         <- x
  res$panel_id  <- panel_id
  res$time_id   <- time_id

  structure(res, class = "xtpcaus")
}


# ============================================================
# PFTY: Panel Fourier Toda-Yamamoto
# ============================================================

#' @keywords internal
.xtpcaus_select_lag_freq <- function(y_vec, x_vec, TT, pmax, kmax,
                                     dmax, ic_fn) {
  best_ic  <- Inf
  best_p   <- 1L
  best_f   <- 1L

  for (p in seq_len(pmax)) {
    for (f in seq_len(kmax)) {
      k_aug <- p + dmax
      n_reg <- TT - k_aug
      if (n_reg < p + 3L) next

      t_idx <- seq(k_aug + 1L, TT)
      # Fourier terms
      sin_f <- sin(2 * pi * f * t_idx / TT)
      cos_f <- cos(2 * pi * f * t_idx / TT)

      # Build regressors: const, lags of y (k_aug), lags of x (k_aug), sin, cos
      X_reg <- matrix(1, nrow = n_reg, ncol = 1L + 2L * k_aug + 2L)
      for (lag in seq_len(k_aug)) {
        X_reg[, 1L + lag]           <- y_vec[t_idx - lag]
        X_reg[, 1L + k_aug + lag]   <- x_vec[t_idx - lag]
      }
      X_reg[, 1L + 2L * k_aug + 1L] <- sin_f
      X_reg[, 1L + 2L * k_aug + 2L] <- cos_f

      y_dep <- y_vec[t_idx]
      beta <- tryCatch(
        solve(crossprod(X_reg), crossprod(X_reg, y_dep)),
        error = function(e) NULL
      )
      if (is.null(beta)) next
      resid <- y_dep - X_reg %*% beta
      sigma2 <- sum(resid^2) / n_reg
      ic_val <- ic_fn(sigma2, ncol(X_reg), n_reg)

      if (ic_val < best_ic) {
        best_ic <- ic_val
        best_p  <- p
        best_f  <- f
      }
    }
  }
  list(p = best_p, f = best_f)
}

#' @keywords internal
.xtpcaus_ic_fn <- function(ic) {
  if (ic == "aic") {
    function(sigma2, k, n) log(sigma2) + 2 * k / n
  } else {
    function(sigma2, k, n) log(sigma2) + log(n) * k / n
  }
}

#' @keywords internal
.xtpcaus_wald_stat <- function(y_vec, x_vec, TT, p, f, dmax) {
  k_aug <- p + dmax
  n_reg <- TT - k_aug
  if (n_reg < p + 3L) return(list(wald = NA_real_, resid = NULL))

  t_idx <- seq(k_aug + 1L, TT)
  sin_f <- sin(2 * pi * f * t_idx / TT)
  cos_f <- cos(2 * pi * f * t_idx / TT)

  n_cols <- 1L + 2L * k_aug + 2L
  X_reg <- matrix(1, nrow = n_reg, ncol = n_cols)
  for (lag in seq_len(k_aug)) {
    X_reg[, 1L + lag]           <- y_vec[t_idx - lag]
    X_reg[, 1L + k_aug + lag]   <- x_vec[t_idx - lag]
  }
  X_reg[, 1L + 2L * k_aug + 1L] <- sin_f
  X_reg[, 1L + 2L * k_aug + 2L] <- cos_f

  y_dep <- y_vec[t_idx]
  XtX <- crossprod(X_reg)
  XtX_inv <- tryCatch(solve(XtX), error = function(e) NULL)
  if (is.null(XtX_inv)) return(list(wald = NA_real_, resid = NULL))

  beta <- XtX_inv %*% crossprod(X_reg, y_dep)
  resid <- y_dep - X_reg %*% beta
  sigma2 <- sum(resid^2) / n_reg

  # Test H0: first p lags of x (cols 2+k_aug to 1+k_aug+p) = 0
  x_lag_cols <- seq(2L + k_aug, 1L + k_aug + p)
  R_mat <- matrix(0, nrow = p, ncol = n_cols)
  for (jj in seq_len(p)) R_mat[jj, x_lag_cols[jj]] <- 1

  Rbeta <- R_mat %*% beta
  RVR   <- R_mat %*% XtX_inv %*% t(R_mat)
  wald  <- as.numeric(t(Rbeta) %*% solve(sigma2 * RVR) %*% Rbeta)
  list(wald = wald, resid = as.numeric(resid), beta = as.numeric(beta),
       sigma2 = sigma2, n_reg = n_reg, t_idx = t_idx, X_reg = X_reg,
       p = p, f = f, k_aug = k_aug, x_lag_cols = x_lag_cols)
}

#' @keywords internal
.xtpcaus_pfty <- function(Y_mat, X_mat, N, TT, pmax, dmax, kmax, ic, nboot) {
  ic_fn <- .xtpcaus_ic_fn(ic)

  ind_wald   <- numeric(N)
  ind_freq   <- integer(N)
  ind_lags   <- integer(N)
  ind_pval_a <- numeric(N)
  ind_pval_b <- numeric(N)
  resid_list <- vector("list", N)
  wald_info  <- vector("list", N)

  for (i in seq_len(N)) {
    sel <- .xtpcaus_select_lag_freq(Y_mat[, i], X_mat[, i],
                                    TT, pmax, kmax, dmax, ic_fn)
    wi  <- .xtpcaus_wald_stat(Y_mat[, i], X_mat[, i], TT,
                               sel$p, sel$f, dmax)
    ind_wald[i]   <- if (is.na(wi$wald)) 0 else wi$wald
    ind_freq[i]   <- sel$f
    ind_lags[i]   <- sel$p
    ind_pval_a[i] <- 1 - stats::pchisq(ind_wald[i], df = sel$p)
    resid_list[[i]] <- wi$resid
    wald_info[[i]]  <- wi
  }

  # Bootstrap p-values (sieve bootstrap under H0)
  boot_wald <- matrix(0, nrow = nboot, ncol = N)
  for (b in seq_len(nboot)) {
    for (i in seq_len(N)) {
      wi <- wald_info[[i]]
      if (is.null(wi$resid)) { boot_wald[b, i] <- 0; next }
      # Resample residuals
      n_r   <- wi$n_reg
      idx_b <- sample.int(n_r, n_r, replace = TRUE)
      e_b   <- wi$resid[idx_b]
      # Reconstruct y under H0 (null: x lags zero)
      beta_null <- wi$beta
      beta_null[wi$x_lag_cols] <- 0
      y_null <- wi$X_reg %*% beta_null + e_b
      # Compute Wald on bootstrapped data
      wi_b <- .xtpcaus_wald_stat_direct(
        y_null, wi$X_reg, wi$p, wi$k_aug, wi$x_lag_cols
      )
      boot_wald[b, i] <- if (is.na(wi_b)) 0 else wi_b
    }
  }

  for (i in seq_len(N)) {
    ind_pval_b[i] <- mean(boot_wald[, i] >= ind_wald[i])
    # Avoid exact 0
    ind_pval_b[i] <- max(ind_pval_b[i], 1 / (nboot + 1))
  }

  # Fisher panel statistic
  fisher_stat <- -2 * sum(log(ind_pval_b))
  fisher_df   <- 2L * N
  fisher_pv   <- 1 - stats::pchisq(fisher_stat, df = fisher_df)

  # Dumitrescu-Hurlin statistics
  wbar <- mean(ind_wald)
  K    <- mean(ind_lags)
  zbar <- sqrt(N / (2 * K)) * (wbar - K)
  zbar_pv <- 2 * (1 - stats::pnorm(abs(zbar)))

  list(
    fisher = fisher_stat, fisher_df = fisher_df, fisher_pv = fisher_pv,
    wbar = wbar, zbar = zbar, zbar_pv = zbar_pv,
    ind_wald = ind_wald, ind_freq = ind_freq,
    ind_pval_a = ind_pval_a, ind_pval_b = ind_pval_b,
    ind_lags = ind_lags,
    pmax = pmax, dmax = dmax, kmax = kmax, ic = ic
  )
}

#' @keywords internal
.xtpcaus_wald_stat_direct <- function(y_dep, X_reg, p, k_aug, x_lag_cols) {
  n_reg  <- length(y_dep)
  n_cols <- ncol(X_reg)
  XtX_inv <- tryCatch(solve(crossprod(X_reg)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NA_real_)
  beta  <- XtX_inv %*% crossprod(X_reg, y_dep)
  resid <- y_dep - X_reg %*% beta
  sigma2 <- sum(resid^2) / n_reg
  R_mat <- matrix(0, nrow = p, ncol = n_cols)
  for (jj in seq_len(p)) R_mat[jj, x_lag_cols[jj]] <- 1
  Rbeta <- R_mat %*% beta
  RVR   <- R_mat %*% XtX_inv %*% t(R_mat)
  wald  <- as.numeric(t(Rbeta) %*% solve(sigma2 * RVR) %*% Rbeta)
  wald
}


# ============================================================
# PQC: Panel Quantile Causality
# ============================================================

#' @keywords internal
.xtpcaus_qreg_coef <- function(y, X, tau) {
  # Barrodale-Roberts quantile regression (L1)
  n  <- length(y)
  k  <- ncol(X)
  # Use stats::optim with Koenker-Bassett loss for small problems
  # For efficiency, implement direct iteratively-reweighted approach
  # Start from OLS
  beta0 <- tryCatch(
    as.numeric(solve(crossprod(X), crossprod(X, y))),
    error = function(e) rep(0, k)
  )
  # IRLS for quantile regression
  for (iter in seq_len(50L)) {
    resid <- as.numeric(y - X %*% beta0)
    w <- pmax(abs(resid), 1e-6)
    W <- diag(1 / w)
    # Tilted weight
    wt <- ifelse(resid >= 0, tau / w, (1 - tau) / w)
    beta1 <- tryCatch(
      as.numeric(solve(t(X) %*% W %*% X, t(X) %*% (W %*% y + wt * w - wt * resid))),
      error = function(e) beta0
    )
    if (max(abs(beta1 - beta0)) < 1e-8) { beta0 <- beta1; break }
    beta0 <- beta1
  }
  beta0
}

#' @keywords internal
.xtpcaus_pqc_wald <- function(y_vec, x_vec, N, TT, p, tau) {
  # Panel quantile VAR with individual fixed effects (within transformation)
  # Stack: outcome = y, regressors = lags of y + lags of x + individual dummies
  n_reg <- TT - p
  if (n_reg < p + 2L) return(list(wald = 0, coef_sum = 0))

  t_idx <- seq(p + 1L, TT)

  # Build stacked (N * n_reg) x (2p + N) design matrix
  n_obs <- N * n_reg
  n_cols <- 2L * p + N  # p lags y + p lags x + N fixed effects
  X_all <- matrix(0, nrow = n_obs, ncol = n_cols)
  y_all <- numeric(n_obs)

  for (i in seq_len(N)) {
    rows_i <- seq((i - 1L) * n_reg + 1L, i * n_reg)
    y_all[rows_i] <- y_vec[t_idx + (i - 1L) * TT]
    for (lag in seq_len(p)) {
      X_all[rows_i, lag]     <- y_vec[t_idx - lag + (i - 1L) * TT]
      X_all[rows_i, p + lag] <- x_vec[t_idx - lag + (i - 1L) * TT]
    }
    X_all[rows_i, 2L * p + i] <- 1  # panel fixed effect
  }

  beta_q <- tryCatch(
    .xtpcaus_qreg_coef(y_all, X_all, tau),
    error = function(e) NULL
  )
  if (is.null(beta_q)) return(list(wald = 0, coef_sum = 0))

  # Extract x-lag coefficients (cols p+1 to 2p)
  x_lag_cols <- seq(p + 1L, 2L * p)
  beta_x <- beta_q[x_lag_cols]

  # Wald test: H0: beta_x = 0
  resid <- y_all - X_all %*% beta_q
  sigma2 <- sum(resid^2) / (n_obs - ncol(X_all))
  XtX_inv <- tryCatch(solve(crossprod(X_all)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(list(wald = 0, coef_sum = sum(beta_x)))

  R_mat <- matrix(0, nrow = p, ncol = n_cols)
  for (jj in seq_len(p)) R_mat[jj, x_lag_cols[jj]] <- 1
  RVR   <- R_mat %*% XtX_inv %*% t(R_mat)
  wald  <- tryCatch(
    as.numeric(t(beta_x) %*% solve(sigma2 * RVR) %*% beta_x),
    error = function(e) 0
  )
  list(wald = max(wald, 0), coef_sum = sum(beta_x))
}

#' @keywords internal
.xtpcaus_pqc <- function(Y_mat, X_mat, N, TT, pmax, dmax, quantiles, nboot) {
  # Select optimal lag via AIC (OLS)
  y_vec <- as.vector(Y_mat)
  x_vec <- as.vector(X_mat)

  best_aic <- Inf
  p_opt <- 1L
  for (p in seq_len(pmax)) {
    n_reg <- TT - p
    if (n_reg < p + 2L) next
    t_idx <- seq(p + 1L, TT)
    n_obs_p <- N * n_reg
    X_ls <- matrix(0, nrow = n_obs_p, ncol = 2L * p + N)
    y_ls <- numeric(n_obs_p)
    for (i in seq_len(N)) {
      rows_i <- seq((i - 1L) * n_reg + 1L, i * n_reg)
      y_ls[rows_i] <- Y_mat[t_idx, i]
      for (lag in seq_len(p)) {
        X_ls[rows_i, lag]     <- Y_mat[t_idx - lag, i]
        X_ls[rows_i, p + lag] <- X_mat[t_idx - lag, i]
      }
      X_ls[rows_i, 2L * p + i] <- 1
    }
    beta_ols <- tryCatch(
      solve(crossprod(X_ls), crossprod(X_ls, y_ls)),
      error = function(e) NULL
    )
    if (is.null(beta_ols)) next
    rss <- sum((y_ls - X_ls %*% beta_ols)^2)
    aic_val <- log(rss / n_obs_p) + 2 * ncol(X_ls) / n_obs_p
    if (aic_val < best_aic) { best_aic <- aic_val; p_opt <- p }
  }

  nq <- length(quantiles)
  wald_xy  <- numeric(nq)
  coef_xy  <- numeric(nq)
  wald_yx  <- numeric(nq)
  coef_yx  <- numeric(nq)
  pval_xy  <- numeric(nq)
  pval_yx  <- numeric(nq)

  # Vectorized Y and X for stacking
  y_stacked <- as.vector(Y_mat)
  x_stacked <- as.vector(X_mat)

  for (q_idx in seq_len(nq)) {
    tau <- quantiles[q_idx]
    res_xy <- .xtpcaus_pqc_wald(y_stacked, x_stacked, N, TT, p_opt, tau)
    res_yx <- .xtpcaus_pqc_wald(x_stacked, y_stacked, N, TT, p_opt, tau)
    wald_xy[q_idx] <- res_xy$wald
    coef_xy[q_idx] <- res_xy$coef_sum
    wald_yx[q_idx] <- res_yx$wald
    coef_yx[q_idx] <- res_yx$coef_sum
  }

  # Bootstrap p-values
  n_reg <- TT - p_opt
  t_idx_b <- seq(p_opt + 1L, TT)
  n_obs_b <- N * n_reg
  X_bs <- matrix(0, nrow = n_obs_b, ncol = 2L * p_opt + N)
  y_bs <- numeric(n_obs_b)
  for (i in seq_len(N)) {
    rows_i <- seq((i - 1L) * n_reg + 1L, i * n_reg)
    y_bs[rows_i] <- Y_mat[t_idx_b, i]
    for (lag in seq_len(p_opt)) {
      X_bs[rows_i, lag]           <- Y_mat[t_idx_b - lag, i]
      X_bs[rows_i, p_opt + lag]   <- X_mat[t_idx_b - lag, i]
    }
    X_bs[rows_i, 2L * p_opt + i] <- 1
  }

  # OLS under H0 (restrict x-lag coefs to 0)
  beta_ols_full <- tryCatch(
    as.numeric(solve(crossprod(X_bs), crossprod(X_bs, y_bs))),
    error = function(e) rep(0, ncol(X_bs))
  )
  beta_null <- beta_ols_full
  x_cols_null <- seq(p_opt + 1L, 2L * p_opt)
  beta_null[x_cols_null] <- 0
  resid_null <- y_bs - X_bs %*% beta_null

  # Bootstrap
  boot_wald_xy <- matrix(0, nrow = nboot, ncol = nq)
  boot_wald_yx <- matrix(0, nrow = nboot, ncol = nq)

  for (b in seq_len(nboot)) {
    idx_b <- sample.int(n_obs_b, n_obs_b, replace = TRUE)
    e_b   <- resid_null[idx_b]
    y_boot <- X_bs %*% beta_null + e_b
    y_boot_vec <- as.numeric(y_boot)

    for (q_idx in seq_len(nq)) {
      tau <- quantiles[q_idx]
      wb <- tryCatch(
        .xtpcaus_pqc_wald(y_boot_vec, x_stacked, N, TT, p_opt, tau)$wald,
        error = function(e) 0
      )
      boot_wald_xy[b, q_idx] <- if (is.na(wb)) 0 else wb
      # For y=>x direction, bootstrap uses same y (x causes y direction uses x as dependent)
      wb2 <- tryCatch(
        .xtpcaus_pqc_wald(x_stacked, y_boot_vec, N, TT, p_opt, tau)$wald,
        error = function(e) 0
      )
      boot_wald_yx[b, q_idx] <- if (is.na(wb2)) 0 else wb2
    }
  }

  for (q_idx in seq_len(nq)) {
    pval_xy[q_idx] <- max(mean(boot_wald_xy[, q_idx] >= wald_xy[q_idx]),
                          1 / (nboot + 1))
    pval_yx[q_idx] <- max(mean(boot_wald_yx[, q_idx] >= wald_yx[q_idx]),
                          1 / (nboot + 1))
  }

  supwald_xy <- max(wald_xy)
  supwald_yx <- max(wald_yx)

  list(
    quantiles = quantiles, p_opt = p_opt,
    wald_xy = wald_xy, coef_xy = coef_xy, pval_xy = pval_xy,
    wald_yx = wald_yx, coef_yx = coef_yx, pval_yx = pval_yx,
    supwald_xy = supwald_xy, supwald_yx = supwald_yx
  )
}


# ============================================================
# S3 Methods
# ============================================================

#' Print Method for xtpcaus Objects
#'
#' @param x An object of class \code{"xtpcaus"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtpcaus <- function(x, ...) {
  message(strrep("-", 65))
  message(sprintf("Panel Causality Test: %s", toupper(x$test)))
  message(sprintf("H0: %s does NOT Granger-cause %s", x$x, x$y))
  message(sprintf("N = %d panels, T = %d periods, nboot = %d",
    x$N, x$TT, x$nboot))
  message(strrep("-", 65))

  if (x$test == "pfty") {
    message(sprintf("Fisher stat: %.4f  df: %d  p-value: %.4f%s",
      x$fisher, x$fisher_df, x$fisher_pv,
      if (x$fisher_pv < 0.01) "***" else if (x$fisher_pv < 0.05) "**" else
      if (x$fisher_pv < 0.10) "*" else ""))
    message(sprintf("W-bar:       %.4f", x$wbar))
    message(sprintf("Z-bar:       %.4f  p-value: %.4f%s",
      x$zbar, x$zbar_pv,
      if (x$zbar_pv < 0.01) "***" else if (x$zbar_pv < 0.05) "**" else
      if (x$zbar_pv < 0.10) "*" else ""))
    message("")
    message("Individual results:")
    message(sprintf("  %-6s %-6s %-8s %-12s %-12s",
      "Panel", "Lag", "Freq", "Wald", "Boot p-val"))
    for (i in seq_len(x$N)) {
      stars <- if (x$ind_pval_b[i] < 0.01) "***" else
                if (x$ind_pval_b[i] < 0.05) "**" else
                if (x$ind_pval_b[i] < 0.10) "*" else ""
      message(sprintf("  %-6d %-6d %-8d %12.4f %8.4f%s",
        i, x$ind_lags[i], x$ind_freq[i],
        x$ind_wald[i], x$ind_pval_b[i], stars))
    }
  } else {
    message(sprintf("Sup-Wald (%s => %s): %.4f", x$x, x$y, x$supwald_xy))
    message(sprintf("Sup-Wald (%s => %s): %.4f", x$y, x$x, x$supwald_yx))
    message(sprintf("Optimal lag: %d", x$p_opt))
    message("")
    message(sprintf("  %-6s %-12s %-10s %-12s %-10s",
      "Tau", paste0("Wald(", x$x, "=>", x$y, ")"),
      "p-val", paste0("Wald(", x$y, "=>", x$x, ")"), "p-val"))
    for (q_idx in seq_along(x$quantiles)) {
      message(sprintf("  %-6.2f %12.4f %10.4f%s %12.4f %10.4f%s",
        x$quantiles[q_idx],
        x$wald_xy[q_idx], x$pval_xy[q_idx],
        if (x$pval_xy[q_idx] < 0.05) "*" else " ",
        x$wald_yx[q_idx], x$pval_yx[q_idx],
        if (x$pval_yx[q_idx] < 0.05) "*" else " "))
    }
  }
  message(strrep("-", 65))
  message("*** p<0.01, ** p<0.05, * p<0.10")
  invisible(x)
}

#' Summary Method for xtpcaus Objects
#'
#' @param object An object of class \code{"xtpcaus"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtpcaus <- function(object, ...) {
  print(object)
  invisible(object)
}

#' Example Panel Data for xtpcaus
#'
#' Returns a small balanced panel dataset (subset of Grunfeld 1958) for
#' use in examples and testing.
#'
#' @return A data frame with columns \code{firm}, \code{year},
#'   \code{invest}, and \code{mvalue}.
#'
#' @examples
#' dat <- grunfeld_panel()
#' head(dat)
#'
#' @export
grunfeld_panel <- function() {
  data.frame(
    firm = rep(1:5, each = 20L),
    year = rep(1935L:1954L, times = 5L),
    invest = c(
      317.6, 391.8, 410.6, 257.7, 330.8, 461.2, 512.0, 448.0, 499.6, 547.5,
      561.2, 688.1, 568.9, 529.2, 555.1, 642.9, 755.9, 891.2, 1304.4, 1486.7,
      40.29, 72.76, 66.26, 65.67, 76.04, 87.99, 100.0, 95.1, 104.4, 118.2,
      114.0, 135.7, 130.7, 169.6, 162.7, 162.0, 190.2, 181.9, 232.8, 256.7,
      209.9, 355.3, 318.9, 267.4, 339.5, 400.3, 419.5, 417.3, 347.2, 364.2,
      361.1, 312.2, 273.1, 264.0, 187.0, 145.8, 208.4, 162.8, 184.8, 152.1,
      33.1, 45.0, 77.2, 44.6, 48.1, 74.4, 113.0, 91.9, 61.3, 56.8,
      93.6, 159.9, 147.2, 146.3, 98.3, 93.5, 135.2, 157.3, 179.5, 189.6,
      40.29, 72.76, 66.26, 65.67, 76.04, 87.99, 100.0, 95.1, 104.4, 118.2,
      114.0, 135.7, 130.7, 169.6, 162.7, 162.0, 190.2, 181.9, 232.8, 256.7
    ),
    mvalue = c(
      2792.7, 2759.9, 2132.0, 1834.1, 1588.0, 1749.4, 1687.2, 2007.7, 2177.3, 2011.6,
      2533.2, 2459.4, 2157.7, 2082.0, 1808.1, 1786.0, 1975.7, 2021.7, 2390.4, 2613.6,
      182.8, 213.3, 206.6, 209.0, 225.5, 226.7, 226.7, 222.3, 227.9, 237.9,
      226.7, 232.6, 262.1, 302.1, 324.3, 335.9, 351.9, 344.9, 365.5, 382.1,
      1362.4, 1807.1, 1952.1, 2080.0, 1968.1, 1795.5, 1666.6, 1634.8, 1640.0, 1671.7,
      1547.5, 1560.0, 1401.4, 1364.0, 1216.8, 1099.8, 1131.3, 1169.7, 1289.5, 1467.1,
      95.3, 107.5, 136.3, 113.4, 141.1, 168.6, 214.8, 218.0, 186.0, 230.3,
      228.7, 293.2, 289.0, 271.7, 246.0, 239.6, 281.0, 265.3, 270.8, 290.5,
      182.8, 213.3, 206.6, 209.0, 225.5, 226.7, 226.7, 222.3, 227.9, 237.9,
      226.7, 232.6, 262.1, 302.1, 324.3, 335.9, 351.9, 344.9, 365.5, 382.1
    ),
    stringsAsFactors = FALSE
  )
}
