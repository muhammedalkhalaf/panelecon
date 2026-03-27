#' Panel Quantile Unit Root Tests
#'
#' Tests for a panel unit root using either the CIPS(tau) quantile test
#' (Yang, Wei & Cai 2022) or the tFR Fourier-smooth-break test
#' (Corakci & Omay 2023).
#'
#' @param data A data frame in long format.
#' @param var Character. Name of the variable to test.
#' @param panel_id Character. Name of the panel identifier variable.
#' @param time_id Character. Name of the time variable.
#' @param test Character. Test type: \code{"cipstau"} for the CIPS(tau)
#'   quantile test (default) or \code{"tfr"} for the Fourier-LST test.
#' @param model Character. Deterministic terms: \code{"intercept"} (default)
#'   for intercept only, \code{"trend"} for intercept plus linear trend.
#'   For the tFR test, also accepts \code{"trendshift"} for Model C.
#' @param quantiles Numeric vector. Quantiles for CIPS(tau) test (values
#'   strictly between 0 and 1). Default is
#'   \code{c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)}.
#' @param maxlag Integer. Maximum lag for ADF augmentation. Default -1 means
#'   automatic selection.
#' @param reps Integer. Number of Monte Carlo replications for CIPS(tau)
#'   p-values. Default is \code{500}.
#' @param bootreps Integer. Number of sieve bootstrap replications for tFR
#'   p-values. Default is \code{500}.
#'
#' @return An object of class \code{"xtpqroot"} containing:
#' \describe{
#'   \item{test}{Character. \code{"cipstau"} or \code{"tfr"}.}
#'   \item{var}{Character. Variable name tested.}
#'   \item{N}{Integer. Number of panel units.}
#'   \item{TT}{Integer. Number of time periods.}
#'   \item{model}{Character. Deterministic specification.}
#'   For CIPS(tau):
#'   \item{cips}{Numeric. Standard OLS CIPS statistic.}
#'   \item{cips_pv}{Numeric. Monte Carlo p-value for CIPS.}
#'   \item{cipstau}{Numeric vector. CIPS(tau) statistics.}
#'   \item{cipstau_pv}{Numeric vector. Monte Carlo p-values.}
#'   \item{quantiles}{Numeric vector. Quantile grid used.}
#'   For tFR:
#'   \item{tfr}{Numeric. tFR panel statistic.}
#'   \item{pvalue}{Numeric. Bootstrap p-value.}
#'   \item{cv01}{Numeric. 1% bootstrap critical value.}
#'   \item{cv05}{Numeric. 5% bootstrap critical value.}
#'   \item{cv10}{Numeric. 10% bootstrap critical value.}
#'   \item{ind_results}{Data frame of per-panel results.}
#' }
#'
#' @references
#' Corakci, A. and Omay, T. (2023). Is there convergence in renewable energy
#' deployment? \emph{Renewable Energy}, 205, 648--662.
#' \doi{10.1016/j.renene.2023.01.060}
#'
#' Pesaran, M.H. (2007). A simple panel unit root test in the presence of
#' cross-section dependence. \emph{Journal of Applied Econometrics}, 22,
#' 265--312. \doi{10.1002/jae.951}
#'
#' Yang, Z., Wei, Z. and Cai, Y. (2022). Quantile unit root inference for
#' panel data with common shocks. \emph{Economics Letters}, 219, 110809.
#' \doi{10.1016/j.econlet.2022.110809}
#'
#' @examples
#' \donttest{
#' dat <- grunfeld_pqroot()
#'
#' # CIPS(tau) test at 3 quantiles
#' res <- xtpqroot(dat, var = "invest", panel_id = "firm",
#'                 time_id = "year", test = "cipstau",
#'                 model = "intercept", maxlag = 2L,
#'                 quantiles = c(0.1, 0.5, 0.9), reps = 100L)
#' print(res)
#'
#' # tFR test
#' res2 <- xtpqroot(dat, var = "invest", panel_id = "firm",
#'                  time_id = "year", test = "tfr",
#'                  model = "intercept", maxlag = 2L, bootreps = 200L)
#' print(res2)
#' }
#'
#' @export
xtpqroot <- function(data, var, panel_id, time_id,
                     test = c("cipstau", "tfr"),
                     model = c("intercept", "trend", "trendshift"),
                     quantiles = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                     maxlag = -1L, reps = 500L, bootreps = 500L) {

  test  <- match.arg(test)
  model <- match.arg(model)

  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (v in c(var, panel_id, time_id)) {
    if (!v %in% names(data)) stop(sprintf("Variable '%s' not found.", v))
  }
  if (test == "cipstau" && model == "trendshift") {
    stop("model = 'trendshift' is only available for test = 'tfr'.")
  }
  if (any(quantiles <= 0) || any(quantiles >= 1)) {
    stop("All 'quantiles' values must be strictly between 0 and 1.")
  }

  data <- data[order(data[[panel_id]], data[[time_id]]), ]
  panels <- sort(unique(data[[panel_id]]))
  N <- length(panels)
  TT <- length(sort(unique(data[[time_id]])))

  if (N < 3L) stop("At least 3 panel units are required.")
  if (TT < 15L) stop("At least 15 time periods are required.")
  if (nrow(data) != N * TT) stop("Panel must be strongly balanced.")

  Y <- matrix(data[[var]], nrow = TT, ncol = N)
  maxlag <- as.integer(maxlag)
  if (maxlag == -1L) {
    maxlag <- max(1L, floor(4 * (TT / 100)^(1 / 4)))
  }

  if (test == "cipstau") {
    res <- .xtpqroot_cipstau(Y, N, TT, model, quantiles, maxlag, reps)
  } else {
    res <- .xtpqroot_tfr(Y, N, TT, model, maxlag, bootreps)
  }

  res$test     <- test
  res$var      <- var
  res$N        <- N
  res$TT       <- TT
  res$model    <- model
  res$panel_id <- panel_id
  res$time_id  <- time_id

  structure(res, class = "xtpqroot")
}


# ============================================================
# CIPS(tau) Test
# ============================================================

#' @keywords internal
.xtpqroot_cadf_tau <- function(y_i, y_bar, dy_bar, p, model, tau) {
  # CADF(tau) for panel unit i at quantile tau
  TT <- length(y_i)
  n_reg <- TT - 1L - p

  if (n_reg < p + 3L) return(NA_real_)

  t_idx <- seq(p + 2L, TT)  # response indices

  # Build regressors
  # y_{i,t-1}, Delta y_{i,t-1},..., Delta y_{i,t-p},
  # y_bar_{t-1}, Delta y_bar_t, Delta y_bar_{t-1},..., Delta y_bar_{t-p}
  n_det <- if (model == "intercept") 1L else 2L
  n_col <- n_det + 1L + p + 1L + p + 1L  # intercept [+trend], y_{i,t-1}, dy lags, ybar_{t-1}, dybar lags

  X <- matrix(0, nrow = n_reg, ncol = n_col)
  col <- 1L
  X[, col] <- 1; col <- col + 1L
  if (model == "trend") { X[, col] <- t_idx; col <- col + 1L }
  X[, col] <- y_i[t_idx - 1L]; col <- col + 1L  # y_{i,t-1}
  dyi_full <- diff(y_i)
  for (lag in seq_len(p)) {
    idx_lag <- t_idx - 1L - lag
    if (all(idx_lag >= 1L)) X[, col] <- dyi_full[idx_lag]
    col <- col + 1L
  }
  X[, col] <- y_bar[t_idx - 1L]; col <- col + 1L  # y_bar_{t-1}
  X[, col] <- dy_bar[t_idx]; col <- col + 1L       # Delta y_bar_t
  for (lag in seq_len(p)) {
    if (col > n_col) break
    dy_bar_lag <- dy_bar[t_idx - lag]
    if (length(dy_bar_lag) == n_reg) { X[, col] <- dy_bar_lag; col <- col + 1L }
  }
  X <- X[, seq_len(col - 1L), drop = FALSE]

  y_dep <- diff(y_i)[t_idx - 1L]  # Delta y_{i,t}

  # Quantile regression t-ratio for y_{i,t-1} coefficient
  beta_q <- tryCatch(
    .xtpqroot_qreg(y_dep, X, tau),
    error = function(e) NULL
  )
  if (is.null(beta_q)) return(NA_real_)

  # t-statistic for rho coefficient (col 2 or 3 depending on model)
  rho_col <- if (model == "intercept") 2L else 3L
  rho_hat <- beta_q[rho_col]

  # Sparsity estimation for standard error
  resid <- y_dep - X %*% beta_q
  h_n <- stats::bw.nrd0(resid) * (n_reg^(-1/5))
  h_n <- max(h_n, 1e-6)
  f_hat <- mean(stats::dnorm(resid / h_n) / h_n)
  if (f_hat < 1e-8) f_hat <- 1e-8

  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NA_real_)

  se_rho <- sqrt(tau * (1 - tau) / (n_reg * f_hat^2) * XtX_inv[rho_col, rho_col])
  if (se_rho < 1e-10) return(NA_real_)

  (rho_hat - 1) / se_rho
}

#' @keywords internal
.xtpqroot_qreg <- function(y, X, tau) {
  # Simple quantile regression via IRLS
  n <- length(y)
  k <- ncol(X)
  beta <- tryCatch(
    as.numeric(solve(crossprod(X), crossprod(X, y))),
    error = function(e) rep(0, k)
  )
  for (iter in seq_len(50L)) {
    r <- as.numeric(y - X %*% beta)
    w <- pmax(abs(r), 1e-6)
    wt <- ifelse(r >= 0, tau / w, (1 - tau) / w)
    XWX <- t(X) %*% diag(1 / w) %*% X
    XWy <- t(X) %*% (diag(1 / w) %*% y + wt * w - wt * r)
    beta1 <- tryCatch(as.numeric(solve(XWX, XWy)), error = function(e) beta)
    if (max(abs(beta1 - beta)) < 1e-8) { beta <- beta1; break }
    beta <- beta1
  }
  beta
}

#' @keywords internal
.xtpqroot_cipstau <- function(Y, N, TT, model, quantiles, maxlag, reps) {
  # Cross-sectional averages
  y_bar  <- rowMeans(Y)
  dy_bar <- c(NA_real_, diff(y_bar))

  # Truncation constants (Pesaran 2007)
  K1 <- 6.19
  K2 <- if (model == "intercept") 2.16 else 2.61

  # Standard OLS CADF per panel
  cadf_stats <- numeric(N)
  for (i in seq_len(N)) {
    cadf_stats[i] <- .xtpqroot_cadf_ols(Y[, i], y_bar, dy_bar, maxlag, model)
  }
  cadf_trunc <- pmax(-K1, pmin(K2, cadf_stats))
  cips_stat  <- mean(cadf_trunc, na.rm = TRUE)

  nq <- length(quantiles)
  cipstau_stats <- matrix(NA_real_, nrow = N, ncol = nq)
  rho_tau       <- matrix(NA_real_, nrow = N, ncol = nq)

  for (q_idx in seq_len(nq)) {
    tau <- quantiles[q_idx]
    for (i in seq_len(N)) {
      cadf_q <- .xtpqroot_cadf_tau(Y[, i], y_bar, dy_bar, maxlag, model, tau)
      cipstau_stats[i, q_idx] <- if (is.na(cadf_q)) cadf_trunc[i] else cadf_q
      # Rho(tau)
      rho_tau[i, q_idx] <- .xtpqroot_rho_tau(Y[, i], y_bar, dy_bar,
                                               maxlag, model, tau)
    }
  }
  # Panel CIPS(tau) = cross-sectional mean of truncated CADF(tau)
  cipstau_trunc <- apply(cipstau_stats, 2, function(col) {
    mean(pmax(-K1, pmin(K2, col)), na.rm = TRUE)
  })

  # Monte Carlo p-values
  mc_dist_cips    <- .xtpqroot_mc_sim(N, TT, model, maxlag, quantiles, reps)
  cips_pv         <- mean(mc_dist_cips$cips <= cips_stat)
  cipstau_pv      <- sapply(seq_len(nq), function(q_idx) {
    mean(mc_dist_cips$cipstau[, q_idx] <= cipstau_trunc[q_idx])
  })

  list(
    cips = cips_stat, cips_pv = cips_pv,
    cipstau = cipstau_trunc, cipstau_pv = cipstau_pv,
    quantiles = quantiles,
    maxlag = maxlag, reps = reps,
    rho_tau = rho_tau,
    cadf_stats = cadf_stats
  )
}

#' @keywords internal
.xtpqroot_cadf_ols <- function(y_i, y_bar, dy_bar, p, model) {
  TT <- length(y_i)
  t_idx <- seq(p + 2L, TT)
  n_reg <- length(t_idx)
  if (n_reg < 3L) return(NA_real_)

  n_det <- if (model == "intercept") 1L else 2L
  X <- matrix(0, nrow = n_reg, ncol = n_det + 1L + p + 1L + p + 1L)
  col <- 1L
  X[, col] <- 1; col <- col + 1L
  if (model == "trend") { X[, col] <- t_idx; col <- col + 1L }
  X[, col] <- y_i[t_idx - 1L]; col <- col + 1L
  dy_i <- diff(y_i)
  for (lag in seq_len(p)) {
    idx_lag <- t_idx - 1L - lag
    if (all(idx_lag >= 1L))
      X[, col] <- dy_i[idx_lag]
    col <- col + 1L
  }
  X[, col] <- y_bar[t_idx - 1L]; col <- col + 1L
  X[, col] <- dy_bar[t_idx]; col <- col + 1L
  for (lag in seq_len(p)) {
    X[, col] <- dy_bar[pmax(t_idx - lag, 1L)]; col <- col + 1L
  }
  X <- X[, seq_len(col - 1L), drop = FALSE]
  y_dep <- dy_i[t_idx - 1L]

  beta <- tryCatch(solve(crossprod(X), crossprod(X, y_dep)),
                   error = function(e) NULL)
  if (is.null(beta)) return(NA_real_)

  resid  <- y_dep - X %*% beta
  sigma2 <- sum(resid^2) / (n_reg - ncol(X))
  if (sigma2 <= 0) return(NA_real_)

  rho_col <- if (model == "intercept") 2L else 3L
  se_rho  <- sqrt(sigma2 * solve(crossprod(X))[rho_col, rho_col])
  if (se_rho < 1e-10) return(NA_real_)
  (beta[rho_col] - 1) / se_rho
}

#' @keywords internal
.xtpqroot_rho_tau <- function(y_i, y_bar, dy_bar, p, model, tau) {
  TT <- length(y_i)
  t_idx <- seq(p + 2L, TT)
  n_reg <- length(t_idx)
  if (n_reg < 3L) return(NA_real_)

  n_det <- if (model == "intercept") 1L else 2L
  X <- matrix(0, nrow = n_reg, ncol = n_det + 1L + p + 1L + p + 1L)
  col <- 1L
  X[, col] <- 1; col <- col + 1L
  if (model == "trend") { X[, col] <- t_idx; col <- col + 1L }
  X[, col] <- y_i[t_idx - 1L]; col <- col + 1L
  dy_i <- diff(y_i)
  for (lag in seq_len(p)) {
    idx_lag <- t_idx - 1L - lag
    if (all(idx_lag >= 1L)) X[, col] <- dy_i[idx_lag]
    col <- col + 1L
  }
  X[, col] <- y_bar[t_idx - 1L]; col <- col + 1L
  X[, col] <- dy_bar[t_idx]; col <- col + 1L
  for (lag in seq_len(p)) {
    X[, col] <- dy_bar[pmax(t_idx - lag, 1L)]; col <- col + 1L
  }
  X <- X[, seq_len(col - 1L), drop = FALSE]
  y_dep <- dy_i[t_idx - 1L]

  beta_q <- tryCatch(.xtpqroot_qreg(y_dep, X, tau), error = function(e) NULL)
  if (is.null(beta_q)) return(NA_real_)

  rho_col <- if (model == "intercept") 2L else 3L
  1 + beta_q[rho_col]  # level rho = 1 + delta (coefficient on y_{t-1})
}

#' @keywords internal
.xtpqroot_mc_sim <- function(N, TT, model, maxlag, quantiles, reps) {
  nq  <- length(quantiles)
  K1  <- 6.19
  K2  <- if (model == "intercept") 2.16 else 2.61

  sim_cips    <- numeric(reps)
  sim_cipstau <- matrix(0, nrow = reps, ncol = nq)

  for (r in seq_len(reps)) {
    # Simulate random walk panel under H0
    eps <- matrix(stats::rnorm(N * TT), nrow = TT, ncol = N)
    Y_sim <- apply(eps, 2, cumsum)
    yb_sim  <- rowMeans(Y_sim)
    dyb_sim <- c(NA_real_, diff(yb_sim))

    cadf_r <- numeric(N)
    for (i in seq_len(N)) {
      cadf_r[i] <- .xtpqroot_cadf_ols(Y_sim[, i], yb_sim, dyb_sim, maxlag, model)
    }
    sim_cips[r] <- mean(pmax(-K1, pmin(K2, cadf_r), na.rm = TRUE))

    for (q_idx in seq_len(nq)) {
      tau <- quantiles[q_idx]
      cadf_qt <- numeric(N)
      for (i in seq_len(N)) {
        cv <- .xtpqroot_cadf_tau(Y_sim[, i], yb_sim, dyb_sim, maxlag, model, tau)
        cadf_qt[i] <- if (is.na(cv)) cadf_r[i] else cv
      }
      sim_cipstau[r, q_idx] <- mean(pmax(-K1, pmin(K2, cadf_qt)), na.rm = TRUE)
    }
  }

  list(cips = sim_cips, cipstau = sim_cipstau)
}


# ============================================================
# tFR Test: Fourier + Logistic Smooth Transition
# ============================================================

#' @keywords internal
.xtpqroot_lst <- function(t_vec, gamma, tau_loc, TT) {
  1 / (1 + exp(-gamma * (t_vec / TT - tau_loc)))
}

#' @keywords internal
.xtpqroot_tfr_unit <- function(y_i, TT, model, maxlag) {
  # Step 1: Search for optimal fractional Fourier frequency k_fr
  k_grid   <- seq(0.1, 5.0, by = 0.1)
  gamma_grid <- c(0.1, 0.5, 1, 3, 5, 10, 30, 50)
  tau_grid   <- seq(0.15, 0.85, by = 0.05)
  t_vec    <- seq_len(TT)

  best_ssr <- Inf
  best_kfr <- 1.0
  best_gam <- 5.0
  best_tau <- 0.5

  for (kfr in k_grid) {
    sin_k <- sin(pi * kfr * t_vec / TT)
    cos_k <- cos(pi * kfr * t_vec / TT)
    for (gam in gamma_grid) {
      for (tau_loc in tau_grid) {
        lst_t <- .xtpqroot_lst(t_vec, gam, tau_loc, TT)

        if (model == "intercept") {
          X_det <- cbind(1, sin_k, cos_k, lst_t)
        } else if (model == "trend") {
          X_det <- cbind(1, t_vec, sin_k, cos_k, lst_t)
        } else {
          X_det <- cbind(1, t_vec, sin_k, cos_k, lst_t, t_vec * lst_t)
        }

        beta_d <- tryCatch(
          as.numeric(solve(crossprod(X_det), crossprod(X_det, y_i))),
          error = function(e) NULL
        )
        if (is.null(beta_d)) next
        resid_d <- y_i - X_det %*% beta_d
        ssr <- sum(resid_d^2)
        if (ssr < best_ssr) {
          best_ssr <- ssr
          best_kfr <- kfr
          best_gam <- gam
          best_tau <- tau_loc
        }
      }
    }
  }

  # Step 2: ADF with Fourier + LST deterministics, BIC lag selection
  sin_k <- sin(pi * best_kfr * t_vec / TT)
  cos_k <- cos(pi * best_kfr * t_vec / TT)
  lst_t <- .xtpqroot_lst(t_vec, best_gam, best_tau, TT)

  if (model == "intercept") {
    det_cols <- cbind(1, sin_k, cos_k, lst_t)
  } else if (model == "trend") {
    det_cols <- cbind(1, t_vec, sin_k, cos_k, lst_t)
  } else {
    det_cols <- cbind(1, t_vec, sin_k, cos_k, lst_t, t_vec * lst_t)
  }

  best_bic <- Inf
  best_p   <- 0L
  best_tstat <- NA_real_

  for (p in seq(0L, maxlag)) {
    t_reg <- seq(p + 2L, TT)
    n_r   <- length(t_reg)
    if (n_r < p + 3L) next

    dy_i <- diff(y_i)
    y_lag1 <- y_i[t_reg - 1L]

    X_base <- det_cols[t_reg, , drop = FALSE]
    X_reg  <- cbind(X_base, y_lag1)
    if (p > 0L) {
      for (lag in seq_len(p)) X_reg <- cbind(X_reg, dy_i[t_reg - lag])
    }

    y_dep <- dy_i[t_reg - 1L]
    beta <- tryCatch(solve(crossprod(X_reg), crossprod(X_reg, y_dep)),
                     error = function(e) NULL)
    if (is.null(beta)) next

    resid  <- y_dep - X_reg %*% beta
    sigma2 <- sum(resid^2) / n_r
    bic_v  <- n_r * log(sigma2) + p * log(n_r)

    if (bic_v < best_bic) {
      best_bic <- bic_v
      best_p   <- as.integer(p)
      # t-statistic for y_{i,t-1} coefficient
      rho_col <- ncol(X_base) + 1L
      se_rho <- sqrt(sigma2 * solve(crossprod(X_reg))[rho_col, rho_col])
      best_tstat <- if (se_rho < 1e-10) NA_real_ else
        (beta[rho_col] - 1) / se_rho
    }
  }

  # Break date
  break_date_frac <- best_tau
  list(
    tstat = best_tstat, kfr = best_kfr, gamma = best_gam,
    tau = best_tau, p_opt = best_p,
    break_date_frac = break_date_frac
  )
}

#' @keywords internal
.xtpqroot_tfr <- function(Y, N, TT, model, maxlag, bootreps) {
  ind_res <- vector("list", N)
  for (i in seq_len(N)) {
    ind_res[[i]] <- .xtpqroot_tfr_unit(Y[, i], TT, model, maxlag)
  }

  tfr_stats <- sapply(ind_res, function(r) {
    if (is.na(r$tstat)) 0 else r$tstat
  })
  tfr_panel <- mean(tfr_stats)

  # Sieve bootstrap for panel-level distribution
  boot_stats <- numeric(bootreps)
  # Estimate AR(1) residuals for each panel
  ar_resids <- vector("list", N)
  ar_coefs  <- numeric(N)
  for (i in seq_len(N)) {
    dy_i <- diff(Y[, i])
    n_d  <- length(dy_i)
    if (n_d < 3L) { ar_resids[[i]] <- dy_i; ar_coefs[i] <- 0; next }
    y_lag <- dy_i[-n_d]
    y_cur <- dy_i[-1L]
    phi <- tryCatch(
      as.numeric(crossprod(y_lag, y_cur) / crossprod(y_lag)),
      error = function(e) 0
    )
    phi <- max(-0.99, min(0.99, phi))
    ar_coefs[i]  <- phi
    ar_resids[[i]] <- y_cur - phi * y_lag
  }

  for (b in seq_len(bootreps)) {
    Y_boot <- matrix(0, nrow = TT, ncol = N)
    for (i in seq_len(N)) {
      e_i <- ar_resids[[i]]
      n_e <- length(e_i)
      idx_b <- sample.int(n_e, TT - 1L, replace = TRUE)
      e_b <- e_i[idx_b]
      dy_boot <- numeric(TT - 1L)
      dy_boot[1L] <- e_b[1L]
      for (t in seq(2L, TT - 1L)) {
        dy_boot[t] <- ar_coefs[i] * dy_boot[t - 1L] + e_b[t]
      }
      Y_boot[, i] <- cumsum(c(0, dy_boot))
    }
    ts_b <- numeric(N)
    for (i in seq_len(N)) {
      r_b <- .xtpqroot_tfr_unit(Y_boot[, i], TT, model, maxlag)
      ts_b[i] <- if (is.na(r_b$tstat)) 0 else r_b$tstat
    }
    boot_stats[b] <- mean(ts_b)
  }

  pvalue <- mean(boot_stats <= tfr_panel)
  cv01   <- stats::quantile(boot_stats, 0.01)
  cv05   <- stats::quantile(boot_stats, 0.05)
  cv10   <- stats::quantile(boot_stats, 0.10)

  ind_df <- data.frame(
    panel      = seq_len(N),
    tstat      = tfr_stats,
    kfr        = sapply(ind_res, "[[", "kfr"),
    gamma      = sapply(ind_res, "[[", "gamma"),
    tau        = sapply(ind_res, "[[", "tau"),
    p_opt      = sapply(ind_res, "[[", "p_opt"),
    break_frac = sapply(ind_res, "[[", "break_date_frac"),
    stringsAsFactors = FALSE
  )

  list(
    tfr = tfr_panel, pvalue = pvalue,
    cv01 = as.numeric(cv01), cv05 = as.numeric(cv05), cv10 = as.numeric(cv10),
    ind_results = ind_df, bootreps = bootreps
  )
}


# ============================================================
# S3 Methods
# ============================================================

#' Print Method for xtpqroot Objects
#'
#' @param x An object of class \code{"xtpqroot"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtpqroot <- function(x, ...) {
  message(strrep("-", 65))
  if (x$test == "cipstau") {
    message("CIPS(tau) Panel Quantile Unit Root Test")
    message("Yang, Wei & Cai (2022, Economics Letters)")
  } else {
    message("tFR Panel Unit Root Test (Fourier + LST)")
    message("Corakci & Omay (2023, Renewable Energy)")
  }
  message(sprintf("Variable: %s | N = %d, T = %d | Model: %s",
    x$var, x$N, x$TT, x$model))
  message(strrep("-", 65))

  if (x$test == "cipstau") {
    cips_stars <- if (is.na(x$cips_pv)) "" else
                  if (x$cips_pv < 0.01) "***" else
                  if (x$cips_pv < 0.05) "**" else
                  if (x$cips_pv < 0.10) "*" else ""
    message(sprintf("CIPS (OLS):   %.4f  p-value: %s%s",
      x$cips,
      if (is.na(x$cips_pv)) "NA" else sprintf("%.4f", x$cips_pv),
      cips_stars))
    message("")
    message(sprintf("%-8s %-12s %-12s %s",
      "Tau", "CIPS(tau)", "p-value", "Decision"))
    message(strrep("-", 50))
    for (q_idx in seq_along(x$quantiles)) {
      pv_q <- x$cipstau_pv[q_idx]
      stars <- if (is.na(pv_q)) "" else
               if (pv_q < 0.01) "***" else
               if (pv_q < 0.05) "**" else
               if (pv_q < 0.10) "*" else ""
      dec <- if (is.na(pv_q)) "NA" else if (pv_q < 0.05) "Reject H0" else "Fail to reject"
      message(sprintf("%-8.2f %-12.4f %-12.4f%s %s",
        x$quantiles[q_idx], x$cipstau[q_idx], x$cipstau_pv[q_idx],
        stars, dec))
    }
  } else {
    stars <- if (is.na(x$pvalue)) "" else
             if (x$pvalue < 0.01) "***" else
             if (x$pvalue < 0.05) "**" else
             if (x$pvalue < 0.10) "*" else ""
    message(sprintf("tFR statistic: %.4f  Bootstrap p-value: %.4f%s",
      x$tfr, x$pvalue, stars))
    message(sprintf("Bootstrap CVs: 1%%=%.4f  5%%=%.4f  10%%=%.4f",
      x$cv01, x$cv05, x$cv10))
    message("")
    message(sprintf("%-8s %-10s %-8s %-8s %-8s %-8s",
      "Panel", "t_i,FR", "k^FR", "gamma", "tau", "ADF p"))
    for (i in seq_len(nrow(x$ind_results))) {
      r <- x$ind_results[i, ]
      message(sprintf("%-8d %-10.4f %-8.3f %-8.3f %-8.3f %-8d",
        r$panel, r$tstat, r$kfr, r$gamma, r$tau, r$p_opt))
    }
  }

  message(strrep("-", 65))
  message("H0: All panels contain a unit root")
  message("H1: Some panels are stationary")
  message("*** p<0.01, ** p<0.05, * p<0.10")
  invisible(x)
}

#' Summary Method for xtpqroot Objects
#'
#' @param object An object of class \code{"xtpqroot"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtpqroot <- function(object, ...) {
  print(object)
  invisible(object)
}

#' Example Panel Data for xtpqroot
#'
#' Returns a small balanced panel dataset for use in examples and testing.
#'
#' @return A data frame with columns \code{firm}, \code{year}, \code{invest},
#'   and \code{mvalue}.
#'
#' @examples
#' dat <- grunfeld_pqroot()
#' head(dat)
#'
#' @export
grunfeld_pqroot <- function() {
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

