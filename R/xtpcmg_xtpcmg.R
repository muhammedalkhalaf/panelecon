#' Panel Cointegrating Polynomial Regressions via FM-OLS
#'
#' Estimates a polynomial cointegrating regression in a panel using either
#' Group-Mean FM-OLS (Wagner & Reichold 2023) or Pooled FM-OLS (de Jong &
#' Wagner 2022). Models the long-run relationship:
#'
#' \deqn{y_{it} = \alpha_i + \beta_1 x_{it} + \beta_2 x_{it}^2 [+ \beta_3 x_{it}^3]
#'        [+ \gamma z_{it}] + u_{it}}
#'
#' where \eqn{x_{it}} and \eqn{z_{it}} are I(1) processes.
#'
#' @param data A data frame in long format.
#' @param y Character. Name of the dependent variable.
#' @param x Character. Name of the polynomial (I(1)) regressor.
#' @param panel_id Character. Name of the panel identifier variable.
#' @param time_id Character. Name of the time variable.
#' @param model Character. Estimator: \code{"mg"} for Group-Mean FM-OLS
#'   (default) or \code{"pmg"} for Pooled FM-OLS.
#' @param q Integer. Polynomial degree: \code{2} (quadratic, default) or
#'   \code{3} (cubic).
#' @param controls Character vector. Names of additional I(1) control
#'   variables. Default is \code{NULL} (no controls).
#' @param trend Integer. Deterministic trend type: \code{1} for demeaning
#'   only (default), \code{2} for demeaning and detrending.
#' @param kernel Character. HAC kernel: \code{"ba"} (Bartlett, default),
#'   \code{"pa"} (Parzen), \code{"qs"} (Quadratic Spectral),
#'   \code{"tr"} (truncated), \code{"bo"} (Bohman).
#' @param bw Character or numeric. Bandwidth for HAC estimation. \code{"And91"}
#'   (default) uses Andrews (1991) automatic selection. A numeric value sets
#'   the bandwidth directly.
#' @param effects Character. For Pooled FM-OLS: \code{"oneway"} (default)
#'   for one-way demeaning or \code{"twoways"} for two-way demeaning.
#' @param corr_rob Logical. For Group-Mean FM-OLS: if \code{TRUE}, use
#'   cross-sectionally robust VCV. Default is \code{FALSE}.
#'
#' @return An object of class \code{"xtpcmg"} containing:
#' \describe{
#'   \item{coefficients}{Named numeric vector of FM-OLS coefficient estimates.}
#'   \item{vcov}{Variance-covariance matrix.}
#'   \item{se}{Standard errors.}
#'   \item{tstat}{t-statistics.}
#'   \item{pvalue}{p-values (two-sided, normal approximation).}
#'   \item{model}{Character. \code{"mg"} or \code{"pmg"}.}
#'   \item{q}{Integer. Polynomial degree.}
#'   \item{N}{Integer. Number of panel units.}
#'   \item{TT}{Integer. Number of time periods.}
#'   \item{y}{Character. Dependent variable name.}
#'   \item{x}{Character. Polynomial variable name.}
#'   \item{tp}{Numeric. Turning point x* (q=2 only; \code{NA} if not applicable).}
#'   \item{tp_se}{Numeric. Delta-method SE for turning point.}
#'   \item{tp_lo}{Numeric. Lower 95% CI bound for turning point.}
#'   \item{tp_hi}{Numeric. Upper 95% CI bound for turning point.}
#'   \item{ind_coef}{Matrix of individual FM-OLS estimates (MG model, N x q).}
#'   \item{swamy_s}{Numeric. Swamy S-statistic for slope homogeneity.}
#'   \item{swamy_p}{Numeric. Swamy test p-value.}
#' }
#'
#' @references
#' Andrews, D.W.K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817--858.
#' \doi{10.2307/2938229}
#'
#' de Jong, R.M. and Wagner, M. (2022). Panel cointegrating polynomial
#' regressions. \emph{Annals of Applied Statistics}, 16(1), 416--442.
#' \doi{10.1214/21-AOAS1536}
#'
#' Wagner, M. and Reichold, K. (2023). Panel cointegrating polynomial
#' regressions. \emph{Econometric Reviews}, 42(9--10), 782--827.
#' \doi{10.1080/07474938.2022.2070522}
#'
#' @examples
#' dat <- grunfeld_cmg()
#'
#' # Group-Mean FM-OLS (quadratic EKC-type model)
#' res <- xtpcmg(dat, y = "invest", x = "mvalue",
#'               panel_id = "firm", time_id = "year",
#'               model = "mg", q = 2L)
#' print(res)
#' summary(res)
#'
#' # Pooled FM-OLS
#' res2 <- xtpcmg(dat, y = "invest", x = "mvalue",
#'                panel_id = "firm", time_id = "year",
#'                model = "pmg", q = 2L)
#' print(res2)
#'
#' @export
xtpcmg <- function(data, y, x, panel_id, time_id,
                   model = c("mg", "pmg"),
                   q = 2L,
                   controls = NULL,
                   trend = 1L,
                   kernel = "ba",
                   bw = "And91",
                   effects = "oneway",
                   corr_rob = FALSE) {

  model <- match.arg(model)
  q     <- as.integer(q)
  if (!q %in% c(2L, 3L)) stop("'q' must be 2 or 3.")
  if (!kernel %in% c("ba", "pa", "qs", "tr", "bo", "da")) {
    stop("'kernel' must be one of: 'ba', 'pa', 'qs', 'tr', 'bo', 'da'.")
  }
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  for (v in c(y, x, panel_id, time_id)) {
    if (!v %in% names(data)) stop(sprintf("Variable '%s' not found.", v))
  }
  if (!is.null(controls)) {
    for (cv in controls) {
      if (!cv %in% names(data)) stop(sprintf("Control variable '%s' not found.", cv))
    }
  }

  data <- data[order(data[[panel_id]], data[[time_id]]), ]
  N    <- length(unique(data[[panel_id]]))
  TT   <- length(unique(data[[time_id]]))
  if (nrow(data) != N * TT) stop("Panel must be strongly balanced.")
  if (N < 2L) stop("At least 2 panel units are required.")
  if (TT < 5L) stop("At least 5 time periods are required.")

  Y_mat <- matrix(data[[y]], nrow = TT, ncol = N)
  X_mat <- matrix(data[[x]], nrow = TT, ncol = N)
  Z_list <- if (!is.null(controls)) {
    lapply(controls, function(cv) matrix(data[[cv]], nrow = TT, ncol = N))
  } else list()

  nc <- length(Z_list)

  if (model == "mg") {
    res <- .xtpcmg_mg(Y_mat, X_mat, Z_list, N, TT, q, trend, kernel, bw, corr_rob)
  } else {
    res <- .xtpcmg_pmg(Y_mat, X_mat, Z_list, N, TT, q, kernel, bw, effects)
  }

  coef_names <- c(x, paste0(x, "^2"))
  if (q == 3L) coef_names <- c(coef_names, paste0(x, "^3"))
  if (nc > 0L) coef_names <- c(coef_names, controls)
  names(res$coefficients) <- coef_names
  colnames(res$vcov)      <- coef_names
  rownames(res$vcov)      <- coef_names

  se_vec <- sqrt(diag(res$vcov))
  tstat  <- res$coefficients / se_vec
  pval   <- 2 * (1 - stats::pnorm(abs(tstat)))

  # Turning point (quadratic)
  tp <- tp_se <- tp_lo <- tp_hi <- tp_z <- tp_p <- NA_real_
  if (q == 2L) {
    b1 <- res$coefficients[1L]
    b2 <- res$coefficients[2L]
    if (abs(b2) > 1e-12) {
      tp    <- -b1 / (2 * b2)
      g1    <- -1 / (2 * b2)
      g2    <- b1 / (2 * b2^2)
      v11   <- res$vcov[1L, 1L]
      v12   <- res$vcov[1L, 2L]
      v22   <- res$vcov[2L, 2L]
      tp_var <- g1^2 * v11 + 2 * g1 * g2 * v12 + g2^2 * v22
      tp_se  <- sqrt(abs(tp_var))
      tp_lo  <- tp - stats::qnorm(0.975) * tp_se
      tp_hi  <- tp + stats::qnorm(0.975) * tp_se
      tp_z   <- if (tp_se > 0) tp / tp_se else NA_real_
      tp_p   <- if (!is.na(tp_z)) 2 * (1 - stats::pnorm(abs(tp_z))) else NA_real_
    }
  }

  # Swamy test (MG model)
  swamy_s <- swamy_df <- swamy_p <- NA_real_
  if (model == "mg" && !is.null(res$ind_coef)) {
    sw <- .xtpcmg_swamy(res$ind_coef, N, q + nc)
    swamy_s  <- sw$stat
    swamy_df <- sw$df
    swamy_p  <- sw$pvalue
  }

  structure(
    list(
      coefficients = res$coefficients,
      vcov         = res$vcov,
      se           = se_vec,
      tstat        = tstat,
      pvalue       = pval,
      model        = model,
      q            = q,
      N            = N,
      TT           = TT,
      y            = y,
      x            = x,
      controls     = controls,
      kernel       = kernel,
      bw_type      = bw,
      tp           = tp,
      tp_se        = tp_se,
      tp_lo        = tp_lo,
      tp_hi        = tp_hi,
      tp_z         = tp_z,
      tp_p         = tp_p,
      ind_coef     = res$ind_coef,
      swamy_s      = swamy_s,
      swamy_df     = swamy_df,
      swamy_p      = swamy_p
    ),
    class = "xtpcmg"
  )
}


# ============================================================
# HAC Kernel Weights
# ============================================================

#' @keywords internal
.xtpcmg_lr_weights <- function(TT, kern, band) {
  w <- numeric(TT - 1L)
  M <- band
  j_max <- 0L

  if (kern == "tr") {
    j_max <- min(floor(M), TT - 1L)
    if (j_max >= 1L) w[seq_len(j_max)] <- 1
  } else if (kern == "ba") {
    j_max <- min(ceiling(M) - 1L, TT - 1L)
    if (j_max >= 1L) for (j in seq_len(j_max)) w[j] <- 1 - j / M
  } else if (kern == "pa") {
    j_max <- min(ceiling(M) - 1L, TT - 1L)
    ulim1 <- min(floor(M / 2), TT - 1L)
    for (j in seq_len(j_max)) {
      jj <- j / M
      if (j <= ulim1) {
        w[j] <- 1 - 6 * jj^2 + 6 * jj^3
      } else {
        w[j] <- 2 * (1 - jj)^3
      }
    }
  } else if (kern == "bo") {
    j_max <- min(ceiling(M) - 1L, TT - 1L)
    if (j_max >= 1L) for (j in seq_len(j_max)) {
      jj <- j / M
      w[j] <- (1 - jj) * cos(pi * jj) + sin(pi * jj) / pi
    }
  } else if (kern == "qs") {
    j_max <- TT - 1L
    sc <- 6 * pi / 5
    for (j in seq_len(j_max)) {
      jj <- j / M
      w[j] <- 25 / (12 * pi^2 * jj^2) * (sin(sc * jj) / (sc * jj) - cos(sc * jj))
    }
  } else if (kern == "da") {
    j_max <- TT - 1L
    for (j in seq_len(j_max)) {
      jj <- j / M
      w[j] <- sin(pi * jj) / (pi * jj)
    }
  }

  list(w = w, j_max = j_max)
}

#' @keywords internal
.xtpcmg_lr_var <- function(u, kern, band, demean = FALSE) {
  TT <- nrow(u)
  m  <- ncol(u)
  if (demean) u <- u - matrix(colMeans(u), nrow = TT, ncol = m, byrow = TRUE)

  wl <- .xtpcmg_lr_weights(TT, kern, band)
  w  <- wl$w
  j_max <- wl$j_max

  Sigma <- crossprod(u) / TT
  Omega <- Sigma
  Delta <- Sigma

  if (j_max >= 1L) {
    for (j in seq_len(j_max)) {
      T1 <- crossprod(u[seq(j + 1L, TT), , drop = FALSE],
                      u[seq(1L, TT - j), , drop = FALSE]) / TT
      T2 <- t(T1)
      Omega <- Omega + w[j] * (T1 + T2)
      Delta <- Delta + w[j] * T1
    }
  }
  list(Omega = Omega, Delta = Delta, Sigma = Sigma)
}

#' @keywords internal
.xtpcmg_and91 <- function(v, kern) {
  TT   <- nrow(v)
  dimv <- ncol(v)
  rho  <- numeric(dimv)
  sig2 <- numeric(dimv)

  for (k in seq_len(dimv)) {
    yy <- v[seq(2L, TT), k]
    xx <- v[seq(1L, TT - 1L), k]
    denom_k <- sum(xx^2)
    r <- if (denom_k > 0) sum(xx * yy) / denom_k else 0
    rho[k]  <- r
    e_k <- yy - xx * r
    sig2[k] <- mean(e_k^2)
  }

  denom <- sum(sig2^2 / (1 - rho)^4)
  if (denom < 1e-15) denom <- 1e-15

  if (kern == "ba") {
    a1 <- sum(4 * rho^2 * sig2^2 / ((1 - rho)^6 * (1 + rho)^2)) / denom
    bw <- ceiling(1.1447 * (a1 * TT)^(1 / 3))
  } else {
    a2 <- sum(4 * rho^2 * sig2^2 / (1 - rho)^8) / denom
    bw <- if (kern == "tr") ceiling(0.6611 * (a2 * TT)^(1 / 5)) else
           if (kern == "pa") ceiling(2.6614 * (a2 * TT)^(1 / 5)) else
           ceiling(1.3221 * (a2 * TT)^(1 / 5))
  }
  max(bw, 1L)
}

#' @keywords internal
.xtpcmg_get_bw <- function(u, v, kern, bw_spec) {
  if (is.numeric(bw_spec) && is.finite(bw_spec)) return(as.numeric(bw_spec))
  .xtpcmg_and91(cbind(u, v), kern)
}

#' @keywords internal
.xtpcmg_demean <- function(y, x, typee) {
  TT <- nrow(x)
  if (typee == 1L) {
    D <- matrix(1, nrow = TT, ncol = 1L)
  } else {
    D <- cbind(1, seq_len(TT))
  }
  P <- diag(TT) - D %*% solve(crossprod(D), t(D))
  list(yt = P %*% y, xt = P %*% x,
       x2t = P %*% x^2,
       x3t = P %*% x^3)
}


# ============================================================
# Group-Mean FM-OLS
# ============================================================

#' @keywords internal
.xtpcmg_mg <- function(Y_mat, X_mat, Z_list, N, TT, q, typee, kern, bw_spec,
                        corr_rob) {
  nc <- length(Z_list)
  p  <- q + nc

  # First difference of x
  dX <- rbind(matrix(0, 1L, N), diff(X_mat))

  betaGM   <- numeric(p)
  V_sum    <- matrix(0, p, p)
  ind_coef <- matrix(0, N, p)
  u_all    <- matrix(0, TT - 1L, N)
  v_all    <- matrix(0, TT - 1L, N)

  # Detrend per panel
  yt_list  <- vector("list", N)
  xt_list  <- vector("list", N)
  x2t_list <- vector("list", N)
  x3t_list <- vector("list", N)
  zt_list  <- vector("list", N * nc)
  for (i in seq_len(N)) {
    dm <- .xtpcmg_demean(Y_mat[-1L, i, drop = FALSE],
                          X_mat[-1L, i, drop = FALSE], typee)
    yt_list[[i]]  <- dm$yt
    xt_list[[i]]  <- dm$xt
    x2t_list[[i]] <- dm$x2t
    x3t_list[[i]] <- dm$x3t
    if (nc > 0L) {
      for (k in seq_len(nc)) {
        zm <- .xtpcmg_demean(Z_list[[k]][-1L, i, drop = FALSE],
                              X_mat[-1L, i, drop = FALSE], typee)
        zt_list[[(k - 1L) * N + i]] <- zm$yt
      }
    }
  }

  TT_fd <- TT - 1L

  # First pass: OLS to get residuals and bandwidths
  for (i in seq_len(N)) {
    if (q == 2L) {
      Xi <- cbind(xt_list[[i]], x2t_list[[i]])
    } else {
      Xi <- cbind(xt_list[[i]], x2t_list[[i]], x3t_list[[i]])
    }
    if (nc > 0L) {
      for (k in seq_len(nc)) Xi <- cbind(Xi, zt_list[[(k - 1L) * N + i]])
    }
    yi  <- yt_list[[i]]
    ols_b <- tryCatch(solve(crossprod(Xi), crossprod(Xi, yi)),
                      error = function(e) rep(0, p))
    u_ols <- as.numeric(yi - Xi %*% ols_b)
    v_i   <- dX[-1L, i] - mean(dX[-1L, i])
    u_all[, i] <- u_ols
    v_all[, i] <- v_i
  }

  # FM-OLS per panel
  for (i in seq_len(N)) {
    u_i <- u_all[, i]
    v_i <- v_all[, i]
    xi_i <- X_mat[-1L, i]
    yi_i <- yt_list[[i]]

    bw_i <- .xtpcmg_get_bw(u_i, v_i, kern, bw_spec)
    lr_i <- .xtpcmg_lr_var(cbind(u_i, v_i), kern, bw_i)

    Omega_uu <- lr_i$Omega[1L, 1L]
    Omega_vv <- lr_i$Omega[2L, 2L]
    Omega_uv <- lr_i$Omega[1L, 2L]
    Delta_vv <- lr_i$Delta[2L, 2L]
    Delta_vu <- lr_i$Delta[2L, 1L]
    Delta_vu_plus <- Delta_vu - Delta_vv * Omega_uv / Omega_vv
    Omega_udotv   <- Omega_uu - Omega_uv^2 / Omega_vv

    if (q == 2L) {
      Xi <- cbind(xt_list[[i]], x2t_list[[i]])
    } else {
      Xi <- cbind(xt_list[[i]], x2t_list[[i]], x3t_list[[i]])
    }
    if (nc > 0L) {
      for (k in seq_len(nc)) Xi <- cbind(Xi, zt_list[[(k - 1L) * N + i]])
    }

    # Bias-corrected y
    yi_plus <- yi_i - v_i * Omega_uv / Omega_vv

    # Bias correction vector (polynomial part)
    if (q == 2L) {
      C_i <- Delta_vu_plus * c(TT_fd, 2 * sum(xi_i))
    } else {
      C_i <- Delta_vu_plus * c(TT_fd, 2 * sum(xi_i), 3 * sum(xi_i^2))
    }
    if (nc > 0L) C_i <- c(C_i, rep(0, nc))

    XtX_inv <- tryCatch(solve(crossprod(Xi)), error = function(e) NULL)
    if (is.null(XtX_inv)) next

    beta_fm <- as.numeric(XtX_inv %*% (crossprod(Xi, yi_plus) - C_i))
    ind_coef[i, ] <- beta_fm
    betaGM <- betaGM + beta_fm

    if (!corr_rob) {
      V_sum <- V_sum + Omega_udotv * XtX_inv
    }
  }

  betaGM <- betaGM / N

  if (!corr_rob) {
    V_out <- V_sum / N^2
  } else {
    # Cross-sectional robust VCV
    bw_all <- .xtpcmg_get_bw(u_all, v_all, kern, bw_spec)
    lr_all <- .xtpcmg_lr_var(cbind(u_all, v_all), kern, bw_all)
    V_sum2 <- matrix(0, p, p)
    for (i in seq_len(N)) {
      if (q == 2L) Xi <- cbind(xt_list[[i]], x2t_list[[i]])
      else Xi <- cbind(xt_list[[i]], x2t_list[[i]], x3t_list[[i]])
      if (nc > 0L) for (k in seq_len(nc)) Xi <- cbind(Xi, zt_list[[(k - 1L) * N + i]])
      Mii <- crossprod(Xi)
      for (j in seq_len(N)) {
        if (q == 2L) Xj <- cbind(xt_list[[j]], x2t_list[[j]])
        else Xj <- cbind(xt_list[[j]], x2t_list[[j]], x3t_list[[j]])
        if (nc > 0L) for (k in seq_len(nc)) Xj <- cbind(Xj, zt_list[[(k - 1L) * N + j]])
        Mij <- crossprod(Xi, Xj)
        Mjj <- crossprod(Xj)
        Ouu_ij <- lr_all$Omega[i, j]
        Mii_inv <- tryCatch(solve(Mii), error = function(e) NULL)
        Mjj_inv <- tryCatch(solve(Mjj), error = function(e) NULL)
        if (!is.null(Mii_inv) && !is.null(Mjj_inv)) {
          V_sum2 <- V_sum2 + Ouu_ij * (Mii_inv %*% Mij %*% Mjj_inv)
        }
      }
    }
    V_out <- V_sum2 / N^2
  }

  list(coefficients = betaGM, vcov = V_out, ind_coef = ind_coef)
}


# ============================================================
# Pooled FM-OLS
# ============================================================

#' @keywords internal
.xtpcmg_simple_demean <- function(y_mat, x_mat, way, q, Z_list, nc) {
  TT <- nrow(x_mat)
  N  <- ncol(x_mat)

  if (way == "oneway") {
    yt  <- y_mat - matrix(colMeans(y_mat), TT, N, byrow = TRUE)
    xt  <- x_mat - matrix(colMeans(x_mat), TT, N, byrow = TRUE)
    x2t <- x_mat^2 - matrix(colMeans(x_mat^2), TT, N, byrow = TRUE)
    x3t <- if (q == 3L) x_mat^3 - matrix(colMeans(x_mat^3), TT, N, byrow = TRUE) else NULL
    zt_list <- lapply(Z_list, function(Z) Z - matrix(colMeans(Z), TT, N, byrow = TRUE))
  } else {
    # Two-way: subtract col means, row means, add grand mean
    yt  <- .xtpcmg_twoway(y_mat)
    xt  <- .xtpcmg_twoway(x_mat)
    x2t <- .xtpcmg_twoway(x_mat^2)
    x3t <- if (q == 3L) .xtpcmg_twoway(x_mat^3) else NULL
    zt_list <- lapply(Z_list, .xtpcmg_twoway)
  }
  list(yt = yt, xt = xt, x2t = x2t, x3t = x3t, zt_list = zt_list)
}

#' @keywords internal
.xtpcmg_twoway <- function(M) {
  TT <- nrow(M); N <- ncol(M)
  M - matrix(colMeans(M), TT, N, byrow = TRUE) -
    matrix(rowMeans(M), TT, N) + mean(M)
}

#' @keywords internal
.xtpcmg_pmg <- function(Y_mat, X_mat, Z_list, N, TT, q, kern, bw_spec, effects) {
  nc <- length(Z_list)
  p  <- q + nc

  dm <- .xtpcmg_simple_demean(Y_mat, X_mat, effects, q, Z_list, nc)
  yt_vec  <- as.vector(dm$yt)
  xt_vec  <- as.vector(dm$xt)
  x2t_vec <- as.vector(dm$x2t)

  Xvec <- if (q == 2L) cbind(xt_vec, x2t_vec) else
    cbind(xt_vec, x2t_vec, as.vector(dm$x3t))
  if (nc > 0L) {
    for (k in seq_len(nc)) Xvec <- cbind(Xvec, as.vector(dm$zt_list[[k]]))
  }

  # LSDV
  XtX_inv <- tryCatch(solve(crossprod(Xvec)), error = function(e) NULL)
  if (is.null(XtX_inv)) {
    return(list(coefficients = rep(0, p),
                vcov = diag(p), ind_coef = NULL))
  }
  beta_lsdv <- as.numeric(XtX_inv %*% crossprod(Xvec, yt_vec))
  u_vec  <- as.numeric(yt_vec - Xvec %*% beta_lsdv)

  # Long-run var per panel, pooled correction
  Sum_Lr <- matrix(0, 2L, 2L)
  Sum_Dr <- matrix(0, 2L, 2L)
  Sum_FM  <- numeric(p)

  for (i in seq_len(N)) {
    idx_i <- seq((i - 1L) * TT + 1L, i * TT)
    u_i   <- u_vec[idx_i]
    vt_i  <- c(X_mat[1L, i], diff(X_mat[, i]))
    bw_i  <- .xtpcmg_get_bw(u_i, vt_i, kern, bw_spec)
    lr_i  <- .xtpcmg_lr_var(cbind(u_i, vt_i), kern, bw_i)
    Sum_Lr <- Sum_Lr + lr_i$Omega
    Sum_Dr <- Sum_Dr + lr_i$Delta
  }

  Lr_mean <- Sum_Lr / N
  Dr_mean <- Sum_Dr / N

  Dr_vu_plus <- Dr_mean[2L, 1L] -
    Dr_mean[2L, 2L] * Lr_mean[2L, 1L] / Lr_mean[2L, 2L]

  for (i in seq_len(N)) {
    xi_i <- X_mat[, i]
    idx_i <- seq((i - 1L) * TT + 1L, i * TT)
    vt_i  <- c(X_mat[1L, i], diff(xi_i))

    if (q == 2L) {
      Mi <- c(TT, 2 * sum(xi_i))
    } else {
      Mi <- c(TT, 2 * sum(xi_i), 3 * sum(xi_i^2))
    }
    if (nc > 0L) Mi <- c(Mi, rep(0, nc))

    C_plus <- Dr_vu_plus * c(Mi[seq_len(q)], rep(0, nc))

    Xi_i  <- Xvec[idx_i, , drop = FALSE]
    yi_plus <- as.vector(dm$yt)[idx_i] -
      Lr_mean[1L, 2L] / Lr_mean[2L, 2L] * vt_i
    Sum_FM <- Sum_FM + as.numeric(crossprod(Xi_i, yi_plus) - C_plus)
  }

  beta_fm <- as.numeric(XtX_inv %*% Sum_FM)

  # VCV (simplified sandwich)
  u_fm  <- yt_vec - Xvec %*% beta_fm
  sigma2 <- sum(u_fm^2) / (N * TT - p)
  VCV    <- sigma2 * XtX_inv

  list(coefficients = beta_fm, vcov = VCV, ind_coef = NULL)
}


# ============================================================
# Swamy Test
# ============================================================

#' @keywords internal
.xtpcmg_swamy <- function(ind_coef, N, K) {
  b_bar <- colMeans(ind_coef)
  V_emp <- matrix(0, K, K)
  for (i in seq_len(N)) {
    d <- ind_coef[i, ] - b_bar
    V_emp <- V_emp + outer(d, d)
  }
  V_emp <- V_emp / (N - 1L)
  V_inv <- tryCatch(solve(V_emp), error = function(e) NULL)
  if (is.null(V_inv)) return(list(stat = NA_real_, df = NA_real_, pvalue = NA_real_))

  S_stat <- 0
  for (i in seq_len(N)) {
    d <- ind_coef[i, ] - b_bar
    S_stat <- S_stat + as.numeric(d %*% V_inv %*% d)
  }
  df_sw <- (N - 1L) * K
  list(stat = S_stat, df = df_sw,
       pvalue = 1 - stats::pchisq(S_stat, df = df_sw))
}


# ============================================================
# S3 Methods
# ============================================================

#' Print Method for xtpcmg Objects
#'
#' @param x An object of class \code{"xtpcmg"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtpcmg <- function(x, ...) {
  model_lab <- if (x$model == "mg") "Group-Mean FM-OLS (Wagner & Reichold 2023)" else
                 "Pooled FM-OLS (de Jong & Wagner 2022)"
  message(strrep("-", 65))
  message(sprintf("Panel Cointegrating Polynomial Regression: %s", model_lab))
  message(sprintf("Dep: %s | Poly: %s (degree %d) | N=%d, T=%d",
    x$y, x$x, x$q, x$N, x$TT))
  message(strrep("-", 65))
  message(sprintf("%-16s %10s %10s %10s %10s",
    "Variable", "Estimate", "Std. Err", "t-stat", "p-value"))
  message(strrep("-", 65))
  for (j in seq_along(x$coefficients)) {
    stars <- if (x$pvalue[j] < 0.01) "***" else if (x$pvalue[j] < 0.05) "**" else
              if (x$pvalue[j] < 0.10) "*" else ""
    message(sprintf("%-16s %10.4f %10.4f %10.4f %10.4f%s",
      names(x$coefficients)[j],
      x$coefficients[j], x$se[j], x$tstat[j], x$pvalue[j], stars))
  }
  message(strrep("-", 65))
  if (!is.na(x$tp)) {
    tp_shape <- if (x$coefficients[2L] < 0) "(Inverted-U)" else "(U-shaped)"
    message(sprintf("Turning point x*: %.4f [95%% CI: %.4f, %.4f] %s",
      x$tp, x$tp_lo, x$tp_hi, tp_shape))
  }
  if (!is.na(x$swamy_s) && x$model == "mg") {
    message(sprintf("Swamy S = %.4f, df = %d, p = %.4f %s",
      x$swamy_s, as.integer(x$swamy_df), x$swamy_p,
      if (!is.na(x$swamy_p) && x$swamy_p < 0.05) "[Heterogeneous]" else "[Homogeneous]"))
  }
  message("*** p<0.01, ** p<0.05, * p<0.10")
  invisible(x)
}

#' Summary Method for xtpcmg Objects
#'
#' @param object An object of class \code{"xtpcmg"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtpcmg <- function(object, ...) {
  print(object)
  if (object$model == "mg" && !is.null(object$ind_coef)) {
    message("")
    message("Individual FM-OLS Estimates (first 5 panels):")
    n_show <- min(5L, nrow(object$ind_coef))
    for (i in seq_len(n_show)) {
      message(sprintf("  Panel %d: %s",
        i,
        paste(sprintf("%.4f", object$ind_coef[i, ]), collapse = "  ")))
    }
  }
  invisible(object)
}

#' Example Panel Data for xtpcmg
#'
#' Returns the Grunfeld (1958) balanced panel dataset for examples.
#'
#' @return A data frame with columns \code{firm}, \code{year},
#'   \code{invest}, and \code{mvalue}.
#'
#' @examples
#' dat <- grunfeld_cmg()
#' head(dat)
#'
#' @export
grunfeld_cmg <- function() {
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
