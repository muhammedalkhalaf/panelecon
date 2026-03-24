#' Fixed Effects Filtered and Vector Decomposition for Panel Data
#'
#' Estimates panel data models with time-invariant regressors using the Fixed
#' Effects Vector Decomposition (FEVD), Fixed Effects Filtered (FEF), or
#' FEF-IV methods. Provides consistent standard errors following Pesaran and
#' Zhou (2016).
#'
#' @param formula A formula of the form \code{y ~ x1 + x2} specifying the
#'   dependent variable and time-varying regressors.
#' @param data A data frame containing the panel data in long format.
#' @param index A character vector of length 2: \code{c("id_var", "time_var")}.
#' @param zinvariants A character vector naming the time-invariant regressors
#'   (variables that do not vary within panel units).
#' @param method Estimation method: \code{"fevd"} (default), \code{"fef"}, or
#'   \code{"fefiv"}.
#' @param instruments A character vector of instrument variable names. Required
#'   when \code{method = "fefiv"}.
#' @param robust Logical. If \code{TRUE}, reports heteroskedasticity-robust
#'   standard errors for the time-varying part. Default is \code{FALSE}.
#'
#' @return An object of class \code{"xtfifevd"} with components:
#'   \describe{
#'     \item{beta_fe}{Numeric vector of FE (within) estimates for time-varying regressors.}
#'     \item{gamma}{Numeric vector of estimates for time-invariant regressors.}
#'     \item{alpha}{Intercept estimate.}
#'     \item{se_beta}{Standard errors for beta_fe.}
#'     \item{se_gamma}{Pesaran-Zhou corrected standard errors for gamma.}
#'     \item{se_alpha}{Standard error for alpha.}
#'     \item{V_gamma_pz}{Pesaran-Zhou variance matrix for gamma (k_z x k_z).}
#'     \item{V_gamma_fevd}{Raw FEVD variance for gamma (only for FEVD method).}
#'     \item{delta}{Stage-3 coefficient on the FEVD error component (should be 1).}
#'     \item{sigma2_e}{Within-unit error variance from FE stage.}
#'     \item{sigma2_u}{Between-unit variance from FE stage.}
#'     \item{N}{Number of observations.}
#'     \item{N_g}{Number of panel units.}
#'     \item{T_avg}{Average time periods per unit.}
#'     \item{method}{Method used: "FEVD", "FEF", or "FEF-IV".}
#'     \item{depvar}{Name of the dependent variable.}
#'     \item{xvars}{Names of time-varying regressors.}
#'     \item{zinvariants}{Names of time-invariant regressors.}
#'   }
#'
#' @details
#' \strong{FEVD} (Plumper & Troeger 2007): Three-stage estimator.
#' Stage 1 fits a standard within (FE) estimator. Stage 2 decomposes the
#' unit fixed effects into a part explained by time-invariant regressors plus
#' a residual. Stage 3 pools all variables. Point estimates are identical to
#' FEF but the raw Stage-3 SEs are inconsistent (Pesaran & Zhou 2016, Remark 4);
#' this package uses the corrected Pesaran-Zhou SEs.
#'
#' \strong{FEF} (Pesaran & Zhou 2016): Two-stage estimator equivalent to FEVD
#' in terms of point estimates but using theoretically justified standard errors.
#'
#' \strong{FEF-IV}: Instrumental variables version of FEF for endogenous
#' time-invariant regressors. Requires \code{instruments} with at least as many
#' variables as \code{zinvariants}.
#'
#' @references
#' Plumper, T. and Troeger, V.E. (2007). Efficient Estimation of Time-Invariant
#' and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
#' Effects. \emph{Political Analysis}, 15(2), 124--139.
#' \doi{10.1093/pan/mpm002}
#'
#' Pesaran, M.H. and Zhou, Q. (2016). Estimation of Time-Invariant Effects in
#' Static Panel Data Models. \emph{Econometric Reviews}, 37(10), 1137--1171.
#' \doi{10.1080/07474938.2015.1032164}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 10; tt <- 15
#' uid  <- rep(1:n, each = tt)
#' tval <- rep(1:tt, times = n)
#' z_i  <- rep(rnorm(n), each = tt)      # time-invariant
#' x_it <- rnorm(n * tt)                  # time-varying
#' y    <- 0.5 * x_it + 0.8 * z_i + rnorm(n * tt, sd = 0.5)
#' dat  <- data.frame(id = uid, time = tval, y = y, x = x_it, z = z_i)
#'
#' res <- xtfifevd(y ~ x, data = dat, index = c("id", "time"),
#'                 zinvariants = "z", method = "fevd")
#' print(res)
#' summary(res)
#' }
#'
#' @export
xtfifevd <- function(formula, data, index, zinvariants,
                     method = c("fevd", "fef", "fefiv"),
                     instruments = NULL,
                     robust = FALSE) {

  method <- match.arg(method)

  ## ── Input validation ────────────────────────────────────────────────────
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.character(index) || length(index) != 2) {
    stop("'index' must be a character vector of length 2.", call. = FALSE)
  }
  if (!all(index %in% names(data))) {
    stop("Variables in 'index' not found in 'data'.", call. = FALSE)
  }
  if (!is.character(zinvariants) || length(zinvariants) < 1) {
    stop("'zinvariants' must be a character vector of at least one variable name.",
         call. = FALSE)
  }
  if (!all(zinvariants %in% names(data))) {
    stop("Some variables in 'zinvariants' not found in 'data'.", call. = FALSE)
  }
  if (method == "fefiv") {
    if (is.null(instruments)) {
      stop("'instruments' must be provided for method = 'fefiv'.", call. = FALSE)
    }
    if (length(instruments) < length(zinvariants)) {
      stop("Number of instruments must be >= number of 'zinvariants'.", call. = FALSE)
    }
    if (!all(instruments %in% names(data))) {
      stop("Some variables in 'instruments' not found in 'data'.", call. = FALSE)
    }
  }

  ivar <- index[1]
  tvar <- index[2]

  ## ── Prepare data ─────────────────────────────────────────────────────────
  mf      <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  depvar  <- names(mf)[1]
  xvars   <- names(mf)[-1]
  k_x     <- length(xvars)
  k_z     <- length(zinvariants)

  if (k_x < 1) stop("At least one time-varying regressor required.", call. = FALSE)

  # Complete cases across all variables
  all_vars <- c(depvar, xvars, zinvariants, ivar, tvar)
  if (!is.null(instruments)) all_vars <- c(all_vars, instruments)
  keep <- stats::complete.cases(data[, intersect(all_vars, names(data)), drop = FALSE])
  data_c <- data[keep, , drop = FALSE]

  panels  <- sort(unique(data_c[[ivar]]))
  N_g     <- length(panels)
  N_obs   <- nrow(data_c)

  if (N_g < 2) stop("At least 2 panel units are required.", call. = FALSE)

  T_vals <- as.integer(table(data_c[[ivar]]))
  T_avg  <- mean(T_vals)

  ## ── Extract matrices ─────────────────────────────────────────────────────
  Y   <- data_c[[depvar]]
  X   <- as.matrix(data_c[, xvars, drop = FALSE])
  Z   <- as.matrix(data_c[, zinvariants, drop = FALSE])
  uid <- data_c[[ivar]]

  ## ── Stage 1: FE (within) estimation ─────────────────────────────────────
  ## Within-demean Y and X
  Yw  <- .within_demean(Y, uid)
  Xw  <- .within_demean_mat(X, uid)
  fe  <- .ols_fit(Yw, Xw)    # coefficients (no intercept needed after demeaning)

  beta_fe <- as.numeric(fe$coef)

  ## FE residuals: u_hat = Y - X*beta_fe (includes unit effects)
  u_hat     <- Y - X %*% beta_fe
  u_bar_i   <- tapply(u_hat, uid, mean)    # unit-level means

  ## Within residuals for sigma2_e
  e_hat     <- Yw - Xw %*% beta_fe
  sigma2_e  <- sum(e_hat^2) / max(N_obs - N_g - k_x, 1)

  ## Between variance
  sigma2_u  <- max(0, var(u_bar_i) - sigma2_e / T_avg)

  ## SE for beta_fe (HC sandwich)
  se_beta   <- .fe_se(Xw, e_hat, N_obs, k_x, robust = robust)

  ## ── Stage 2: Decompose unit effects ─────────────────────────────────────
  ## Get one obs per panel: panel-level Z and u_bar
  tag_idx  <- !duplicated(uid)
  Z_cross  <- Z[tag_idx, , drop = FALSE]
  uid_cross <- uid[tag_idx]
  u_cross  <- u_bar_i[as.character(uid_cross)]

  ## OLS: u_bar_i ~ z_i + intercept  →  gamma_hat, chat_i
  stage2   <- .ols_fit(as.numeric(u_cross), cbind(1, Z_cross))
  gamma_hat <- stage2$coef[-1]
  alpha_hat <- stage2$coef[1]
  chat_i    <- as.numeric(u_cross) - cbind(1, Z_cross) %*% stage2$coef

  names(gamma_hat) <- zinvariants

  ## ── Pesaran-Zhou variance (Eq. 17) ──────────────────────────────────────
  ## Q_zz, V_zz, Q_zxbar
  Xbar_i   <- do.call(rbind, lapply(panels, function(p) {
    colMeans(X[uid == p, , drop = FALSE])
  }))

  Zc    <- scale(Z_cross, scale = FALSE)
  Xbc   <- scale(Xbar_i,  scale = FALSE)
  Qzz   <- crossprod(Zc) / N_g
  Vzz   <- crossprod(Zc * as.numeric(chat_i)^0.5, Zc * as.numeric(chat_i)^0.5) / N_g
  # Correct V_zz: sum_i chat_i^2 * zc_i zc_i'
  Vzz   <- matrix(0, k_z, k_z)
  for (i in seq_len(N_g)) {
    Vzz <- Vzz + chat_i[i]^2 * outer(Zc[i, ], Zc[i, ])
  }
  Vzz   <- Vzz / N_g
  Qzxbar <- crossprod(Zc, Xbc) / N_g

  ## FE variance (with robust option)
  if (robust) {
    V_beta <- .fe_var_robust(Xw, e_hat, N_obs, k_x)
  } else {
    V_beta <- sigma2_e * solve(crossprod(Xw))
  }

  ## V_gamma_pz = (1/N_g) * Qzz^{-1} * [Vzz + Qzxbar * N_g*V_beta * Qzxbar'] * Qzz^{-1}
  Qzz_inv    <- solve(Qzz)
  mid        <- Vzz + Qzxbar %*% (N_g * V_beta) %*% t(Qzxbar)
  V_gamma_pz <- (1 / N_g) * Qzz_inv %*% mid %*% Qzz_inv
  se_gamma   <- sqrt(pmax(diag(V_gamma_pz), 0))
  names(se_gamma) <- zinvariants

  ## Alpha SE from stage 2 OLS
  se_alpha <- sqrt(stage2$var_coef[1, 1])

  ## ── FEVD delta and raw variance ──────────────────────────────────────────
  delta        <- NULL
  V_gamma_fevd <- NULL
  if (method == "fevd") {
    h_i   <- as.numeric(chat_i)
    h_full <- h_i[match(uid, uid_cross)]
    stage3 <- .ols_fit(Y, cbind(X, Z, h_full, 1))
    delta  <- stage3$coef[k_x + k_z + 1]

    ## Raw FEVD SE for gamma (inconsistent, for diagnostic only)
    k3    <- k_x + k_z + 2  # X, Z, h, intercept
    V3    <- stage3$var_coef
    V_gamma_fevd <- V3[(k_x + 1):(k_x + k_z), (k_x + 1):(k_x + k_z), drop = FALSE]
  }

  ## ── FEF-IV: 2SLS stage 2 ─────────────────────────────────────────────────
  if (method == "fefiv") {
    R_cross  <- as.matrix(data_c[tag_idx, instruments, drop = FALSE])
    iv_fit   <- .tsls_fit(as.numeric(u_cross), Z_cross, R_cross)
    gamma_hat <- iv_fit$coef
    alpha_hat <- iv_fit$alpha
    names(gamma_hat) <- zinvariants

    ## IV residuals
    upsilon_i <- as.numeric(u_cross) - Z_cross %*% gamma_hat - alpha_hat

    ## PZ-IV variance (Eq. 51)
    V_gamma_pz <- .pz_iv_variance(Z_cross, R_cross, Xbar_i, upsilon_i,
                                   N_g, k_z, V_beta)
    se_gamma   <- sqrt(pmax(diag(V_gamma_pz), 0))
    names(se_gamma) <- zinvariants
    se_alpha   <- sqrt(iv_fit$var_alpha)
  }

  ## ── Assemble output ───────────────────────────────────────────────────────
  method_str <- switch(method,
    fevd  = "FEVD",
    fef   = "FEF",
    fefiv = "FEF-IV"
  )

  out <- list(
    beta_fe      = beta_fe,
    gamma        = gamma_hat,
    alpha        = alpha_hat,
    se_beta      = se_beta,
    se_gamma     = se_gamma,
    se_alpha     = se_alpha,
    V_gamma_pz   = V_gamma_pz,
    V_gamma_fevd = V_gamma_fevd,
    delta        = delta,
    sigma2_e     = sigma2_e,
    sigma2_u     = sigma2_u,
    N            = N_obs,
    N_g          = N_g,
    T_avg        = T_avg,
    method       = method_str,
    depvar       = depvar,
    xvars        = xvars,
    zinvariants  = zinvariants,
    instruments  = instruments,
    robust       = robust
  )
  class(out) <- "xtfifevd"
  out
}


## ── Internal helpers ─────────────────────────────────────────────────────

#' @keywords internal
.within_demean <- function(y, uid) {
  gm <- tapply(y, uid, mean)
  as.numeric(y) - as.numeric(gm[as.character(uid)])
}

#' @keywords internal
.within_demean_mat <- function(X, uid) {
  uid_levels <- sort(unique(uid))
  gm <- do.call(rbind, lapply(uid_levels, function(p) {
    matrix(colMeans(X[uid == p, , drop = FALSE]), nrow = 1L)
  }))
  rownames(gm) <- as.character(uid_levels)
  X - gm[as.character(uid), , drop = FALSE]
}

#' @keywords internal
.ols_fit <- function(y, X) {
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(X)
  qr_d    <- qr(X)
  coef    <- as.numeric(qr.coef(qr_d, y))
  resid   <- y - as.numeric(X %*% coef)
  sigma2  <- sum(resid^2) / max(n - p, 1)
  XtX_inv <- tryCatch(chol2inv(chol(crossprod(X))), error = function(e) solve(crossprod(X)))
  var_c   <- sigma2 * XtX_inv
  list(coef = as.numeric(coef), resid = as.numeric(resid),
       sigma2 = sigma2, var_coef = var_c)
}

#' @keywords internal
.fe_se <- function(Xw, e_hat, N_obs, k_x, robust = FALSE) {
  if (robust) {
    V <- .fe_var_robust(Xw, e_hat, N_obs, k_x)
  } else {
    sig2 <- sum(e_hat^2) / max(N_obs - k_x, 1)
    V    <- sig2 * tryCatch(solve(crossprod(Xw)), error = function(e) diag(k_x))
  }
  sqrt(pmax(diag(V), 0))
}

#' @keywords internal
.fe_var_robust <- function(Xw, e_hat, N_obs, k_x) {
  XtX_inv <- tryCatch(solve(crossprod(Xw)), error = function(e) diag(k_x))
  meat    <- crossprod(Xw * e_hat)
  XtX_inv %*% meat %*% XtX_inv
}

#' @keywords internal
.tsls_fit <- function(y, Z, R) {
  Rmat <- cbind(1, R)
  Zmat <- cbind(1, Z)
  # First stage: Z ~ R
  Zhat <- Rmat %*% stats::lm.fit(Rmat, Zmat)$coefficients
  # Second stage: y ~ Zhat
  fit2 <- stats::lm.fit(Zhat, y)
  coef_full <- as.numeric(fit2$coefficients)
  alpha <- coef_full[1]
  gamma <- coef_full[-1]
  # IV residuals
  resid <- y - Zmat %*% coef_full
  n  <- length(y)
  p  <- ncol(Zmat)
  sig2 <- sum(resid^2) / max(n - p, 1)
  XtX_inv <- tryCatch(solve(crossprod(Zhat)), error = function(e) diag(p))
  var_c   <- sig2 * XtX_inv
  list(coef = gamma, alpha = alpha, resid = resid,
       var_alpha = var_c[1, 1], var_gamma = var_c[-1, -1, drop = FALSE])
}

#' @keywords internal
.pz_iv_variance <- function(Z_cross, R_cross, Xbar_i, upsilon_i, N_g, k_z, V_beta) {
  k_iv  <- ncol(R_cross)
  Zc    <- scale(Z_cross, scale = FALSE)
  Rc    <- scale(R_cross, scale = FALSE)
  Xbc   <- scale(Xbar_i,  scale = FALSE)

  Qzr   <- crossprod(Zc, Rc) / N_g
  Qrr   <- crossprod(Rc) / N_g
  Qrxbar <- crossprod(Rc, Xbc) / N_g

  Vrr   <- matrix(0, k_iv, k_iv)
  for (i in seq_len(N_g)) {
    Vrr <- Vrr + upsilon_i[i]^2 * outer(Rc[i, ], Rc[i, ])
  }
  Vrr <- Vrr / N_g

  Qrr_inv <- solve(Qrr)
  Hzr     <- solve(Qzr %*% Qrr_inv %*% t(Qzr)) %*% Qzr %*% Qrr_inv
  mid_iv  <- Vrr + Qrxbar %*% (N_g * V_beta) %*% t(Qrxbar)
  (1 / N_g) * Hzr %*% mid_iv %*% t(Hzr)
}


## ── S3 methods ───────────────────────────────────────────────────────────

#' Print method for xtfifevd objects
#'
#' @param x An object of class \code{"xtfifevd"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtfifevd <- function(x, ...) {
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("  %s Estimation Results (Pesaran-Zhou corrected SEs)\n", x$method))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  Dep. variable : %s\n", x$depvar))
  cat(sprintf("  Observations  : %d\n", x$N))
  cat(sprintf("  Groups        : %d\n", x$N_g))
  cat(sprintf("  T (avg)       : %.1f\n", x$T_avg))
  cat(sprintf("  Method        : %s\n", x$method))
  cat(strrep("-", 70), "\n\n")

  ## Time-varying part (FE)
  cat("  Time-varying regressors (FE estimates):\n")
  cat(sprintf("  %-18s  %12s  %10s  %8s  %s\n",
              "Variable", "Estimate", "Std. Error", "z value", "p-value"))
  cat(strrep("-", 70), "\n")
  for (j in seq_along(x$xvars)) {
    z_val <- x$beta_fe[j] / max(x$se_beta[j], 1e-14)
    p_val <- 2 * (1 - stats::pnorm(abs(z_val)))
    cat(sprintf("  %-18s  %12.6f  %10.6f  %8.3f  %.4f  %s\n",
                x$xvars[j], x$beta_fe[j], x$se_beta[j], z_val, p_val,
                .fifevd_stars(p_val)))
  }
  cat("\n")

  ## Time-invariant part (PZ corrected)
  cat(sprintf("  Time-invariant regressors (%s estimates, PZ-corrected SEs):\n",
              x$method))
  cat(sprintf("  %-18s  %12s  %10s  %8s  %s\n",
              "Variable", "Estimate", "Std. Error", "z value", "p-value"))
  cat(strrep("-", 70), "\n")
  for (j in seq_along(x$zinvariants)) {
    z_val <- x$gamma[j] / max(x$se_gamma[j], 1e-14)
    p_val <- 2 * (1 - stats::pnorm(abs(z_val)))
    cat(sprintf("  %-18s  %12.6f  %10.6f  %8.3f  %.4f  %s\n",
                x$zinvariants[j], x$gamma[j], x$se_gamma[j], z_val, p_val,
                .fifevd_stars(p_val)))
  }
  # Intercept
  z_a <- x$alpha / max(x$se_alpha, 1e-14)
  p_a <- 2 * (1 - stats::pnorm(abs(z_a)))
  cat(sprintf("  %-18s  %12.6f  %10.6f  %8.3f  %.4f  %s\n",
              "_cons", x$alpha, x$se_alpha, z_a, p_a,
              .fifevd_stars(p_a)))
  cat(strrep("-", 70), "\n")
  cat("  *** p<0.01, ** p<0.05, * p<0.10\n")
  if (!is.null(x$delta)) {
    cat(sprintf("  Stage-3 delta: %.6f (expected: 1.000000)\n", x$delta))
  }
  cat("\n")
  invisible(x)
}

#' Summary method for xtfifevd objects
#'
#' @param object An object of class \code{"xtfifevd"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtfifevd <- function(object, ...) {
  print(object, ...)
  if (!is.null(object$V_gamma_fevd)) {
    cat("  Comparison: Raw FEVD SEs (inconsistent) vs PZ-corrected SEs:\n")
    cat(sprintf("  %-18s  %12s  %12s  %8s\n",
                "Variable", "PZ SE", "Raw SE", "Ratio"))
    cat(strrep("-", 60), "\n")
    for (j in seq_along(object$zinvariants)) {
      se_pz  <- object$se_gamma[j]
      se_raw <- sqrt(pmax(object$V_gamma_fevd[j, j], 0))
      ratio  <- if (se_raw > 1e-14) se_pz / se_raw else NA_real_
      cat(sprintf("  %-18s  %12.6f  %12.6f  %8.2f\n",
                  object$zinvariants[j], se_pz, se_raw, ratio))
    }
    cat(strrep("-", 60), "\n")
    cat("  Note: Raw FEVD SEs ignore generated-regressor uncertainty\n")
    cat("        (Pesaran & Zhou 2016, Remark 4).\n\n")
  }
  cat(sprintf("  sigma^2_e (idiosyncratic): %.6f\n", object$sigma2_e))
  cat(sprintf("  sigma^2_u (unit effects) : %.6f\n", object$sigma2_u))
  invisible(object)
}

#' @keywords internal
.fifevd_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.10) return("*")
  return("")
}

