#' Tests of No Cross-Sectional Dependence in Panel Quantile Regressions
#'
#' Tests the null hypothesis of no cross-sectional error dependence (CSD) in
#' panel quantile regressions. Implements the T_tau and T-tilde_tau statistics
#' of Demetrescu, Hosseinkouchack and Rodrigues (2023).
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...}. Required for
#'   \code{mode = "pooled"} (default) and \code{mode = "individual"}.
#'   Not used when \code{residuals} is provided.
#' @param data A data frame containing the panel data in long format. Required
#'   unless \code{residuals} is provided.
#' @param index A character vector of length 2: \code{c("id_var", "time_var")}.
#'   Required unless \code{residuals} is provided.
#' @param quantiles A numeric vector of quantile levels, each strictly between
#'   0 and 1.
#' @param mode Estimation mode: \code{"pooled"} (default, pooled FE-QR),
#'   \code{"individual"} (per-unit QR), or \code{"residuals"} (provide
#'   pre-computed residuals via the \code{residuals} argument).
#' @param residuals A list (or named list) of numeric vectors or a matrix with
#'   one column per quantile, containing pre-computed QR residuals. Only used
#'   when \code{mode = "residuals"}.
#' @param bandwidth Numeric. KDE bandwidth for sparsity estimation. If
#'   \code{NULL} (default), uses \eqn{0.35 (NT)^{-0.2}}.
#' @param correction Logical. If \code{TRUE} (default), reports the
#'   bias-corrected T-tilde statistic in addition to T_tau.
#'
#' @return An object of class \code{"xtcsdq"} with components:
#'   \describe{
#'     \item{T_tau}{Numeric vector of T_tau statistics (one per quantile).}
#'     \item{Ttilde_tau}{Numeric vector of bias-corrected T-tilde_tau statistics.}
#'     \item{pval_T}{p-values for T_tau.}
#'     \item{pval_Ttilde}{p-values for T-tilde_tau.}
#'     \item{fhat}{KDE density estimates at zero (one per quantile).}
#'     \item{M_K}{Portmanteau statistic (average of T_tau over quantiles).}
#'     \item{Mtilde_K}{Bias-corrected portmanteau statistic.}
#'     \item{pval_M}{p-value for M_K.}
#'     \item{pval_Mc}{p-value for Mtilde_K.}
#'     \item{quantiles}{Quantile levels used.}
#'     \item{N}{Number of cross-sectional units.}
#'     \item{TT}{Number of time periods.}
#'     \item{bandwidth}{KDE bandwidth used.}
#'   }
#'
#' @details
#' The T_tau statistic (Equation 3 in Demetrescu et al., 2023) tests for CSD
#' by examining pairwise correlations of demeaned QR residuals across units.
#' Under the null of no CSD, T_tau is asymptotically standard normal.
#'
#' The bias-corrected version T-tilde_tau (Equation 5) subtracts two correction
#' terms that account for the estimation uncertainty in the QR slope and the
#' sparsity at the quantile. Reject H0 for large positive values.
#'
#' The portmanteau statistic \eqn{M_K = K^{-1} \sum_{q=1}^K T_\tau^{(q)}}
#' aggregates across K quantile levels.
#'
#' The KDE bandwidth defaults to \eqn{h = 0.35 (NT)^{-0.2}} as in the
#' original paper.
#'
#' @references
#' Demetrescu, M., Hosseinkouchack, M. and Rodrigues, P.M.M. (2023).
#' Testing for No Cross-Sectional Error Dependence in Panel Quantile
#' Regressions.
#' \emph{Ruhr Economic Papers}, No. 1041.
#' \doi{10.4419/96973002}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 8; tt <- 20
#' dat <- data.frame(
#'   id   = rep(1:n, each = tt),
#'   time = rep(1:tt, times = n),
#'   y    = rnorm(n * tt),
#'   x1   = rnorm(n * tt)
#' )
#' res <- xtcsdq(y ~ x1, data = dat, index = c("id", "time"),
#'               quantiles = c(0.25, 0.5, 0.75))
#' print(res)
#' summary(res)
#' }
#'
#' @export
xtcsdq <- function(formula = NULL, data = NULL, index = NULL,
                   quantiles,
                   mode = c("pooled", "individual", "residuals"),
                   residuals = NULL,
                   bandwidth = NULL,
                   correction = TRUE) {

  mode <- match.arg(mode)

  ## ── Validate inputs ──────────────────────────────────────────────────────
  if (!is.numeric(quantiles) || any(quantiles <= 0) || any(quantiles >= 1)) {
    stop("'quantiles' must be numeric with all values strictly between 0 and 1.",
         call. = FALSE)
  }
  quantiles <- sort(unique(quantiles))
  K         <- length(quantiles)

  if (mode == "residuals") {
    if (is.null(residuals)) {
      stop("'residuals' must be provided when mode = 'residuals'.", call. = FALSE)
    }
    # Normalize residuals to a matrix: T x K
    if (is.data.frame(residuals) || is.matrix(residuals)) {
      resid_mat <- as.matrix(residuals)
    } else if (is.list(residuals)) {
      resid_mat <- do.call(cbind, residuals)
    } else if (is.numeric(residuals) && K == 1L) {
      resid_mat <- matrix(residuals, ncol = 1L)
    } else {
      stop("'residuals' must be a matrix, data frame, or list of vectors.",
           call. = FALSE)
    }
    if (ncol(resid_mat) != K) {
      stop(sprintf("Number of columns in 'residuals' (%d) must equal length of 'quantiles' (%d).",
                   ncol(resid_mat), K), call. = FALSE)
    }
    # For residuals mode we need N and T; require index + data OR infer from dim
    if (is.null(data) || is.null(index)) {
      stop("'data' and 'index' are required even in 'residuals' mode to determine N and T.",
           call. = FALSE)
    }
  } else {
    if (is.null(formula) || is.null(data) || is.null(index)) {
      stop("'formula', 'data', and 'index' are required for mode = '", mode, "'.",
           call. = FALSE)
    }
    if (!inherits(formula, "formula")) {
      stop("'formula' must be a formula object.", call. = FALSE)
    }
    if (!is.character(index) || length(index) != 2) {
      stop("'index' must be a character vector of length 2.", call. = FALSE)
    }
    if (!all(index %in% names(data))) {
      stop("Variables in 'index' not found in 'data'.", call. = FALSE)
    }
  }

  ## ── Panel structure ───────────────────────────────────────────────────────
  ivar <- index[1]
  tvar <- index[2]

  panels <- sort(unique(data[[ivar]]))
  N      <- length(panels)
  times  <- sort(unique(data[[tvar]]))
  TT     <- length(times)

  if (N < 2) stop("At least 2 cross-sectional units are required.", call. = FALSE)
  if (TT < 3) stop("At least 3 time periods are required.", call. = FALSE)

  ## Check balance
  obs_count <- as.integer(table(data[[ivar]]))
  if (length(unique(obs_count)) > 1) {
    stop("Panel must be balanced (all units must have the same number of observations).",
         call. = FALSE)
  }

  ## ── KDE bandwidth ─────────────────────────────────────────────────────────
  if (is.null(bandwidth) || bandwidth <= 0) {
    bw_used <- 0.35 * (N * TT)^(-0.2)
  } else {
    bw_used <- bandwidth
  }

  ## ── Get/arrange residuals ─────────────────────────────────────────────────
  # resid_list[[q]] is an N x TT matrix of demeaned residuals for quantile q
  resid_list <- vector("list", K)

  if (mode == "residuals") {
    ## Residuals provided: reshape to N x TT per quantile
    for (qi in seq_len(K)) {
      r_vec <- resid_mat[, qi]
      Rmat  <- .reshape_resid(r_vec, data, ivar, tvar, panels, N, TT)
      resid_list[[qi]] <- Rmat
    }
  } else if (mode == "individual") {
    mf_vars <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
    depv    <- names(mf_vars)[1]
    indv    <- names(mf_vars)[-1]
    for (qi in seq_len(K)) {
      raw_resid <- rep(NA_real_, nrow(data))
      for (pi in seq_along(panels)) {
        idx      <- which(data[[ivar]] == panels[pi])
        sub_data <- data[idx, , drop = FALSE]
        sub_mf   <- tryCatch(
          stats::model.frame(formula, data = sub_data, na.action = stats::na.omit),
          error = function(e) NULL
        )
        if (is.null(sub_mf) || nrow(sub_mf) < length(indv) + 2) next
        Y_i <- sub_mf[[1]]
        X_i <- as.matrix(sub_mf[, -1, drop = FALSE])
        fit <- tryCatch(.qr_irls(Y_i, cbind(1, X_i), quantiles[qi]),
                        error = function(e) NULL)
        if (is.null(fit)) next
        raw_resid[idx] <- Y_i - cbind(1, X_i) %*% fit
      }
      Rmat <- .reshape_resid(raw_resid, data, ivar, tvar, panels, N, TT)
      resid_list[[qi]] <- Rmat
    }
  } else {
    ## pooled FE-QR
    mf_vars <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
    depv    <- names(mf_vars)[1]
    indv    <- names(mf_vars)[-1]
    data_s  <- data[order(data[[ivar]], data[[tvar]]), , drop = FALSE]
    # Create dummy FE dummies (within-demeaning approach: subtract unit means)
    for (qi in seq_len(K)) {
      Y_vec <- mf_vars[[1]]
      X_mat <- as.matrix(mf_vars[, -1, drop = FALSE])
      # Within transformation
      uid  <- data[[ivar]]
      Yw   <- Y_vec - tapply(Y_vec, uid, mean)[as.character(uid)]
      group_means <- do.call(rbind, lapply(split(seq_len(nrow(X_mat)), uid), function(i) { matrix(colMeans(X_mat[i, , drop = FALSE]), 1) }))
  rownames(group_means) <- names(split(seq_len(nrow(X_mat)), uid))
  Xw   <- X_mat - group_means[as.character(uid), , drop = FALSE]

      fit  <- tryCatch(.qr_irls(as.numeric(Yw), cbind(1, Xw), quantiles[qi]),
                       error = function(e) NULL)
      if (is.null(fit)) {
        raw_resid <- rep(NA_real_, nrow(data))
      } else {
        raw_resid <- as.numeric(Yw) - cbind(1, Xw) %*% fit
      }
      Rmat <- .reshape_resid(raw_resid, data, ivar, tvar, panels, N, TT)
      resid_list[[qi]] <- Rmat
    }
  }

  ## ── Compute test statistics per quantile ──────────────────────────────────
  T_tau_vec    <- numeric(K)
  Ttilde_vec   <- numeric(K)
  pval_T_vec   <- numeric(K)
  pval_Tt_vec  <- numeric(K)
  fhat_vec     <- numeric(K)

  for (qi in seq_len(K)) {
    Rmat <- resid_list[[qi]]
    if (is.null(Rmat) || any(is.na(Rmat))) {
      T_tau_vec[qi]   <- NA_real_
      Ttilde_vec[qi]  <- NA_real_
      pval_T_vec[qi]  <- NA_real_
      pval_Tt_vec[qi] <- NA_real_
      fhat_vec[qi]    <- NA_real_
      next
    }

    ## Demean: subtract unit means
    unit_means <- rowMeans(Rmat)
    Rmd        <- Rmat - unit_means  # N x TT demeaned

    ## Standardize by unit SD for KDE
    unit_sd <- apply(Rmd, 1, stats::sd)
    unit_sd[unit_sd < 1e-14] <- 1e-14
    Rstd    <- Rmat / unit_sd   # N x TT (raw / sigma_i)

    ## KDE estimate of density at zero
    all_std <- as.numeric(Rstd)
    h       <- bw_used
    fhat    <- mean(stats::dnorm(all_std / h) / h)
    fhat_vec[qi] <- fhat

    ## Compute T_tau via pairwise correlations on demeaned residuals
    T_tau_val <- .compute_T_tau(Rmd, N, TT)

    ## Bias corrections
    corr1 <- sqrt(N * (N - 1) / (2 * TT))
    tau_q  <- quantiles[qi]
    corr2  <- (tau_q * (1 - tau_q)) / max(fhat^2, 1e-14) *
              sqrt(N * (N - 1) / TT)
    Ttilde_val <- T_tau_val - corr1 - corr2

    T_tau_vec[qi]   <- T_tau_val
    Ttilde_vec[qi]  <- Ttilde_val
    pval_T_vec[qi]  <- 1 - stats::pnorm(T_tau_val)
    pval_Tt_vec[qi] <- 1 - stats::pnorm(Ttilde_val)
  }

  ## ── Portmanteau M_K ───────────────────────────────────────────────────────
  M_K      <- mean(T_tau_vec, na.rm = TRUE)
  Mtilde_K <- mean(Ttilde_vec, na.rm = TRUE)
  pval_M   <- 1 - stats::pnorm(M_K)
  pval_Mc  <- 1 - stats::pnorm(Mtilde_K)

  ## ── Output ────────────────────────────────────────────────────────────────
  out <- list(
    T_tau      = T_tau_vec,
    Ttilde_tau = Ttilde_vec,
    pval_T     = pval_T_vec,
    pval_Ttilde = pval_Tt_vec,
    fhat       = fhat_vec,
    M_K        = M_K,
    Mtilde_K   = Mtilde_K,
    pval_M     = pval_M,
    pval_Mc    = pval_Mc,
    quantiles  = quantiles,
    N          = N,
    TT         = TT,
    bandwidth  = bw_used,
    correction = correction,
    mode       = mode
  )
  class(out) <- "xtcsdq"
  out
}


## ── Internal helpers ─────────────────────────────────────────────────────

#' @keywords internal
.reshape_resid <- function(r_vec, data, ivar, tvar, panels, N, TT) {
  Rmat <- matrix(NA_real_, nrow = N, ncol = TT)
  times <- sort(unique(data[[tvar]]))
  for (pi in seq_along(panels)) {
    idx <- which(data[[ivar]] == panels[pi])
    sub <- data[idx, c(ivar, tvar), drop = FALSE]
    sub <- sub[order(sub[[tvar]]), , drop = FALSE]
    ti  <- match(sub[[tvar]], times)
    Rmat[pi, ti] <- r_vec[idx[order(data[idx, tvar])]]
  }
  Rmat
}

#' @keywords internal
.compute_T_tau <- function(Rmd, N, TT) {
  ## T_tau = sum_{i<j} (T * rho_ij^2 - 1) / sqrt(N(N-1))
  ## rho_ij = cor of demeaned resids i and j
  rho2_sum <- 0
  for (i in seq_len(N - 1)) {
    u_i  <- Rmd[i, ]
    ss_i <- sum(u_i^2)
    for (j in (i + 1):N) {
      u_j    <- Rmd[j, ]
      ss_j   <- sum(u_j^2)
      cov_ij <- sum(u_i * u_j)
      denom  <- ss_i * ss_j
      if (denom < 1e-30) next
      rho_sq   <- cov_ij^2 / denom
      rho2_sum <- rho2_sum + (TT * rho_sq - 1)
    }
  }
  rho2_sum / sqrt(N * (N - 1))
}

#' @keywords internal
.qr_irls <- function(y, X, tau, maxit = 200L, tol = 1e-8) {
  n    <- length(y)
  beta <- stats::lm.fit(X, y)$coefficients
  for (iter in seq_len(maxit)) {
    resid    <- y - X %*% beta
    w        <- pmax(abs(resid), tol)
    w        <- ifelse(resid >= 0, tau / w, (1 - tau) / w)
    W        <- diag(as.numeric(w), n)
    XtWX     <- crossprod(X, W %*% X)
    XtWy     <- crossprod(X, W %*% y)
    beta_new <- tryCatch(solve(XtWX, XtWy), error = function(e) beta)
    if (max(abs(beta_new - beta)) < tol) { beta <- beta_new; break }
    beta <- beta_new
  }
  as.numeric(beta)
}


## ── S3 methods ───────────────────────────────────────────────────────────

#' Print method for xtcsdq objects
#'
#' @param x An object of class \code{"xtcsdq"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtcsdq <- function(x, ...) {
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("  XTCSDQ: CSD Test in Panel Quantile Regressions\n")
  cat("  Demetrescu, Hosseinkouchack & Rodrigues (2023)\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("  N (units)  : %d\n", x$N))
  cat(sprintf("  T (periods): %d\n", x$TT))
  cat(sprintf("  Quantiles  : %s\n", paste(x$quantiles, collapse = ", ")))
  cat(sprintf("  Bandwidth  : %.6f\n", x$bandwidth))
  cat(sprintf("  Mode       : %s\n", x$mode))
  cat(strrep("-", 70), "\n")
  cat("  H0: No cross-sectional error dependence in panel QR\n")
  cat("  (Reject for large positive values)\n")
  cat(strrep("-", 70), "\n")

  if (x$correction) {
    cat(sprintf("  %-8s  %10s  %10s  %10s  %10s\n",
                "tau", "T_tau", "T~_tau", "p(T)", "p(T~)"))
  } else {
    cat(sprintf("  %-8s  %10s  %10s\n", "tau", "T_tau", "p(T)"))
  }
  cat(strrep("-", 70), "\n")

  for (i in seq_along(x$quantiles)) {
    if (x$correction) {
      cat(sprintf("  %-8.3f  %10.4f  %10.4f  %10.4f  %10.4f  %s\n",
                  x$quantiles[i], x$T_tau[i], x$Ttilde_tau[i],
                  x$pval_T[i], x$pval_Ttilde[i],
                  .csd_stars(x$pval_Ttilde[i])))
    } else {
      cat(sprintf("  %-8.3f  %10.4f  %10.4f  %s\n",
                  x$quantiles[i], x$T_tau[i], x$pval_T[i],
                  .csd_stars(x$pval_T[i])))
    }
  }

  if (length(x$quantiles) > 1) {
    cat(strrep("-", 70), "\n")
    if (x$correction) {
      cat(sprintf("  %-8s  %10.4f  %10.4f  %10.4f  %10.4f  %s\n",
                  "M_K", x$M_K, x$Mtilde_K, x$pval_M, x$pval_Mc,
                  .csd_stars(x$pval_Mc)))
    } else {
      cat(sprintf("  %-8s  %10.4f  %10.4f  %s\n",
                  "M_K", x$M_K, x$pval_M, .csd_stars(x$pval_M)))
    }
  }
  cat(strrep("-", 70), "\n")
  cat("  *** p<0.01, ** p<0.05, * p<0.10\n\n")

  invisible(x)
}

#' Summary method for xtcsdq objects
#'
#' @param object An object of class \code{"xtcsdq"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtcsdq <- function(object, ...) {
  print(object, ...)
  cat("  KDE density estimates at zero (fhat):\n")
  for (i in seq_along(object$quantiles)) {
    cat(sprintf("    tau=%.2f: fhat=%.6f\n", object$quantiles[i], object$fhat[i]))
  }
  cat("\n")
  invisible(object)
}

#' @keywords internal
.csd_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.10) return("*")
  return("")
}

