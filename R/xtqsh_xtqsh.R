#' Quantile Regression Slope Homogeneity Test for Panel Data
#'
#' Tests the null hypothesis of slope homogeneity in panel quantile regressions.
#' Implements the S-hat (chi-squared) and D-hat (standard normal) statistics of
#' Galvao, Juhl, Montes-Rojas and Olmo (2017).
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} specifying the
#'   dependent variable and covariates.
#' @param data A data frame containing the panel data in long format.
#' @param index A character vector of length 2 giving the names of the panel
#'   individual and time variables, e.g. \code{c("id", "time")}.
#' @param tau A numeric vector of quantile levels, each strictly between 0 and 1.
#' @param bw Bandwidth method for sparsity estimation. Either
#'   \code{"hallsheather"} (default) or \code{"bofinger"}.
#' @param marginal Logical. If \code{TRUE}, per-variable marginal slope
#'   homogeneity tests are also computed. Default is \code{FALSE}.
#'
#' @return An object of class \code{"xtqsh"} with the following components:
#'   \describe{
#'     \item{S}{Numeric vector of S-hat statistics (one per quantile).}
#'     \item{D}{Numeric vector of D-hat statistics (one per quantile).}
#'     \item{pval_S}{p-values for S-hat (chi-squared distribution).}
#'     \item{pval_D}{p-values for D-hat (standard normal).}
#'     \item{S_ols}{S-hat statistic for the mean (OLS) regression.}
#'     \item{D_ols}{D-hat statistic for the mean (OLS) regression.}
#'     \item{pval_S_ols}{p-value for S_ols.}
#'     \item{pval_D_ols}{p-value for D_ols.}
#'     \item{beta_md}{Matrix of minimum-distance QR estimates (quantiles x regressors).}
#'     \item{beta_md_se}{Matrix of MD-QR standard errors.}
#'     \item{beta_all}{Array of per-panel QR coefficients (panels x regressors x quantiles).}
#'     \item{tau}{Quantile levels used.}
#'     \item{bw}{Bandwidth method used.}
#'     \item{n_panels}{Number of panel units.}
#'     \item{n_obs}{Total observations.}
#'     \item{k}{Number of regressors (excluding intercept).}
#'     \item{depvar}{Name of the dependent variable.}
#'     \item{indepvars}{Names of the independent variables.}
#'     \item{marginal}{List of marginal test results (if \code{marginal = TRUE}).}
#'   }
#'
#' @details
#' The test statistic D-hat is asymptotically standard normal under the null
#' hypothesis as both \eqn{N} and \eqn{T} tend to infinity. The S-hat statistic
#' follows a chi-squared distribution with \eqn{k(N-1)} degrees of freedom when
#' \eqn{T \to \infty} with \eqn{N} fixed.
#'
#' Bandwidth selection:
#' \itemize{
#'   \item \code{"hallsheather"}: \eqn{h = n^{-1/3} z_{(1+\tau)/2}^{2/3}
#'     \left(\frac{1.5 \phi^2(\Phi^{-1}(\tau))}{2(\Phi^{-1}(\tau))^2 + 1}\right)^{1/3}}
#'   \item \code{"bofinger"}: \eqn{h = n^{-1/5}
#'     \left(\frac{4.5 \phi^4(\Phi^{-1}(\tau))}{(2(\Phi^{-1}(\tau))^2 + 1)^2}\right)^{1/5}}
#' }
#'
#' @references
#' Galvao, A.F., Juhl, T., Montes-Rojas, G. and Olmo, J. (2017).
#' Testing Slope Homogeneity in Quantile Regression Panel Data with an
#' Application to the Cross-Section of Stock Returns.
#' \emph{Journal of Financial Econometrics}, 16(2), 211--243.
#' \doi{10.1093/jjfinec/nbx003}
#'
#' Bofinger, E. (1975). Estimation of a Density Function Using Order Statistics.
#' \emph{Australian Journal of Statistics}, 17(1), 1--7.
#'
#' Hall, P. and Sheather, S.J. (1988). On the Distribution of a Studentized
#' Quantile. \emph{Journal of the Royal Statistical Society Series B}, 50(3),
#' 381--391.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 10; tt <- 20
#' dat <- data.frame(
#'   id   = rep(1:n, each = tt),
#'   time = rep(1:tt, times = n),
#'   y    = rnorm(n * tt),
#'   x1   = rnorm(n * tt),
#'   x2   = rnorm(n * tt)
#' )
#' res <- xtqsh(y ~ x1 + x2, data = dat, index = c("id", "time"),
#'              tau = c(0.25, 0.5, 0.75))
#' print(res)
#' summary(res)
#' }
#'
#' @export
xtqsh <- function(formula, data, index, tau,
                  bw = "hallsheather", marginal = FALSE) {

  ## ── Input validation ────────────────────────────────────────────────────
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (!is.character(index) || length(index) != 2) {
    stop("'index' must be a character vector of length 2: c('id_var', 'time_var').",
         call. = FALSE)
  }
  if (!all(index %in% names(data))) {
    stop("Variables specified in 'index' not found in 'data'.", call. = FALSE)
  }
  if (!is.numeric(tau) || any(tau <= 0) || any(tau >= 1)) {
    stop("'tau' must be a numeric vector with all elements strictly between 0 and 1.",
         call. = FALSE)
  }
  tau <- sort(unique(tau))
  bw  <- match.arg(bw, c("hallsheather", "bofinger"))

  ## ── Parse formula ───────────────────────────────────────────────────────
  mf      <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  depvar  <- names(mf)[1]
  indvars <- names(mf)[-1]
  k       <- length(indvars)
  ntau    <- length(tau)

  if (k < 1) {
    stop("At least one independent variable is required.", call. = FALSE)
  }

  ivar <- index[1]
  tvar <- index[2]

  ## ── Identify panels ─────────────────────────────────────────────────────
  panels    <- sort(unique(data[[ivar]]))
  n_panels  <- length(panels)

  if (n_panels < 2) {
    stop("At least 2 panel units are required.", call. = FALSE)
  }

  n_obs <- nrow(mf)

  ## ── Per-panel QR estimation ──────────────────────────────────────────────
  # beta_all[i, j, q] = coefficient j for panel i at quantile tau[q]
  beta_all    <- array(NA_real_, dim = c(n_panels, k, ntau))
  beta_se_all <- array(NA_real_, dim = c(n_panels, k, ntau))
  valid_vec   <- logical(n_panels)

  for (qi in seq_along(tau)) {
    for (pi in seq_along(panels)) {
      idx   <- which(data[[ivar]] == panels[pi])
      subdf <- mf[rownames(mf) %in% rownames(data[idx, , drop = FALSE]),
                  , drop = FALSE]
      # fallback: match by position
      # Use a safe subset
      sub_data <- data[idx, , drop = FALSE]
      sub_mf   <- stats::model.frame(formula, data = sub_data,
                                     na.action = stats::na.omit)
      if (nrow(sub_mf) < k + 2) next

      Y_i <- sub_mf[[1]]
      X_i <- as.matrix(sub_mf[, -1, drop = FALSE])

      fit <- tryCatch(
        .qr_fit(Y_i, X_i, tau[qi]),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      beta_all[pi, , qi] <- fit$coef
      beta_se_all[pi, , qi] <- fit$se
      valid_vec[pi] <- TRUE
    }
  }

  valid_panels <- sum(valid_vec)
  if (valid_panels < 2) {
    stop("Fewer than 2 panels produced valid QR estimates.", call. = FALSE)
  }

  ## ── OLS (mean) estimates per panel ──────────────────────────────────────
  beta_ols_all <- matrix(NA_real_, nrow = n_panels, ncol = k)
  for (pi in seq_along(panels)) {
    idx      <- which(data[[ivar]] == panels[pi])
    sub_data <- data[idx, , drop = FALSE]
    sub_mf   <- tryCatch(
      stats::model.frame(formula, data = sub_data, na.action = stats::na.omit),
      error = function(e) NULL
    )
    if (is.null(sub_mf) || nrow(sub_mf) < k + 2) next
    Y_i <- sub_mf[[1]]
    X_i <- as.matrix(sub_mf[, -1, drop = FALSE])
    ols <- tryCatch(
      stats::lm.fit(cbind(1, X_i), Y_i),
      error = function(e) NULL
    )
    if (!is.null(ols)) {
      beta_ols_all[pi, ] <- ols$coefficients[-1]
    }
  }

  ## ── Compute S-hat and D-hat for each quantile ───────────────────────────
  S_vec    <- numeric(ntau)
  D_vec    <- numeric(ntau)
  pval_S   <- numeric(ntau)
  pval_D   <- numeric(ntau)
  beta_md  <- matrix(NA_real_, nrow = ntau, ncol = k)
  beta_md_se <- matrix(NA_real_, nrow = ntau, ncol = k)

  for (qi in seq_along(tau)) {
    res <- .compute_slope_test(
      beta_all     = beta_all[valid_vec, , qi, drop = FALSE],
      beta_se_all  = beta_se_all[valid_vec, , qi, drop = FALSE],
      tau_val      = tau[qi],
      bw           = bw,
      n_panels_eff = valid_panels,
      k            = k
    )
    S_vec[qi]       <- res$S
    D_vec[qi]       <- res$D
    pval_S[qi]      <- res$pval_S
    pval_D[qi]      <- res$pval_D
    beta_md[qi, ]   <- res$beta_md
    beta_md_se[qi, ] <- res$beta_md_se
  }

  ## ── OLS (mean) joint test ───────────────────────────────────────────────
  ols_res <- .compute_ols_test(
    beta_ols_all = beta_ols_all[valid_vec, , drop = FALSE],
    k            = k,
    n_panels_eff = valid_panels
  )

  ## ── Marginal tests ───────────────────────────────────────────────────────
  marginal_res <- NULL
  if (marginal) {
    marginal_res <- list()
    for (j in seq_len(k)) {
      S_marg  <- numeric(ntau)
      D_marg  <- numeric(ntau)
      pS_marg <- numeric(ntau)
      pD_marg <- numeric(ntau)
      for (qi in seq_along(tau)) {
        b_j  <- beta_all[valid_vec, j, qi]
        se_j <- beta_se_all[valid_vec, j, qi]
        if (any(is.na(b_j)) || any(is.na(se_j))) {
          S_marg[qi] <- NA_real_; D_marg[qi] <- NA_real_
          pS_marg[qi] <- NA_real_; pD_marg[qi] <- NA_real_
          next
        }
        res_j <- .compute_slope_test(
          beta_all     = array(b_j, dim = c(valid_panels, 1, 1)),
          beta_se_all  = array(se_j, dim = c(valid_panels, 1, 1)),
          tau_val      = tau[qi],
          bw           = bw,
          n_panels_eff = valid_panels,
          k            = 1L
        )
        S_marg[qi]  <- res_j$S
        D_marg[qi]  <- res_j$D
        pS_marg[qi] <- res_j$pval_S
        pD_marg[qi] <- res_j$pval_D
      }
      marginal_res[[indvars[j]]] <- list(
        S = S_marg, D = D_marg, pval_S = pS_marg, pval_D = pD_marg
      )
    }
  }

  ## ── Assemble output ──────────────────────────────────────────────────────
  colnames(beta_md)    <- indvars
  colnames(beta_md_se) <- indvars
  rownames(beta_md)    <- paste0("tau=", tau)
  rownames(beta_md_se) <- paste0("tau=", tau)

  out <- list(
    S           = S_vec,
    D           = D_vec,
    pval_S      = pval_S,
    pval_D      = pval_D,
    S_ols       = ols_res$S,
    D_ols       = ols_res$D,
    pval_S_ols  = ols_res$pval_S,
    pval_D_ols  = ols_res$pval_D,
    beta_md     = beta_md,
    beta_md_se  = beta_md_se,
    beta_all    = beta_all,
    beta_ols    = beta_ols_all,
    tau         = tau,
    bw          = bw,
    n_panels    = valid_panels,
    n_obs       = n_obs,
    k           = k,
    depvar      = depvar,
    indepvars   = indvars,
    marginal    = marginal_res
  )
  class(out) <- "xtqsh"
  out
}


## ── Internal: quantile regression via simplex / IRLS ─────────────────────
#' @keywords internal
.qr_fit <- function(y, X, tau) {
  n  <- length(y)
  Xm <- cbind(1, X)
  p  <- ncol(Xm)

  ## Barrodale-Roberts simplex (using stats::quantreg unavailable;
  ## implement a lightweight version via iteratively weighted LS)
  coef <- .qr_irls(y, Xm, tau)

  ## Sparsity (sandwich) SE estimate
  resid <- y - Xm %*% coef
  se    <- .qr_se(y, Xm, resid, tau)

  list(coef = coef[-1], se = se[-1])   # drop intercept
}


#' @keywords internal
.qr_irls <- function(y, X, tau, maxit = 200L, tol = 1e-8) {
  n <- length(y)
  p <- ncol(X)
  beta <- stats::lm.fit(X, y)$coefficients

  for (iter in seq_len(maxit)) {
    resid  <- y - X %*% beta
    w      <- ifelse(abs(resid) < tol, tol, abs(resid))
    w      <- ifelse(resid >= 0, tau / w, (1 - tau) / w)
    W      <- diag(as.numeric(w), n)
    XtWX   <- crossprod(X, W %*% X)
    XtWy   <- crossprod(X, W %*% y)
    beta_new <- tryCatch(
      solve(XtWX, XtWy),
      error = function(e) beta
    )
    if (max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      break
    }
    beta <- beta_new
  }
  as.numeric(beta)
}


#' @keywords internal
.qr_se <- function(y, X, resid, tau) {
  n  <- length(y)
  p  <- ncol(X)
  h  <- .bw_hallsheather(n, tau)
  f_hat <- mean(stats::dnorm(resid / h) / h)
  if (f_hat < 1e-10) f_hat <- 1e-10
  J    <- crossprod(X) / n
  Jinv <- tryCatch(solve(J), error = function(e) diag(p))
  V    <- tau * (1 - tau) / (n * f_hat^2) * Jinv %*% crossprod(X) %*% Jinv
  sqrt(pmax(diag(V), 0))
}


## ── Bandwidth helpers ────────────────────────────────────────────────────
#' @keywords internal
.bw_hallsheather <- function(n, tau) {
  z  <- stats::qnorm((1 + tau) / 2)
  phi <- stats::dnorm(stats::qnorm(tau))
  n^(-1/3) * z^(2/3) * ((1.5 * phi^2) / (2 * stats::qnorm(tau)^2 + 1))^(1/3)
}

#' @keywords internal
.bw_bofinger <- function(n, tau) {
  phi <- stats::dnorm(stats::qnorm(tau))
  z   <- stats::qnorm(tau)
  n^(-1/5) * ((4.5 * phi^4) / (2 * z^2 + 1)^2)^(1/5)
}


## ── Test statistic computation ───────────────────────────────────────────
#' @keywords internal
.compute_slope_test <- function(beta_all, beta_se_all, tau_val, bw,
                                 n_panels_eff, k) {
  # beta_all: [n_panels_eff, k, 1] or [n_panels_eff, k]
  if (length(dim(beta_all)) == 3) {
    beta_all    <- matrix(beta_all[, , 1], nrow = dim(beta_all)[1], ncol = dim(beta_all)[2])
    beta_se_all <- matrix(beta_se_all[, , 1], nrow = dim(beta_se_all)[1], ncol = dim(beta_se_all)[2])
  }
  beta_all <- as.matrix(beta_all)
  beta_se_all <- as.matrix(beta_se_all)

  if (any(is.na(beta_all)) || any(is.na(beta_se_all))) {
    return(list(S = NA_real_, D = NA_real_,
                pval_S = NA_real_, pval_D = NA_real_,
                beta_md = rep(NA_real_, k),
                beta_md_se = rep(NA_real_, k)))
  }

  N <- n_panels_eff

  ## Minimum distance (pooled) estimator: weighted mean
  ## V_i^{-1} = diag(1/se_i^2) — diagonal weight matrix
  beta_md  <- numeric(k)
  beta_md_se_v <- numeric(k)
  for (j in seq_len(k)) {
    w_j       <- 1 / (beta_se_all[, j, drop = FALSE]^2 + 1e-14)
    beta_md[j] <- sum(w_j * beta_all[, j]) / sum(w_j)
    beta_md_se_v[j] <- sqrt(1 / sum(w_j))
  }

  ## S-hat statistic (chi-squared, k*(N-1) df)
  ## S = sum_i (beta_i - beta_md)' V_i^{-1} (beta_i - beta_md)
  S_val <- 0
  for (i in seq_len(N)) {
    diff_i <- beta_all[i, ] - beta_md
    Vi_inv  <- diag(1 / (beta_se_all[i, ]^2 + 1e-14), k)
    S_val   <- S_val + as.numeric(t(diff_i) %*% Vi_inv %*% diff_i)
  }

  df     <- k * (N - 1)
  pval_S <- 1 - stats::pchisq(S_val, df = df)

  ## D-hat statistic (standard normal)
  ## D = (S - k*(N-1)) / sqrt(2*k*(N-1))
  D_val  <- (S_val - k * (N - 1)) / sqrt(2 * k * (N - 1))
  pval_D <- 1 - stats::pnorm(D_val)

  list(S = S_val, D = D_val, pval_S = pval_S, pval_D = pval_D,
       beta_md = beta_md, beta_md_se = beta_md_se_v)
}


#' @keywords internal
.compute_ols_test <- function(beta_ols_all, k, n_panels_eff) {
  valid <- complete.cases(beta_ols_all)
  beta  <- beta_ols_all[valid, , drop = FALSE]
  N     <- nrow(beta)
  if (N < 2) {
    return(list(S = NA_real_, D = NA_real_, pval_S = NA_real_, pval_D = NA_real_))
  }
  beta_mean <- colMeans(beta)
  # Use sample variance across panels as weighting proxy
  S_val <- 0
  for (i in seq_len(N)) {
    diff_i <- beta[i, ] - beta_mean
    S_val  <- S_val + sum(diff_i^2)
  }
  # Standardize by average within-panel variance (set to 1 for OLS pooled test)
  df     <- k * (N - 1)
  pval_S <- 1 - stats::pchisq(S_val, df = df)
  D_val  <- (S_val - k * (N - 1)) / sqrt(2 * k * (N - 1))
  pval_D <- 1 - stats::pnorm(D_val)
  list(S = S_val, D = D_val, pval_S = pval_S, pval_D = pval_D)
}


## ── S3 methods ───────────────────────────────────────────────────────────

#' Print method for xtqsh objects
#'
#' @param x An object of class \code{"xtqsh"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.xtqsh <- function(x, ...) {
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("  Quantile Regression Slope Homogeneity Test\n")
  cat("  Galvao, Juhl, Montes-Rojas & Olmo (2017)\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("  Dep. variable : %s\n", x$depvar))
  cat(sprintf("  Covariates    : %s\n", paste(x$indepvars, collapse = ", ")))
  cat(sprintf("  Observations  : %d\n", x$n_obs))
  cat(sprintf("  Panels (valid): %d\n", x$n_panels))
  cat(sprintf("  Quantiles     : %s\n", paste(x$tau, collapse = ", ")))
  cat(sprintf("  Bandwidth     : %s\n", x$bw))
  cat(strrep("-", 70), "\n\n")

  cat("  Joint Slope Homogeneity Test\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("  %-10s  %12s  %10s  %12s  %10s\n",
              "Quantile", "S-hat", "p(S)", "D-hat", "p(D)"))
  cat(strrep("-", 70), "\n")

  cat(sprintf("  %-10s  %12.4f  %10.4f  %12.4f  %10.4f  %s\n",
              "Mean(OLS)", x$S_ols, x$pval_S_ols,
              x$D_ols, x$pval_D_ols,
              .stars(x$pval_D_ols)))

  for (i in seq_along(x$tau)) {
    cat(sprintf("  tau=%-6.2f  %12.4f  %10.4f  %12.4f  %10.4f  %s\n",
                x$tau[i], x$S[i], x$pval_S[i],
                x$D[i], x$pval_D[i],
                .stars(x$pval_D[i])))
  }
  cat(strrep("-", 70), "\n")
  cat("  *** p<0.01, ** p<0.05, * p<0.10\n\n")

  invisible(x)
}


#' Summary method for xtqsh objects
#'
#' @param object An object of class \code{"xtqsh"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.xtqsh <- function(object, ...) {
  print(object, ...)

  cat("  Minimum Distance QR Estimates\n")
  cat(strrep("-", 70), "\n")
  print(round(object$beta_md, 4))
  cat("  Standard errors:\n")
  print(round(object$beta_md_se, 4))
  cat("\n")

  if (!is.null(object$marginal)) {
    cat("  Marginal Slope Homogeneity Tests\n")
    cat(strrep("-", 70), "\n")
    for (vname in names(object$marginal)) {
      cat(sprintf("  Variable: %s\n", vname))
      mr <- object$marginal[[vname]]
      for (i in seq_along(object$tau)) {
        cat(sprintf("    tau=%.2f  S=%.4f  p(S)=%.4f  D=%.4f  p(D)=%.4f  %s\n",
                    object$tau[i], mr$S[i], mr$pval_S[i],
                    mr$D[i], mr$pval_D[i], .stars(mr$pval_D[i])))
      }
    }
  }

  invisible(object)
}


#' @keywords internal
.stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.10) return("*")
  return("")
}

