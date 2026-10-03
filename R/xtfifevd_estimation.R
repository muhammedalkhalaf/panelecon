#' Internal estimation functions
#'
#' @name xtfifevd-estimation-internal
#' @md
#' @keywords internal
#' @noRd
#' @importFrom stats as.formula coef fitted lm.fit na.omit residuals var
#'   model.frame model.matrix terms complete.cases
NULL


#' Parse the xtfifevd formula
#'
#' Splits `y ~ x1 + x2 | z1 + z2` at the top-level `|` of the right-hand
#' side and builds the design matrices with [stats::model.frame()] and
#' [stats::model.matrix()], so that transformations (`log(x)`, `I(x^2)`,
#' `x:w`, `log(y)`) and factors (with their contrasts) are handled exactly as
#' in [stats::lm()]. The intercept column is removed from every part; the
#' intercept of the model is estimated in stage 2.
#'
#' @keywords internal
#' @noRd
.parse_formula <- function(formula, data, id, time, instruments = NULL,
                           na.action = na.omit) {

  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("Formula must be of the form: y ~ x1 + x2 | z1 + z2")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.character(id) || length(id) != 1L ||
      !is.character(time) || length(time) != 1L) {
    stop("'id' and 'time' must be single character strings naming columns of 'data'")
  }
  missing_idx <- setdiff(c(id, time), names(data))
  if (length(missing_idx) > 0) {
    stop("Variables not found in data: ", paste(missing_idx, collapse = ", "))
  }

  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  if (!(is.call(rhs) && identical(rhs[[1L]], as.name("|")))) {
    stop("Formula must include '|' to separate time-varying (x) and ",
         "time-invariant (z) variables.\n",
         "Example: y ~ x1 + x2 | z1 + z2")
  }
  x_part <- rhs[[2L]]
  z_part <- rhs[[3L]]
  has_bar <- function(e) {
    if (!is.call(e)) return(FALSE)
    if (identical(e[[1L]], as.name("|"))) return(TRUE)
    any(vapply(as.list(e)[-1L], has_bar, logical(1)))
  }
  if (has_bar(x_part) || has_bar(z_part)) {
    stop("Formula must contain exactly one '|' separating the time-varying ",
         "(x) and time-invariant (z) parts, e.g. y ~ x1 + x2 | z1 + z2")
  }

  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()

  y_formula <- as.formula(call("~", lhs), env = env)
  x_formula <- as.formula(call("~", x_part), env = env)
  z_formula <- as.formula(call("~", z_part), env = env)

  iv_formula <- NULL
  if (!is.null(instruments)) {
    if (!inherits(instruments, "formula") || length(instruments) != 2L) {
      stop("'instruments' must be a one-sided formula, e.g. ~ iv1 + iv2")
    }
    iv_formula <- instruments
  }

  # Check that all raw variables exist before touching model.frame so that
  # the error message names the offending variables.
  all_vars <- unique(c(all.vars(y_formula), all.vars(x_formula),
                       all.vars(z_formula),
                       if (!is.null(iv_formula)) all.vars(iv_formula)))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  n_all <- nrow(data)

  # Build design matrices on the full data (NA rows kept), then apply
  # na.action jointly so that NaN produced by a transformation (e.g. log of
  # a negative value) is treated like any other missing value.
  mm_part <- function(f) {
    mf <- model.frame(f, data = data, na.action = stats::na.pass)
    mm <- model.matrix(attr(mf, "terms"), mf)
    keep <- colnames(mm) != "(Intercept)"
    mm <- mm[, keep, drop = FALSE]
    attr(mm, "assign") <- NULL
    attr(mm, "contrasts") <- NULL
    mm
  }

  y_mf <- model.frame(y_formula, data = data, na.action = stats::na.pass)
  y <- y_mf[[1L]]
  if (!is.numeric(y)) stop("The response must be numeric")
  y <- as.numeric(y)
  y_name <- deparse(lhs, width.cutoff = 500L)

  X <- mm_part(x_formula)
  Z <- mm_part(z_formula)
  IV <- if (!is.null(iv_formula)) mm_part(iv_formula) else NULL

  if (ncol(X) == 0L) {
    stop("At least one time-varying regressor (before '|') is required")
  }
  if (ncol(Z) == 0L) {
    stop("At least one time-invariant regressor (after '|') is required")
  }

  panel_id <- data[[id]]
  time_id <- data[[time]]

  # Joint NA handling via na.action on a data frame of the used columns
  combined <- data.frame(y = y, X, Z, .id = panel_id, .time = time_id,
                         check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(IV)) combined <- cbind(combined, IV)
  rownames(combined) <- seq_len(n_all)
  combined_na <- na.action(combined)
  keep <- as.integer(rownames(combined_na))
  if (length(keep) == 0L) stop("No complete observations left after na.action")

  y <- y[keep]
  X <- X[keep, , drop = FALSE]
  Z <- Z[keep, , drop = FALSE]
  if (!is.null(IV)) IV <- IV[keep, , drop = FALSE]
  panel_id <- panel_id[keep]
  time_id <- time_id[keep]

  # Drop unused factor levels of the id (e.g. panels removed by na.action)
  if (is.factor(panel_id)) panel_id <- droplevels(panel_id)

  # Duplicated (id, time) pairs are not a valid panel
  if (anyDuplicated(data.frame(panel_id, time_id))) {
    stop("Duplicated (id, time) pairs found: each panel unit must be ",
         "observed at most once per period")
  }

  # Sort by id, time
  ord <- order(panel_id, time_id)
  y <- y[ord]
  X <- X[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  if (!is.null(IV)) IV <- IV[ord, , drop = FALSE]
  panel_id <- panel_id[ord]
  time_id <- time_id[ord]

  # Panel info
  panels <- unique(panel_id)
  group <- match(panel_id, panels)
  N_g <- length(panels)
  N <- length(y)
  panel_sizes <- tabulate(group, nbins = N_g)
  T_bar <- mean(panel_sizes)
  T_min <- min(panel_sizes)
  T_max <- max(panel_sizes)
  balanced <- (T_min == T_max)

  # Validate time-invariant variables (on the model-matrix columns).
  # Rarely changing variables are allowed with a warning; the package then
  # uses their unit means (see .unit_means).
  for (j in seq_len(ncol(Z))) {
    rng <- tapply(Z[, j], group, function(v) diff(range(v)))
    if (any(rng > 1e-8 * max(1, max(abs(Z[, j]))))) {
      warning("Variable '", colnames(Z)[j], "' is not strictly ",
              "time-invariant within panels; its unit (panel) mean is used ",
              "as the time-invariant regressor. Results may be unreliable.",
              call. = FALSE)
    }
  }

  list(
    y = y,
    X = X,
    Z = Z,
    IV = IV,
    panel_id = panel_id,
    time_id = time_id,
    group = group,
    y_name = y_name,
    x_names = colnames(X),
    z_names = colnames(Z),
    iv_names = if (!is.null(IV)) colnames(IV) else NULL,
    id_name = id,
    time_name = time,
    panels = panels,
    N = N,
    N_g = N_g,
    T_i = panel_sizes,
    T_bar = T_bar,
    T_min = T_min,
    T_max = T_max,
    balanced = balanced
  )
}


#' Unit (panel) means of the columns of a matrix
#' @keywords internal
#' @noRd
.unit_means <- function(M, group, N_g) {
  M <- as.matrix(M)
  out <- matrix(NA_real_, N_g, ncol(M))
  for (j in seq_len(ncol(M))) {
    out[, j] <- as.vector(tapply(M[, j], group, mean))
  }
  colnames(out) <- colnames(M)
  out
}


#' Fixed Effects estimation (Stage 1)
#'
#' Within (FE) regression of y on X. Returns the FE estimate of beta, the
#' within residuals e_it, the time averages of the FE residuals
#' ubar_i = ybar_i - xbar_i' beta (Pesaran and Zhou 2018, eq. 3), the
#' Arellano-type robust covariance of beta (eq. 18) and the classical
#' covariance sigma_e^2 (X'MX)^-1.
#'
#' @keywords internal
#' @noRd
.fe_stage1 <- function(y, X, group, N_g) {

  N <- length(y)
  k_x <- ncol(X)

  y_bar_i <- as.vector(tapply(y, group, mean))
  X_bar_i <- .unit_means(X, group, N_g)

  # Within transformation
  y_within <- y - y_bar_i[group]
  X_within <- X - X_bar_i[group, , drop = FALSE]

  # FE regression (demeaned)
  fit_fe <- lm.fit(X_within, y_within)
  beta <- coef(fit_fe)
  if (anyNA(beta)) {
    stop("The within design matrix is rank deficient (columns: ",
         paste(colnames(X)[is.na(beta)], collapse = ", "),
         "). Time-invariant variables cannot appear before '|'.")
  }
  names(beta) <- colnames(X)

  # Residuals
  e_it <- as.vector(residuals(fit_fe))       # within residuals
  u_it <- as.vector(y - X %*% beta)          # alpha + alpha_i + z'gamma + e_it

  # Time-averaged FE residuals ubar_i (PZ eq. 3, averaged over t)
  u_bar <- as.vector(y_bar_i - X_bar_i %*% beta)

  # Idiosyncratic variance
  df_residual <- N - N_g - k_x
  sigma2_e <- sum(e_it^2) / df_residual

  # Bread: (X'MX)^{-1}
  XtX_inv <- solve(crossprod(X_within))

  # Meat: sum over panels of (x_i' e_i)(x_i' e_i)' (PZ eq. 18, Arellano)
  scores <- matrix(0, N_g, k_x)
  for (j in seq_len(k_x)) {
    scores[, j] <- as.vector(tapply(X_within[, j] * e_it, group, sum))
  }
  meat <- crossprod(scores)

  V_beta_robust <- XtX_inv %*% meat %*% XtX_inv
  V_beta_classical <- sigma2_e * XtX_inv
  dimnames(V_beta_robust) <- dimnames(V_beta_classical) <-
    list(colnames(X), colnames(X))

  list(
    beta = beta,
    e_it = e_it,
    u_it = u_it,
    u_bar = u_bar,
    y_bar_i = y_bar_i,
    X_bar_i = X_bar_i,
    sigma2_e = sigma2_e,
    V_beta_robust = V_beta_robust,
    V_beta_classical = V_beta_classical,
    df_residual = df_residual
  )
}


#' FEF Stage 2 estimation (Pesaran and Zhou 2018, eq. 4 and 5)
#'
#' Unit-level OLS of ubar_i on an intercept and z_i. `Z_panel` holds one row
#' per panel unit (the unit means of the Z columns).
#'
#' @keywords internal
#' @noRd
.fef_stage2 <- function(u_bar, Z_panel) {

  Z_aug <- cbind(1, Z_panel)
  fit <- lm.fit(Z_aug, u_bar)

  coefs <- coef(fit)
  if (anyNA(coefs)) {
    stop("The stage 2 design matrix (intercept and time-invariant variables) ",
         "is rank deficient")
  }
  alpha <- unname(coefs[1L])
  gamma <- coefs[-1L]
  names(gamma) <- colnames(Z_panel)

  # Residuals: chat_i = ubar_i - alpha - z_i' gamma  (PZ eq. 20)
  chat <- as.vector(residuals(fit))

  list(
    gamma = gamma,
    alpha = alpha,
    chat = chat,
    Z_panel = Z_panel
  )
}


#' FEF-IV Stage 2 estimation (Pesaran and Zhou 2018, eq. 48)
#' @keywords internal
#' @noRd
.fef_iv_stage2 <- function(u_bar, Z_panel, IV_panel) {

  N_g <- length(u_bar)
  k_z <- ncol(Z_panel)
  k_iv <- ncol(IV_panel)

  if (k_iv < k_z) {
    stop("Number of instruments (", k_iv, ") must be >= number of ",
         "endogenous z-variables (", k_z, ")")
  }

  # Centred unit-level quantities
  Zc <- sweep(Z_panel, 2, colMeans(Z_panel))
  Rc <- sweep(IV_panel, 2, colMeans(IV_panel))
  uc <- u_bar - mean(u_bar)

  Qzr <- crossprod(Zc, Rc) / N_g
  Qrr <- crossprod(Rc) / N_g
  Qru <- crossprod(Rc, uc) / N_g

  # H_zr = (Q_zr Q_rr^{-1} Q_zr')^{-1} Q_zr Q_rr^{-1}
  Qrr_inv <- solve(Qrr)
  Hzr <- solve(Qzr %*% Qrr_inv %*% t(Qzr)) %*% Qzr %*% Qrr_inv

  gamma <- as.vector(Hzr %*% Qru)                     # eq. (48)
  names(gamma) <- colnames(Z_panel)
  alpha <- mean(u_bar) - sum(colMeans(Z_panel) * gamma)

  # IV residuals: xi_i = ubar_i - alpha - z_i' gamma_iv
  upsilon <- as.vector(u_bar - alpha - Z_panel %*% gamma)

  list(
    gamma = gamma,
    alpha = alpha,
    upsilon = upsilon,
    Z_panel = Z_panel,
    IV_panel = IV_panel,
    Hzr = Hzr
  )
}


#' Full Pesaran-Zhou covariance matrix of (beta, gamma, intercept)
#'
#' From PZ eq. (A.11) and (A.12),
#'   gamma_hat - gamma = G_z v0 - G_b (beta_hat - beta),
#' and from eq. (5), alpha_hat = ubar - zbar' gamma_hat, so
#'   alpha_hat - alpha = w' v0 / N - c' (beta_hat - beta),
#' where v0_i = alpha_i + epsbar_i (estimated by the stage 2 residual),
#' w_i = 1 - N * zbar transposed times column i of G_z and c = xbar - N^{-1} sum_i xbar_i w_i
#' (equivalently c = xbar - G_b' zbar). Stacking theta_2 = (gamma, alpha):
#'   Var(theta_2) = S diag(resid^2) S' + L Var(beta) L',
#'   Cov(theta_2, beta) = L Var(beta),
#' with S = rbind(G_z, w'/N) and L = rbind(-G_b, -c'). The gamma block is
#' exactly PZ eq. (17) (FEF) or eq. (51) (FEF-IV); the cross term between
#' the unit-level scores and beta_hat is dropped, as in eq. (17).
#'
#' @param G_z k_z x N matrix mapping the unit-level errors to gamma:
#'   Qzz^{-1} (z_i - zbar)/N for FEF, Hzr (r_i - rbar)/N for FEF-IV.
#' @param G_b k_z x k_x matrix: Qzz^{-1} Qzxbar (FEF) or Hzr Qrxbar (FEF-IV).
#' @keywords internal
#' @noRd
.pz_full_vcov <- function(G_z, G_b, zbar, xbar, resid, V_beta, names_all) {

  N_g <- length(resid)
  k_z <- nrow(G_z)
  k_x <- ncol(G_b)

  w <- as.vector(1 - crossprod(zbar, G_z) * N_g)     # length N_g
  cvec <- as.vector(xbar - crossprod(G_b, zbar))     # length k_x

  S <- rbind(G_z, w / N_g)                           # (k_z + 1) x N_g
  L <- rbind(-G_b, -cvec)                            # (k_z + 1) x k_x

  V22 <- S %*% (resid^2 * t(S)) + L %*% V_beta %*% t(L)
  V21 <- L %*% V_beta

  k <- k_x + k_z + 1
  V <- matrix(0, k, k)
  V[1:k_x, 1:k_x] <- V_beta
  V[(k_x + 1):k, (k_x + 1):k] <- V22
  V[(k_x + 1):k, 1:k_x] <- V21
  V[1:k_x, (k_x + 1):k] <- t(V21)
  V <- (V + t(V)) / 2
  dimnames(V) <- list(names_all, names_all)
  V
}


#' Pesaran-Zhou variance of gamma_FEF (eq. 17) and full covariance
#' @keywords internal
#' @noRd
.pz_variance_fef <- function(Z_panel, X_bar_i, chat, V_beta, N_g,
                             names_all) {

  Z_mean <- colMeans(Z_panel)
  X_mean <- colMeans(X_bar_i)
  Zc <- sweep(Z_panel, 2, Z_mean)
  Xc <- sweep(X_bar_i, 2, X_mean)

  Qzz <- crossprod(Zc) / N_g                 # eq. (8)
  Qzxbar <- crossprod(Zc, Xc) / N_g          # eq. (9)
  Qzz_inv <- solve(Qzz)

  G_z <- Qzz_inv %*% t(Zc) / N_g             # k_z x N
  G_b <- Qzz_inv %*% Qzxbar                  # k_z x k_x

  V <- .pz_full_vcov(G_z, G_b, Z_mean, X_mean, chat, V_beta, names_all)

  k_x <- ncol(X_bar_i)
  k_z <- ncol(Z_panel)
  V_gamma <- V[(k_x + 1):(k_x + k_z), (k_x + 1):(k_x + k_z), drop = FALSE]
  list(V = V, V_gamma = V_gamma)
}


#' Pesaran-Zhou variance of gamma_FEF-IV (eq. 51) and full covariance
#' @keywords internal
#' @noRd
.pz_variance_fef_iv <- function(Z_panel, IV_panel, X_bar_i, upsilon, Hzr,
                                V_beta, N_g, names_all) {

  Z_mean <- colMeans(Z_panel)
  R_mean <- colMeans(IV_panel)
  X_mean <- colMeans(X_bar_i)
  Rc <- sweep(IV_panel, 2, R_mean)
  Xc <- sweep(X_bar_i, 2, X_mean)

  Qrxbar <- crossprod(Rc, Xc) / N_g

  G_z <- Hzr %*% t(Rc) / N_g                 # k_z x N
  G_b <- Hzr %*% Qrxbar                      # k_z x k_x

  V <- .pz_full_vcov(G_z, G_b, Z_mean, X_mean, upsilon, V_beta, names_all)

  k_x <- ncol(X_bar_i)
  k_z <- ncol(Z_panel)
  V_gamma <- V[(k_x + 1):(k_x + k_z), (k_x + 1):(k_x + k_z), drop = FALSE]
  list(V = V, V_gamma = V_gamma)
}


#' Variance of the unexplained unit effect
#'
#' The stage 2 residual chat_i estimates alpha_i + epsbar_i (the part of the
#' unit effect not explained by z) plus the time average of the idiosyncratic
#' error, whose variance is sigma_e^2 / T_i. The unit-effect variance is
#' therefore estimated as the residual variance of stage 2 minus
#' sigma_e^2 mean(1 / T_i), truncated at zero.
#'
#' @keywords internal
#' @noRd
.sigma2_u <- function(resid, k_z, sigma2_e, T_i) {
  N_g <- length(resid)
  s2 <- sum(resid^2) / max(N_g - k_z - 1, 1)
  max(s2 - sigma2_e * mean(1 / T_i), 0)
}


#' Shared assembly of the result list for FEF and FEF-IV
#' @keywords internal
#' @noRd
.assemble <- function(parsed, stage1, stage2, pz, Z_panel, vcov_beta) {

  k_x <- length(stage1$beta)
  k_z <- length(stage2$gamma)

  coefficients <- c(stage1$beta, stage2$gamma, "_cons" = stage2$alpha)

  resid2 <- if (!is.null(stage2$chat)) stage2$chat else stage2$upsilon

  # Fitted values use the unit-mean z (expanded to observations)
  Z_it <- Z_panel[parsed$group, , drop = FALSE]
  fitted <- as.vector(parsed$X %*% stage1$beta + Z_it %*% stage2$gamma +
                        stage2$alpha)

  list(
    coefficients = coefficients,
    vcov = pz$V,
    beta = stage1$beta,
    gamma = stage2$gamma,
    intercept = stage2$alpha,
    residuals = stage1$e_it,
    fitted.values = fitted,
    sigma2_e = stage1$sigma2_e,
    sigma2_u = .sigma2_u(resid2, k_z, stage1$sigma2_e, parsed$T_i),
    N = parsed$N,
    N_g = parsed$N_g,
    T_bar = parsed$T_bar,
    balanced = parsed$balanced,
    k_x = k_x,
    k_z = k_z,
    x_names = parsed$x_names,
    z_names = parsed$z_names,
    y_name = parsed$y_name,
    vcov_beta = vcov_beta,
    V_beta_robust = stage1$V_beta_robust,
    V_beta_classical = stage1$V_beta_classical,
    V_gamma_pz = pz$V_gamma,
    stage2_residuals = resid2,
    Z_panel = Z_panel,
    X_bar_i = stage1$X_bar_i,
    u_bar = stage1$u_bar
  )
}


#' FEF estimator (Pesaran and Zhou 2018, Section 3.1)
#' @keywords internal
#' @noRd
.estimate_fef <- function(parsed, vcov_beta = c("robust", "classical")) {

  vcov_beta <- match.arg(vcov_beta)

  # Stage 1: Fixed effects
  stage1 <- .fe_stage1(parsed$y, parsed$X, parsed$group, parsed$N_g)
  V_beta <- if (vcov_beta == "robust") stage1$V_beta_robust else
    stage1$V_beta_classical

  # Unit-level z (unit means; equal to z_i when z is time-invariant)
  Z_panel <- .unit_means(parsed$Z, parsed$group, parsed$N_g)

  # Stage 2: FEF regression
  stage2 <- .fef_stage2(stage1$u_bar, Z_panel)

  names_all <- c(names(stage1$beta), names(stage2$gamma), "_cons")

  # Pesaran-Zhou variance (eq. 17) and full covariance (eq. A.11 and 5)
  pz <- .pz_variance_fef(Z_panel, stage1$X_bar_i, stage2$chat, V_beta,
                         parsed$N_g, names_all)

  .assemble(parsed, stage1, stage2, pz, Z_panel, vcov_beta)
}


#' FEVD estimator (Plumper and Troeger 2007, three stages)
#'
#' Stage 1: within FE regression. Stage 2: unit-level regression of the
#' time-averaged FE residuals on an intercept and z (h_i are the residuals).
#' Stage 3: pooled OLS of y on an intercept, x, z and h_i (PT eq. 7).
#' Standard errors for beta, gamma and the intercept are the Pesaran and Zhou
#' (2018) ones; the naive stage 3 OLS standard errors are returned in
#' `stage3` for reference only.
#'
#' @keywords internal
#' @noRd
.estimate_fevd <- function(parsed, vcov_beta = c("robust", "classical")) {

  vcov_beta <- match.arg(vcov_beta)

  # Stages 1 and 2 and the Pesaran-Zhou covariance are those of FEF
  res <- .estimate_fef(parsed, vcov_beta = vcov_beta)

  # Stage 2 residuals h_i (PT eq. 6, with the intercept of PZ eq. 21),
  # expanded to the observation level
  h_i <- res$stage2_residuals
  h_it <- h_i[parsed$group]
  Z_it <- res$Z_panel[parsed$group, , drop = FALSE]

  # Stage 3: pooled OLS of y on 1, x, z, h (PT eq. 7)
  W <- cbind("(Intercept)" = 1, parsed$X, Z_it, h = h_it)
  fit3 <- lm.fit(W, parsed$y)
  coef3 <- coef(fit3)
  names(coef3) <- colnames(W)
  e3 <- as.vector(residuals(fit3))
  df3 <- parsed$N - ncol(W)
  s2_3 <- sum(e3^2) / df3
  ok <- !is.na(coef3)
  V3 <- matrix(NA_real_, ncol(W), ncol(W), dimnames = list(colnames(W), colnames(W)))
  V3[ok, ok] <- s2_3 * solve(crossprod(W[, ok, drop = FALSE]))

  k_x <- res$k_x
  k_z <- res$k_z
  beta3 <- coef3[2:(1 + k_x)]
  gamma3 <- coef3[(2 + k_x):(1 + k_x + k_z)]
  cons3 <- unname(coef3[1L])
  delta <- unname(coef3["h"])

  # Reported coefficients: stage 3 values. They are identical to FEF and
  # delta = 1 (Pesaran and Zhou 2018, Proposition 3): the within residuals
  # e_it sum to zero within each unit, so (a, beta_FE, gamma_FEF, 1) solves
  # the stage 3 normal equations exactly, also in unbalanced panels.
  res$fef <- list(coefficients = res$coefficients, beta = res$beta,
                  gamma = res$gamma, intercept = res$intercept)
  res$coefficients <- c(beta3, gamma3, "_cons" = cons3)
  res$beta <- beta3
  res$gamma <- gamma3
  res$intercept <- cons3
  res$delta <- delta
  res$fitted.values <- as.vector(parsed$X %*% beta3 + Z_it %*% gamma3 + cons3)
  res$stage3 <- list(
    coefficients = coef3,
    vcov_naive = V3,
    se_naive = sqrt(diag(V3)),
    sigma2 = s2_3,
    df_residual = df3
  )
  res$h_i <- h_i
  res
}


#' FEF-IV estimator (Pesaran and Zhou 2018, Section 4.2)
#' @keywords internal
#' @noRd
.estimate_fef_iv <- function(parsed, vcov_beta = c("robust", "classical")) {

  vcov_beta <- match.arg(vcov_beta)

  if (is.null(parsed$IV)) {
    stop("FEF-IV requires instruments. Provide 'instruments' argument.")
  }

  # Stage 1: Fixed effects
  stage1 <- .fe_stage1(parsed$y, parsed$X, parsed$group, parsed$N_g)
  V_beta <- if (vcov_beta == "robust") stage1$V_beta_robust else
    stage1$V_beta_classical

  Z_panel <- .unit_means(parsed$Z, parsed$group, parsed$N_g)
  IV_panel <- .unit_means(parsed$IV, parsed$group, parsed$N_g)

  # Stage 2: FEF-IV regression (eq. 48)
  stage2 <- .fef_iv_stage2(stage1$u_bar, Z_panel, IV_panel)

  names_all <- c(names(stage1$beta), names(stage2$gamma), "_cons")

  # Pesaran-Zhou variance for IV (eq. 51) and full covariance
  pz <- .pz_variance_fef_iv(Z_panel, IV_panel, stage1$X_bar_i, stage2$upsilon,
                            stage2$Hzr, V_beta, parsed$N_g, names_all)

  res <- .assemble(parsed, stage1, stage2, pz, Z_panel, vcov_beta)
  res$k_iv <- ncol(parsed$IV)
  res$iv_names <- parsed$iv_names
  res$IV_panel <- IV_panel
  res
}
