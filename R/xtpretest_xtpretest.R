#' Comprehensive Panel Data Pre-Testing Suite
#'
#' Performs a full battery of panel data pre-tests: Hsiao (2014) homogeneity
#' F-tests, robust (HC1) versions, Swamy (1970) parameter heterogeneity test,
#' cross-sectional dependence (Pesaran 2004), and panel summary statistics.
#'
#' @param data A \code{data.frame} in long format containing all variables.
#' @param formula A two-sided formula of the form \code{y ~ x1 + x2 + ...}
#'   specifying the dependent and independent variables.
#' @param index Character vector of length 2: \code{c("panel_id", "time_id")}.
#' @param tests Character vector. Which modules to run. Possible values:
#'   \code{"summary"}, \code{"hsiao"}, \code{"robust"}, \code{"heterogeneity"},
#'   \code{"csd"}.  Default \code{"ALL"} runs everything.
#' @param level Numeric. Significance level for decisions (default 0.05).
#'
#' @return A list (invisibly) with components:
#'   \describe{
#'     \item{\code{summary}}{Panel summary statistics.}
#'     \item{\code{hsiao}}{Hsiao homogeneity F-test results.}
#'     \item{\code{robust}}{Robust HC1 F-test results.}
#'     \item{\code{swamy}}{Swamy heterogeneity test results.}
#'     \item{\code{csd}}{Cross-sectional dependence test results.}
#'     \item{\code{recommendation}}{Character. Suggested estimator.}
#'   }
#'
#' @references
#' Hsiao, C. (2014). \emph{Analysis of Panel Data} (3rd ed.).
#' Cambridge University Press.
#' \doi{10.1017/CBO9781139839327}
#'
#' Swamy, P. A. V. B. (1970). Efficient inference in a random coefficient
#' regression model. \emph{Econometrica}, 38(2), 311-323.
#' \doi{10.2307/1909405}
#'
#' Pesaran, M. H. (2004). General diagnostic tests for cross section
#' dependence in panels. \emph{Cambridge Working Paper in Economics}, No. 0435.
#' \doi{10.2139/ssrn.572504}
#'
#' @examples
#' set.seed(10)
#' n <- 5; t <- 10
#' df <- data.frame(
#'   id   = rep(1:n, each = t),
#'   time = rep(1:t, times = n),
#'   y    = rnorm(n * t),
#'   x1   = rnorm(n * t)
#' )
#' res <- xtpretest(df, y ~ x1, index = c("id", "time"),
#'                  tests = c("hsiao", "csd"))
#'
#' @export
xtpretest <- function(data,
                      formula,
                      index,
                      tests = "ALL",
                      level = 0.05) {

  ## ---- Validate --------------------------------------------------------
  if (!inherits(data, "data.frame")) stop("'data' must be a data.frame.")
  if (length(index) != 2L)          stop("'index' must be a character vector of length 2.")
  pid_col <- index[1L]; tid_col <- index[2L]
  if (!all(c(pid_col, tid_col) %in% names(data)))
    stop("index columns not found in 'data'.")

  tests <- toupper(tests)
  if (identical(tests, "ALL")) {
    tests <- c("SUMMARY", "HSIAO", "ROBUST", "HETEROGENEITY", "CSD")
  }

  ## ---- Extract variables -----------------------------------------------
  mf   <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  lhs  <- names(mf)[1L]
  rhs  <- names(mf)[-1L]
  y    <- data[[lhs]]
  X    <- data[, rhs, drop = FALSE]
  pid  <- data[[pid_col]]
  tid  <- data[[tid_col]]

  panels <- sort(unique(pid))
  N      <- length(panels)
  T_avg  <- nrow(data) / N
  K      <- length(rhs)

  .xtp_header(lhs, rhs, pid_col, tid_col, N, T_avg, nrow(data))

  out <- list()

  ## ---- MODULE: Summary ------------------------------------------------
  if ("SUMMARY" %in% tests) {
    out$summary <- .xtp_summary(data, formula, pid_col, tid_col, level)
  }

  ## ---- MODULE: Hsiao --------------------------------------------------
  if ("HSIAO" %in% tests) {
    out$hsiao <- .xtp_hsiao(y, X, pid, panels, K, level)
  }

  ## ---- MODULE: Robust -------------------------------------------------
  if ("ROBUST" %in% tests) {
    out$robust <- .xtp_robust(y, X, pid, data, pid_col, rhs, level)
  }

  ## ---- MODULE: Heterogeneity ------------------------------------------
  if ("HETEROGENEITY" %in% tests) {
    out$swamy <- .xtp_swamy(y, X, pid, panels, K)
  }

  ## ---- MODULE: CSD ----------------------------------------------------
  if ("CSD" %in% tests) {
    out$csd <- .xtp_csd_pesaran(y, X, pid, tid, panels, level)
  }

  ## ---- Recommendation -------------------------------------------------
  out$recommendation <- .xtp_recommend(out, level)
  .xtp_print_summary(out, level)

  invisible(out)
}


## ==========================================================================
## Internal helpers
## ==========================================================================

#' @keywords internal
.xtp_header <- function(lhs, rhs, pid, tid, N, T_avg, NT) {
  line <- strrep("-", 72)
  message(line)
  message("  xtpretest -- Comprehensive Panel Data Pre-Testing Suite v1.0.0")
  message(line)
  message("  Dependent  : ", lhs)
  message("  Regressors : ", paste(rhs, collapse = ", "))
  message("  Panel id   : ", pid, "  |  Time id: ", tid)
  message("  N panels   : ", N, "  |  Avg T: ", round(T_avg, 1),
          "  |  Obs: ", NT)
  message(line)
}

#' @keywords internal
.xtp_summary <- function(data, formula, pid_col, tid_col, level) {
  mf   <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  vars <- names(mf)
  pid  <- data[[pid_col]]
  out  <- list()

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Panel Summary Statistics")
  message(line)
  message(sprintf("  %-16s %-10s %-10s %-10s %-10s %-10s",
                  "Variable", "Overall SD", "Between SD", "Within SD", "Min", "Max"))
  message(line)

  for (v in vars) {
    x       <- data[[v]]
    o_sd    <- stats::sd(x, na.rm = TRUE)
    pmeans  <- tapply(x, pid, mean, na.rm = TRUE)
    b_sd    <- stats::sd(pmeans, na.rm = TRUE)
    grand   <- mean(x, na.rm = TRUE)
    within  <- x - pmeans[match(pid, names(pmeans))] + grand
    w_sd    <- stats::sd(within, na.rm = TRUE)
    message(sprintf("  %-16s %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f",
                    utils::head(v, 1L),
                    o_sd, b_sd, w_sd,
                    min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
    out[[v]] <- list(overall_sd = o_sd, between_sd = b_sd, within_sd = w_sd)
  }
  message(line)
  out
}

#' @keywords internal
.xtp_hsiao <- function(y, X, pid, panels, K, level) {
  N <- length(panels)
  ## S1: unrestricted (individual OLS)
  S1 <- 0; df1 <- 0
  for (p in panels) {
    idx  <- pid == p
    if (sum(idx) < K + 2L) next
    yi   <- y[idx]; Xi <- X[idx, , drop = FALSE]
    mat  <- cbind(1, as.matrix(Xi))
    fit  <- tryCatch(stats::lm.fit(mat, yi), error = function(e) NULL)
    if (!is.null(fit)) {
      S1  <- S1  + sum(fit$residuals^2)
      df1 <- df1 + (sum(idx) - K - 1L)
    }
  }

  ## S2: FE (common slopes, different intercepts)
  fe_mat <- cbind(stats::model.matrix(~ factor(pid) - 1), as.matrix(X))
  fit2   <- tryCatch(stats::lm.fit(fe_mat, y), error = function(e) NULL)
  S2     <- if (!is.null(fit2)) sum(fit2$residuals^2) else NA_real_
  df2    <- if (!is.null(fit2)) length(y) - ncol(fe_mat) else NA_integer_

  ## S3: pooled OLS
  mat3 <- cbind(1, as.matrix(X))
  fit3 <- tryCatch(stats::lm.fit(mat3, y), error = function(e) NULL)
  S3   <- if (!is.null(fit3)) sum(fit3$residuals^2) else NA_real_
  df3  <- if (!is.null(fit3)) length(y) - ncol(mat3) else NA_integer_

  ## F-statistics
  F1_nd <- (N - 1L) * (K + 1L); F1_dd <- df1
  F1  <- if (df1 > 0L && !is.na(S3)) ((S3 - S1) / F1_nd) / (S1 / F1_dd) else NA_real_
  F1p <- if (!is.na(F1)) stats::pf(F1, F1_nd, F1_dd, lower.tail = FALSE) else NA_real_

  F2_nd <- (N - 1L) * K; F2_dd <- df1
  F2  <- if (df1 > 0L && !is.na(S2)) ((S2 - S1) / F2_nd) / (S1 / F2_dd) else NA_real_
  F2p <- if (!is.na(F2)) stats::pf(F2, F2_nd, F2_dd, lower.tail = FALSE) else NA_real_

  F3_nd <- N - 1L; F3_dd <- df2
  F3  <- if (!is.na(df2) && df2 > 0L && !is.na(S3)) ((S3 - S2) / F3_nd) / (S2 / F3_dd) else NA_real_
  F3p <- if (!is.na(F3)) stats::pf(F3, F3_nd, F3_dd, lower.tail = FALSE) else NA_real_

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Hsiao (2014) Homogeneity Tests")
  message(line)
  message("  RSS Unrestricted (S1): ", formatC(S1, digits = 4, format = "f"))
  message("  RSS FE (S2):           ", formatC(S2, digits = 4, format = "f"))
  message("  RSS Pooled (S3):       ", formatC(S3, digits = 4, format = "f"))
  message(line)
  .xtp_print_ftest("F1: Overall homogeneity",   F1, F1_nd, F1_dd, F1p, level)
  .xtp_print_ftest("F2: Slope homogeneity",     F2, F2_nd, F2_dd, F2p, level)
  .xtp_print_ftest("F3: Intercept homogeneity", F3, F3_nd, F3_dd, F3p, level)
  message(line)

  list(F1 = F1, F1p = F1p, F1_nd = F1_nd, F1_dd = F1_dd,
       F2 = F2, F2p = F2p, F2_nd = F2_nd, F2_dd = F2_dd,
       F3 = F3, F3p = F3p, F3_nd = F3_nd, F3_dd = F3_dd,
       S1 = S1, S2 = S2, S3 = S3)
}

#' @keywords internal
.xtp_print_ftest <- function(label, F, df1, df2, pval, level) {
  sig <- if (!is.na(pval) && pval < level) " [REJECT]" else " [accept]"
  message(sprintf("  %-38s F=%8.4f  df=(%d,%d)  p=%6.4f%s",
                  label,
                  ifelse(is.na(F), 0, F),
                  df1, df2,
                  ifelse(is.na(pval), 1, pval),
                  sig))
}

#' @keywords internal
.xtp_robust <- function(y, X, pid, data, pid_col, rhs, level) {
  ## HC1 robust Wald tests via OLS with sandwich SE
  fe_dummy <- stats::model.matrix(~ factor(pid) - 1)
  int_terms <- do.call(cbind, lapply(rhs, function(v) {
    fe_dummy * data[[v]]
  }))
  colnames(int_terms) <- paste0("int_", seq_len(ncol(int_terms)))

  ## Full model
  full_mat <- cbind(1, as.matrix(X), fe_dummy[, -1L, drop = FALSE], int_terms)
  fit_full <- tryCatch(stats::lm.fit(full_mat, y), error = function(e) NULL)
  if (is.null(fit_full)) {
    message("  MODULE: Robust tests skipped (collinearity).")
    return(NULL)
  }

  ## HC1 sandwich
  hc1 <- .xtp_hc1(full_mat, fit_full$residuals)

  ## Robust F (Wald) for interaction terms (slope homogeneity)
  K_int <- ncol(int_terms)
  K_fe  <- ncol(fe_dummy) - 1L

  rF2 <- .xtp_wald(fit_full, hc1, seq(ncol(X) + 2L, ncol(X) + 1L + K_fe + K_int, 1L)[(K_fe + 1L):(K_fe + K_int)])
  rF3 <- .xtp_wald(fit_full, hc1, seq(ncol(X) + 2L, ncol(X) + 1L + K_fe, 1L))

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Robust (HC1) Homogeneity Tests")
  message(line)
  if (!is.null(rF2)) .xtp_print_ftest("Robust F2: Slope homogeneity", rF2$F, rF2$df1, rF2$df2, rF2$p, level)
  if (!is.null(rF3)) .xtp_print_ftest("Robust F3: Intercept homogeneity", rF3$F, rF3$df1, rF3$df2, rF3$p, level)
  message(line)

  list(rF2 = rF2, rF3 = rF3)
}

#' @keywords internal
.xtp_hc1 <- function(X, resid) {
  n <- nrow(X); k <- ncol(X)
  e2 <- resid^2 * (n / (n - k))
  Xe <- X * e2
  solve(crossprod(X)) %*% crossprod(Xe) %*% solve(crossprod(X))
}

#' @keywords internal
.xtp_wald <- function(fit, vcov, idx) {
  b   <- fit$coefficients[idx]
  V   <- vcov[idx, idx, drop = FALSE]
  df1 <- length(idx)
  df2 <- length(fit$residuals) - length(fit$coefficients)
  Wst <- tryCatch(as.numeric(t(b) %*% solve(V) %*% b) / df1,
                  error = function(e) NA_real_)
  pval <- if (!is.na(Wst)) stats::pf(Wst, df1, df2, lower.tail = FALSE) else NA_real_
  list(F = Wst, df1 = df1, df2 = df2, p = pval)
}

#' @keywords internal
.xtp_swamy <- function(y, X, pid, panels, K) {
  N <- length(panels)
  ## FE estimator
  fe_mat <- cbind(stats::model.matrix(~ factor(pid) - 1), as.matrix(X))
  fit_fe <- tryCatch(stats::lm.fit(fe_mat, y), error = function(e) NULL)
  if (is.null(fit_fe)) return(NULL)
  b_fe   <- utils::tail(fit_fe$coefficients, K)

  chi2 <- 0
  for (p in panels) {
    idx <- pid == p
    if (sum(idx) < K + 2L) next
    yi  <- y[idx]; Xi <- as.matrix(X[idx, , drop = FALSE])
    mat <- cbind(1, Xi)
    fit <- tryCatch(stats::lm(yi ~ Xi), error = function(e) NULL)
    if (is.null(fit)) next
    bi  <- stats::coef(fit)[-1L]
    Vi  <- stats::vcov(fit)[-1L, -1L, drop = FALSE]
    d   <- bi - b_fe
    chi2 <- chi2 + tryCatch(as.numeric(t(d) %*% solve(Vi) %*% d),
                             error = function(e) 0)
  }
  df   <- K * (N - 1L)
  pval <- stats::pchisq(chi2, df = df, lower.tail = FALSE)

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Swamy (1970) Parameter Heterogeneity")
  message(line)
  message(sprintf("  chi2(%d) = %.4f   p-value = %.4f  %s",
                  df, chi2, pval,
                  if (pval < 0.05) "[REJECT: heterogeneity]" else "[accept: homogeneity]"))
  message(line)

  list(chi2 = chi2, df = df, pval = pval)
}

#' @keywords internal
.xtp_csd_pesaran <- function(y, X, pid, tid, panels, level) {
  ## Run FE on the full model to get residuals
  fe_mat <- cbind(stats::model.matrix(~ factor(pid) - 1), as.matrix(X))
  fit_fe <- tryCatch(stats::lm.fit(fe_mat, y), error = function(e) NULL)
  if (is.null(fit_fe)) {
    message("  MODULE: CSD skipped.")
    return(NULL)
  }
  resid <- fit_fe$residuals
  N     <- length(panels)
  tvals <- sort(unique(tid))
  T_len <- length(tvals)

  ## Build residual matrix T x N
  R <- matrix(NA_real_, nrow = T_len, ncol = N)
  for (j in seq_along(panels)) {
    idx <- pid == panels[j]
    ti  <- match(tid[idx], tvals)
    R[ti, j] <- resid[idx]
  }

  ## Pairwise correlations
  cd_sum  <- 0; rho_sum <- 0; npairs <- 0
  for (i in seq_len(N - 1L)) {
    for (j in (i + 1L):N) {
      valid <- !is.na(R[, i]) & !is.na(R[, j])
      Tij   <- sum(valid)
      if (Tij < 3L) next
      rij   <- stats::cor(R[valid, i], R[valid, j])
      cd_sum  <- cd_sum  + sqrt(Tij) * rij
      rho_sum <- rho_sum + abs(rij)
      npairs  <- npairs + 1L
    }
  }
  cd_stat <- sqrt(2 / (N * (N - 1L))) * cd_sum
  pval    <- 2 * stats::pnorm(abs(cd_stat), lower.tail = FALSE)
  avg_rho <- rho_sum / max(npairs, 1L)

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Cross-Sectional Dependence (Pesaran 2004)")
  message(line)
  message(sprintf("  CD statistic   = %.4f", cd_stat))
  message(sprintf("  p-value        = %.4f", pval))
  message(sprintf("  Avg |rho_ij|   = %.4f", avg_rho))
  if (pval < level) {
    message("  -> Cross-sectional dependence detected. Consider CCE/CCEMG.")
  } else {
    message("  -> No significant cross-sectional dependence.")
  }
  message(line)

  list(cd_stat = cd_stat, pval = pval, avg_rho = avg_rho)
}

#' @keywords internal
.xtp_recommend <- function(out, level) {
  rec <- "Pooled OLS"

  if (!is.null(out$hsiao)) {
    if (!is.na(out$hsiao$F2p) && out$hsiao$F2p < level) rec <- "Mean Group (MG) / CCEMG"
    else if (!is.na(out$hsiao$F3p) && out$hsiao$F3p < level) rec <- "Fixed Effects (FE)"
  }
  if (!is.null(out$swamy)) {
    if (!is.na(out$swamy$pval) && out$swamy$pval < level && rec == "Pooled OLS")
      rec <- "Fixed Effects (FE)"
  }
  if (!is.null(out$csd)) {
    if (!is.na(out$csd$pval) && out$csd$pval < level) {
      if (rec == "Mean Group (MG) / CCEMG") rec <- "CCEMG (cross-sectional dependence)"
      else if (rec == "Fixed Effects (FE)") rec <- "FE with Driscoll-Kraay SE"
    }
  }
  rec
}

#' @keywords internal
.xtp_print_summary <- function(out, level) {
  line <- strrep("-", 72)
  message(line)
  message("  OVERALL PRE-TEST SUMMARY")
  message(line)
  if (!is.null(out$hsiao)) {
    message("  Hsiao F1 (overall)   : p=", round(out$hsiao$F1p, 3L))
    message("  Hsiao F2 (slopes)    : p=", round(out$hsiao$F2p, 3L))
    message("  Hsiao F3 (intercepts): p=", round(out$hsiao$F3p, 3L))
  }
  if (!is.null(out$swamy)) {
    message("  Swamy test           : p=", round(out$swamy$pval, 3L))
  }
  if (!is.null(out$csd)) {
    message("  CSD (Pesaran)        : p=", round(out$csd$pval, 3L))
  }
  message("  >>> Recommended estimator: ", out$recommendation)
  message(line)
}
