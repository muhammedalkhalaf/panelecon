#' Methods for xtfifevd objects
#'
#' @md
#' @name xtfifevd-methods
#' @description
#' S3 methods for objects of class `"xtfifevd"` returned by [xtfifevd()],
#' [fevd()], [fef()] and [fef_iv()]. `summary()` reports the coefficients
#' with the Pesaran and Zhou (2018) standard errors taken from the full
#' covariance matrix; `coef()`, `vcov()` and `confint()` return the
#' coefficient vector, the full covariance matrix and normal-based
#' confidence intervals.
#'
#' @param x An object of class `"xtfifevd"` (or `"summary.xtfifevd"` for
#'   the print method of the summary).
#' @param object An object of class `"xtfifevd"`.
#' @param digits Number of significant digits to print.
#' @param signif.stars Logical; whether to print significance stars.
#' @param parm Coefficients to include in the confidence intervals
#'   (names or indices); all coefficients by default.
#' @param level Confidence level. Default is 0.95.
#' @param ... Further arguments passed to [stats::printCoefmat()] by the
#'   summary print method and ignored otherwise.
#'
#' @return `print.xtfifevd()` and `print.summary.xtfifevd()` return their
#'   first argument invisibly. `summary.xtfifevd()` returns an object of
#'   class `"summary.xtfifevd"`, a list with the coefficient table and the
#'   model information. `coef.xtfifevd()` returns the named coefficient
#'   vector, `vcov.xtfifevd()` the full covariance matrix and
#'   `confint.xtfifevd()` a matrix of confidence limits.
#'
#' @seealso [xtfifevd()]
#' @importFrom stats coef vcov pnorm qnorm printCoefmat confint
#' @importFrom utils packageVersion
NULL


#' @rdname xtfifevd-methods
#' @export
print.xtfifevd <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  
  cat("\n")
  cat(x$method, "Estimation Results\n")
  cat(rep("-", 50), sep = "")
  cat("\n\n")
  
  cat("Dep. variable:", x$y_name, "\n")
  cat("Method:       ", x$method, "\n")
  cat("Observations: ", x$N, "\n")
  cat("Groups:       ", x$N_g, "\n")
  cat("T (average):  ", round(x$T_bar, 2), "\n")
  cat("\n")
  
  cat("Coefficients:\n")
  print(round(x$coefficients, digits))
  if (!is.null(x$delta)) {
    cat(sprintf("Stage 3 coefficient on h_i (delta): %s\n",
                format(round(x$delta, digits))))
  }
  cat("\n")
  
  invisible(x)
}


#' @rdname xtfifevd-methods
#' @export
summary.xtfifevd <- function(object, ...) {
  
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z_val <- coefs / se
  p_val <- 2 * pnorm(-abs(z_val))
  
  coef_table <- cbind(
    Estimate = coefs,
    `Std. Error` = se,
    `z value` = z_val,
    `Pr(>|z|)` = p_val
  )
  
  # Add significance stars
  stars <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*",
                                ifelse(p_val < 0.1, ".", ""))))
  
  result <- list(
    call = object$call,
    method = object$method,
    coefficients = coef_table,
    stars = stars,
    N = object$N,
    N_g = object$N_g,
    T_bar = object$T_bar,
    sigma2_e = object$sigma2_e,
    sigma2_u = object$sigma2_u,
    vcov_beta = object$vcov_beta,
    delta = object$delta,
    balanced = object$balanced,
    y_name = object$y_name,
    x_names = object$x_names,
    z_names = object$z_names,
    k_x = object$k_x,
    k_z = object$k_z
  )
  
  class(result) <- "summary.xtfifevd"
  result
}


#' @rdname xtfifevd-methods
#' @export
print.summary.xtfifevd <- function(x, digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"),
                                   ...) {
  
  cat("\n")
  cat(rep("=", 70), sep = "")
  cat("\n")
  ver <- tryCatch(as.character(utils::packageVersion("panelecon")),
                  error = function(e) "")
  cat(sprintf("%-50s panelecon %s\n", paste(x$method, "Estimation Results"), ver))
  cat(rep("=", 70), sep = "")
  cat("\n")
  
  cat("Dep. variable:  ", x$y_name, "\n")
  cat("Method:         ", x$method, "\n")
  cat("Variance:        Pesaran and Zhou (2018), beta vcov:", x$vcov_beta, "\n")
  cat(sprintf("Observations:    %-10d Groups:     %d\n", x$N, x$N_g))
  cat(sprintf("T (average):     %-10.2f\n", x$T_bar))
  cat(rep("-", 70), sep = "")
  cat("\n\n")
  
  # Print coefficients
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               P.values = TRUE, has.Pvalue = TRUE, ...)
  
  cat("\n")
  cat(rep("-", 70), sep = "")
  cat("\n")
  cat("Time-varying (FE):     ", paste(x$x_names, collapse = ", "), "\n")
  cat("Time-invariant:        ", paste(x$z_names, collapse = ", "), "\n")
  cat(sprintf("sigma_e: %.4f    sigma_u (unexplained unit effect): %.4f\n",
              sqrt(x$sigma2_e), sqrt(x$sigma2_u)))
  if (!is.null(x$delta)) {
    cat(sprintf("FEVD stage 3 coefficient on h_i (delta): %.6f", x$delta))
    cat("  [equals 1 by construction]\n")
    cat("Naive stage 3 OLS SEs are too small for the time-invariant\n",
        "coefficients; see ?xtfifevd.\n", sep = "")
  }
  cat(rep("=", 70), sep = "")
  cat("\n")
  
  invisible(x)
}


#' @rdname xtfifevd-methods
#' @export
coef.xtfifevd <- function(object, ...) {
  object$coefficients
}


#' @rdname xtfifevd-methods
#' @export
vcov.xtfifevd <- function(object, ...) {
  object$vcov
}


#' @rdname xtfifevd-methods
#' @export
confint.xtfifevd <- function(object, parm, level = 0.95, ...) {
  
  cf <- coef(object)
  se <- sqrt(diag(vcov(object)))
  
  if (missing(parm)) {
    parm <- names(cf)
  } else if (is.numeric(parm)) {
    parm <- names(cf)[parm]
  }
  
  a <- (1 - level) / 2
  fac <- qnorm(1 - a)
  
  ci <- cbind(cf[parm] - fac * se[parm],
              cf[parm] + fac * se[parm])
  
  pct <- paste0(format(100 * c(a, 1 - a), digits = 3), "%")
  colnames(ci) <- pct
  rownames(ci) <- parm
  
  ci
}
