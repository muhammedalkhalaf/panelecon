#' Panel Fixed Effects Estimation for Time-Invariant Variables
#'
#' @md
#' @description
#' Estimates panel models with time-invariant regressors using FEVD, FEF, or
#' FEF-IV methods. Standard fixed effects estimation cannot identify
#' coefficients on time-invariant variables; these methods decompose or filter
#' the unit effects to recover these coefficients.
#'
#' @param formula A formula of the form `y ~ x1 + x2 | z1 + z2` where terms
#'   before `|` are time-varying and terms after `|` are time-invariant. Both
#'   parts are processed with [stats::model.matrix()], so transformations
#'   (`log(x)`, `I(x^2)`, `x:w`, `log(y)` on the left-hand side) and factors
#'   (with their contrasts) are handled as in [stats::lm()]. Intercept columns
#'   are removed; the model intercept is estimated in stage 2.
#' @param data A data frame containing the variables.
#' @param id Character string naming the panel (individual) identifier variable.
#' @param time Character string naming the time identifier variable.
#' @param method Estimation method: `"fevd"` (default), `"fef"`, or `"fef_iv"`.
#' @param instruments For `method = "fef_iv"`, a one-sided formula specifying
#'   instrumental variables, e.g., `~ iv1 + iv2`.
#' @param vcov_beta Covariance estimator for the stage 1 (within) coefficients
#'   `beta`: `"robust"` (default) is the Arellano-type panel-robust matrix of
#'   Pesaran and Zhou (2018, eq. 18), valid under heteroskedasticity and
#'   serial correlation within panels; `"classical"` is the homoskedastic
#'   \eqn{\hat\sigma_e^2 (X'MX)^{-1}}. The chosen matrix is also used inside
#'   the Pesaran and Zhou variance of `gamma` and of the intercept.
#' @param na.action How to handle missing values. Default is `na.omit`.
#'
#' @return An object of class `"xtfifevd"` containing:
#' \describe{
#'   \item{coefficients}{Named vector of all coefficients (beta, gamma, `_cons`).
#'     For `method = "fevd"` these are the stage 3 pooled OLS estimates, which
#'     equal the FEF estimates (see Details).}
#'   \item{vcov}{Full variance-covariance matrix of `coefficients` from the
#'     Pesaran and Zhou (2018) derivation, including the covariances between
#'     `beta`, `gamma` and the intercept (see Details).}
#'   \item{beta}{Coefficients on time-varying variables}
#'   \item{gamma}{Coefficients on time-invariant variables}
#'   \item{intercept}{Overall intercept}
#'   \item{delta}{(FEVD only) Stage 3 coefficient on the unexplained unit
#'     effect \eqn{h_i}; equal to 1 by construction.}
#'   \item{stage3}{(FEVD only) List with the full stage 3 pooled OLS
#'     coefficient vector (`coefficients`, including `h`), the naive OLS
#'     covariance matrix `vcov_naive` and standard errors `se_naive`. These
#'     naive standard errors are known to be too small for the
#'     time-invariant coefficients and are returned for reference only.}
#'   \item{fef}{(FEVD only) The stage 1 and 2 (FEF) coefficients.}
#'   \item{residuals}{Idiosyncratic (within) residuals from stage 1}
#'   \item{fitted.values}{Fitted values \eqn{x_{it}'\hat\beta + \bar z_i'\hat\gamma + \hat\alpha}
#'     (unit effects excluded)}
#'   \item{sigma2_e}{Variance of the idiosyncratic error}
#'   \item{sigma2_u}{Variance of the unexplained unit effect: the stage 2
#'     residual variance minus \eqn{\hat\sigma_e^2 \, \mathrm{mean}(1/T_i)},
#'     truncated at zero}
#'   \item{vcov_beta}{The `vcov_beta` choice used}
#'   \item{V_beta_robust, V_beta_classical}{Both covariance matrices of `beta`}
#'   \item{V_gamma_pz}{The `gamma` block of `vcov` (Pesaran and Zhou eq. 17 or 51)}
#'   \item{stage2_residuals}{Unit-level stage 2 residuals}
#'   \item{N}{Total number of observations}
#'   \item{N_g}{Number of groups (panels)}
#'   \item{T_bar}{Average time periods per panel}
#'   \item{balanced}{Logical, whether the panel is balanced}
#'   \item{method}{Estimation method used}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' ## Model
#' The panel model is:
#' \deqn{y_{it} = \alpha + \alpha_i + x_{it}'\beta + z_i'\gamma + \varepsilon_{it}}
#'
#' where \eqn{x_{it}} are time-varying regressors, \eqn{z_i} are time-invariant
#' regressors, and \eqn{\alpha_i} are individual effects that may be
#' correlated with \eqn{x_{it}}.
#'
#' ## Stage 1 (all methods)
#' Within (fixed effects) regression of \eqn{y_{it}} on \eqn{x_{it}} yields
#' \eqn{\hat{\beta}} and the time-averaged FE residuals
#' \eqn{\bar u_i = \bar y_i - \bar x_i'\hat\beta} (Pesaran and Zhou 2018, eq. 3).
#'
#' ## Stage 2
#' \itemize{
#'   \item **FEF**: unit-level OLS of \eqn{\bar u_i} on an intercept and
#'     \eqn{z_i} (Pesaran and Zhou eq. 4 and 5).
#'   \item **FEF-IV**: unit-level 2SLS using instruments \eqn{r_i}
#'     (Pesaran and Zhou eq. 48).
#'   \item **FEVD**: as FEF; the residuals \eqn{h_i} are the unexplained part
#'     of the unit effect (Plumper and Troeger 2007, eq. 6).
#' }
#'
#' ## Stage 3 (FEVD only)
#' Pooled OLS of \eqn{y_{it}} on an intercept, \eqn{x_{it}}, \eqn{z_i} and
#' \eqn{h_i} (Plumper and Troeger 2007, eq. 7). Plumper and Troeger's stage 2
#' equation (5) is printed without an intercept, but their \eqn{\hat u_i}
#' (eq. 4) contains the model constant and their stage 3 equation (7) has an
#' intercept; the package therefore includes an intercept in stage 2, which
#' makes the stage 3 estimates identical to FEF and the coefficient
#' \eqn{\delta} on \eqn{h_i} identically equal to 1 (Pesaran and Zhou 2018,
#' Proposition 3). The identity holds in balanced and unbalanced panels:
#' with \eqn{(a, \hat\beta, \hat\gamma, 1)} the stage 3 residuals are the
#' within residuals, which sum to zero within every unit and are orthogonal
#' to \eqn{x_{it}}, so they satisfy the stage 3 normal equations exactly.
#' Without the stage 2 intercept the FEVD estimator is in general biased
#' (Pesaran and Zhou 2018, Section 3.4). Plumper and Troeger are silent on
#' unbalanced panels; the package runs stage 2 at the unit level without
#' weights, one observation per panel unit (Pesaran and Zhou's FEF), rather
#' than at the observation level, where units would be weighted by
#' \eqn{T_i}.
#'
#' The naive stage 3 OLS standard errors of Plumper and Troeger are too
#' small for the time-invariant coefficients because they ignore that
#' \eqn{h_i} is a generated regressor (Breusch, Ward, Nguyen and Kompas
#' 2010, Theorem 3; Greene 2011; Pesaran and Zhou 2018). In a Monte Carlo
#' check by the package author (400 replications, N = 200, T = 8, two
#' time-varying and two time-invariant regressors, AR(1) errors with
#' coefficient 0.8 and heteroskedastic across units, x correlated with the
#' unit effects) their 95 percent coverage for \eqn{\gamma} was 35 to 40
#' percent, against 92 to 95 percent for the Pesaran and Zhou standard
#' errors. They are returned in `stage3$se_naive` for reference only and are
#' never used for inference.
#'
#' ## Rarely changing variables
#' If a variable after `|` varies within panels, a warning is issued and its
#' unit (panel) mean is used as \eqn{z_i} in all stages. This is the package's
#' choice; Plumper and Troeger (2007) do not specify how rarely changing
#' variables enter stage 2.
#'
#' ## Variance estimation
#' The `gamma` block of `vcov` is Pesaran and Zhou (2018) equation 17 (FEF and
#' FEVD) or equation 51 (FEF-IV), which account for the estimation
#' uncertainty of \eqn{\hat\beta} through the matrix `vcov_beta`. The `beta`
#' block is `vcov_beta` itself (robust eq. 18 by default). The remaining
#' blocks follow from Pesaran and Zhou eq. (A.11),
#' \deqn{\hat\gamma - \gamma = Q_{zz}^{-1}\left[N^{-1}\sum_i (z_i - \bar z) v_i - Q_{z\bar x}(\hat\beta - \beta)\right],}
#' so that \eqn{Cov(\hat\gamma, \hat\beta) = -Q_{zz}^{-1} Q_{z\bar x} Var(\hat\beta)},
#' and from eq. (5), \eqn{\hat\alpha = \bar u - \bar z'\hat\gamma}, by the delta
#' method with \eqn{c = \bar x - Q_{z\bar x}' Q_{zz}^{-1} \bar z}:
#' \eqn{Var(\hat\alpha) = N^{-2}\sum_i w_i^2 \hat v_i^2 + c' Var(\hat\beta) c},
#' \eqn{w_i = 1 - \bar z' Q_{zz}^{-1}(z_i - \bar z)}, with the corresponding
#' covariances. For FEF-IV the same derivation is used with
#' \eqn{Q_{zz}^{-1}(z_i - \bar z)} replaced by \eqn{H_{zr}(r_i - \bar r)}
#' and \eqn{Q_{z\bar x}} by \eqn{Q_{r\bar x}}. As in Pesaran and Zhou's
#' eq. (17), the cross term between the unit-level scores and
#' \eqn{\hat\beta} (the term defined in their eq. 15, negligible under
#' their condition 16) is dropped in \eqn{Var(\hat\alpha)} and in \eqn{Cov(\hat\gamma, \hat\beta)}.
#'
#' @examples
#' # Simulate panel data
#' set.seed(123)
#' N <- 100  # panels
#' T <- 10   # time periods
#' n <- N * T
#'
#' # Generate data
#' id <- rep(1:N, each = T)
#' time <- rep(1:T, N)
#' alpha_i <- rep(rnorm(N), each = T)  # Fixed effects
#' z <- rep(rnorm(N), each = T)        # Time-invariant
#' x <- rnorm(n)                        # Time-varying
#' y <- 1 + 2 * x + 0.5 * z + alpha_i + rnorm(n, sd = 0.5)
#'
#' data <- data.frame(id = id, time = time, y = y, x = x, z = z)
#'
#' # Estimate with different methods
#' fit_fevd <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time")
#' summary(fit_fevd)
#' fit_fevd$delta   # equals 1 by construction
#'
#' fit_fef <- xtfifevd(y ~ x | z, data = data, id = "id", time = "time",
#'                     method = "fef")
#' summary(fit_fef)
#'
#' # Transformations in the formula are allowed
#' fit_log <- fef(y ~ x + I(x^2) | z, data = data, id = "id", time = "time")
#' coef(fit_log)
#'
#' @references
#' Breusch, T., Ward, M. B., Nguyen, H. and Kompas, T. (2010). On the
#' fixed-effects vector decomposition. MPRA Paper No. 21452.
#' \url{https://mpra.ub.uni-muenchen.de/21452/}
#'
#' Greene, W. H. (2011). Fixed Effects Vector Decomposition: A Magical Solution
#' to the Problem of Time-Invariant Variables in Fixed Effects Models?
#' \emph{Political Analysis}, 19(2), 135-146.
#' \doi{10.1093/pan/mpq034}
#'
#' Plumper, T. and Troeger, V. E. (2007). Efficient Estimation of Time-Invariant
#' and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
#' Effects. \emph{Political Analysis}, 15(2), 124-139.
#' \doi{10.1093/pan/mpm002}
#'
#' Pesaran, M. H. and Zhou, Q. (2018). Estimation of time-invariant effects in
#' static panel data models. \emph{Econometric Reviews}, 37(10), 1137-1171.
#' \doi{10.1080/07474938.2016.1222225}
#'
#' @seealso [fevd()], [fef()], [fef_iv()], [bw_ratio()]
#'
#' @export
xtfifevd <- function(formula, data, id, time,
                     method = c("fevd", "fef", "fef_iv"),
                     instruments = NULL,
                     vcov_beta = c("robust", "classical"),
                     na.action = na.omit) {

  call <- match.call()
  method <- match.arg(method)
  vcov_beta <- match.arg(vcov_beta)

  # Parse the formula: y ~ x1 + x2 | z1 + z2
  parsed <- .parse_formula(formula, data, id, time, instruments, na.action)

  # Dispatch to appropriate estimator
  result <- switch(method,
                   "fevd" = .estimate_fevd(parsed, vcov_beta),
                   "fef" = .estimate_fef(parsed, vcov_beta),
                   "fef_iv" = .estimate_fef_iv(parsed, vcov_beta))

  result$call <- call
  result$formula <- formula
  result$method <- toupper(method)
  result$method <- gsub("_", "-", result$method)

  class(result) <- "xtfifevd"
  result
}


#' @describeIn xtfifevd FEVD estimation (3-stage, Plumper and Troeger 2007,
#'   with Pesaran and Zhou 2018 standard errors)
#' @export
fevd <- function(formula, data, id, time, vcov_beta = c("robust", "classical"),
                 na.action = na.omit) {
  xtfifevd(formula, data, id, time, method = "fevd",
           vcov_beta = match.arg(vcov_beta), na.action = na.action)
}


#' @describeIn xtfifevd FEF estimation (2-stage, Pesaran and Zhou 2018)
#' @export
fef <- function(formula, data, id, time, vcov_beta = c("robust", "classical"),
                na.action = na.omit) {
  xtfifevd(formula, data, id, time, method = "fef",
           vcov_beta = match.arg(vcov_beta), na.action = na.action)
}


#' @describeIn xtfifevd FEF-IV estimation with instruments
#' @export
fef_iv <- function(formula, data, id, time, instruments,
                   vcov_beta = c("robust", "classical"),
                   na.action = na.omit) {
  if (missing(instruments) || is.null(instruments)) {
    stop("FEF-IV requires instruments. Provide 'instruments' argument.")
  }
  xtfifevd(formula, data, id, time, method = "fef_iv",
           instruments = instruments, vcov_beta = match.arg(vcov_beta),
           na.action = na.action)
}
