#' Missing Data Detection and Imputation for Panel Data
#'
#' Detects, diagnoses, and imputes missing values in panel (longitudinal)
#' data sets.  The function can produce summary tables (Module 1), test the
#' missingness mechanism (Module 2), impute a target variable (Module 3),
#' and run a cross-method sensitivity analysis (Module 4).
#'
#' @param data A \code{data.frame} in long format.
#' @param vars Character vector of variable names to analyse.  If \code{NULL}
#'   (default), all numeric columns except the index are used.
#' @param index Character vector of length 2: \code{c("panel_id", "time_id")}.
#' @param detect Logical. Run Module 1 (detection tables, default \code{TRUE}).
#' @param test Logical. Run Module 2 (MCAR/MAR mechanism tests,
#'   default \code{FALSE}).
#' @param impute Character or \code{NULL}. If a method name is given, run
#'   Module 3 (imputation).  Supported methods:
#'   \code{"mean"}, \code{"median"}, \code{"locf"}, \code{"nocb"},
#'   \code{"linear"}, \code{"spline"}, \code{"pmm"}, \code{"hotdeck"},
#'   \code{"knn"}, \code{"rf"}, \code{"em"}.
#' @param target Character. Name of the variable to impute (required when
#'   \code{impute} is not \code{NULL}).
#' @param new_var Character. Name of the output imputed variable
#'   (default \code{"\{target\}_imp"}).
#' @param sensitivity Logical. Run Module 4 (sensitivity analysis across
#'   all imputation methods, default \code{FALSE}).
#' @param knn_k Integer. Number of neighbours for KNN imputation (default 5).
#'
#' @return A list (invisibly) with components:
#'   \describe{
#'     \item{\code{detect}}{Summary statistics per variable/panel/period.}
#'     \item{\code{test}}{MCAR and MAR test results.}
#'     \item{\code{imputed}}{The \code{data} frame augmented with the imputed
#'       column (when imputation is requested).}
#'     \item{\code{impute_stats}}{Summary comparing original vs imputed.}
#'     \item{\code{sensitivity}}{Sensitivity analysis results.}
#'   }
#'
#' @references
#' Little, R. J. A. (1988). A test of missing completely at random for
#' multivariate data with missing values.
#' \emph{Journal of the American Statistical Association}, 83(404), 1198-1202.
#' \doi{10.1080/01621459.1988.10478714}
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   id   = rep(1:4, each = 8),
#'   time = rep(1:8, times = 4),
#'   y    = c(rnorm(32))
#' )
#' # introduce some NAs
#' df$y[c(3, 11, 20)] <- NA
#' res <- xtmispanel(df, vars = "y", index = c("id", "time"), detect = TRUE)
#'
#' @export
xtmispanel <- function(data,
                       vars        = NULL,
                       index,
                       detect      = TRUE,
                       test        = FALSE,
                       impute      = NULL,
                       target      = NULL,
                       new_var     = NULL,
                       sensitivity = FALSE,
                       knn_k       = 5L) {

  ## ---- Validate --------------------------------------------------------
  if (!inherits(data, "data.frame")) stop("'data' must be a data.frame.")
  if (length(index) != 2L) stop("'index' must be a character vector of length 2.")
  pid_col <- index[1L]; tid_col <- index[2L]
  if (!all(c(pid_col, tid_col) %in% names(data)))
    stop("index columns not found in 'data'.")

  ## Determine analysis variables
  if (is.null(vars)) {
    num_cols <- vapply(data, is.numeric, logical(1L))
    vars <- setdiff(names(data)[num_cols], c(pid_col, tid_col))
  }
  vars <- intersect(vars, names(data))
  if (length(vars) == 0L) stop("No valid numeric variables found.")

  pid <- data[[pid_col]]
  tid <- data[[tid_col]]

  .xmp_header(pid_col, tid_col, pid, tid, data)
  out <- list()

  ## ---- Module 1: Detection --------------------------------------------
  if (detect) {
    out$detect <- .xmp_detect(data, vars, pid, tid, pid_col, tid_col)
  }

  ## ---- Module 2: Mechanism tests --------------------------------------
  if (test) {
    out$test <- .xmp_tests(data, vars, pid, tid)
  }

  ## ---- Module 3: Imputation -------------------------------------------
  if (!is.null(impute)) {
    if (is.null(target) || !target %in% names(data))
      stop("'target' must be the name of a variable in 'data'.")
    if (is.null(new_var)) new_var <- paste0(target, "_imp")
    res_imp <- .xmp_impute(data, target, pid, tid, method = impute, knn_k = knn_k)
    data[[new_var]] <- res_imp
    out$imputed      <- data
    out$impute_stats <- .xmp_compare(data[[target]], data[[new_var]],
                                     target, new_var, impute)
  }

  ## ---- Module 4: Sensitivity ------------------------------------------
  if (sensitivity) {
    if (is.null(target) || !target %in% names(data))
      stop("'target' must be set for sensitivity analysis.")
    out$sensitivity <- .xmp_sensitivity(data, target, pid, tid, knn_k = knn_k)
  }

  invisible(out)
}


## ==========================================================================
## Internal helpers
## ==========================================================================

#' @keywords internal
.xmp_header <- function(pid_col, tid_col, pid, tid, data) {
  N  <- length(unique(pid))
  T2 <- length(unique(tid))
  message(strrep("-", 72))
  message("  xtmispanel -- Missing Data Diagnostics for Panel Data  v1.0.0")
  message(strrep("-", 72))
  message("  Panel id: ", pid_col, "  |  Time id: ", tid_col)
  message("  Panels: ", N, "  |  Time periods: ", T2,
          "  |  Rows: ", nrow(data))
  message(strrep("-", 72))
}

#' @keywords internal
.xmp_detect <- function(data, vars, pid, tid, pid_col, tid_col) {
  line <- strrep("-", 72)
  message(line)
  message("  TABLE 1: Missing Data Summary by Variable")
  message(line)
  message(sprintf("  %-16s %8s %8s %7s %10s %10s  %s",
                  "Variable", "N_Total", "N_Miss", "%Miss", "Mean", "SD", "Status"))
  message(line)

  out_vars <- list()
  total_miss <- 0L; total_obs <- 0L

  for (v in vars) {
    x     <- data[[v]]
    ntot  <- length(x)
    nmiss <- sum(is.na(x))
    pct   <- if (ntot > 0) nmiss / ntot * 100 else 0
    vmean <- mean(x, na.rm = TRUE)
    vsd   <- stats::sd(x, na.rm = TRUE)
    if (is.na(vsd)) vsd <- 0
    status <- if (pct == 0) "Complete"
    else if (pct <= 5) "Low"
    else if (pct <= 20) "Moderate"
    else if (pct <= 50) "High"
    else "Severe"
    message(sprintf("  %-16s %8d %8d %6.1f%% %10.3f %10.3f  %s",
                    substr(v, 1L, 16L), ntot, nmiss, pct, vmean, vsd, status))
    out_vars[[v]] <- list(n = ntot, n_miss = nmiss, pct = pct,
                          mean = vmean, sd = vsd, status = status)
    total_miss <- total_miss + nmiss
    total_obs  <- total_obs  + ntot
  }
  message(line)
  message(sprintf("  %-16s %8d %8d %6.1f%%",
                  "Overall", total_obs, total_miss,
                  if (total_obs > 0) total_miss / total_obs * 100 else 0))
  message(line)

  ## Per-panel summary
  panels  <- sort(unique(pid))
  message(line)
  message("  TABLE 2: Missing Data by Panel")
  message(line)
  message(sprintf("  %-14s %8s %8s %7s  %s",
                  "Panel", "N_Obs", "N_Miss", "%Miss", "Status"))
  message(line)
  out_panels <- list()
  for (p in panels) {
    idx   <- pid == p
    pmiss <- sum(is.na(unlist(data[idx, vars, drop = FALSE])))
    pobs  <- sum(idx) * length(vars)
    ppct  <- if (pobs > 0) pmiss / pobs * 100 else 0
    pstat <- if (ppct == 0) "Complete"
    else if (ppct <= 10) "Low"
    else if (ppct <= 30) "Moderate"
    else "High"
    message(sprintf("  %-14s %8d %8d %6.1f%%  %s",
                    substr(as.character(p), 1L, 14L), pobs, pmiss, ppct, pstat))
    out_panels[[as.character(p)]] <- list(pobs = pobs, pmiss = pmiss, pct = ppct)
  }
  message(line)

  list(variables = out_vars, panels = out_panels,
       total_missing = total_miss, total_obs = total_obs)
}

#' @keywords internal
.xmp_tests <- function(data, vars, pid, tid) {
  line <- strrep("-", 72)
  message(line)
  message("  MODULE 2: Missing Data Mechanism Tests")
  message(line)

  ## Approximate Little MCAR test
  message("  Test 1: Approximate Little MCAR Test")
  chi2_total <- 0; df_total <- 0

  for (v in vars) {
    if (sum(is.na(data[[v]])) == 0L) next
    gmean <- mean(data[[v]], na.rm = TRUE)
    gvar  <- stats::var(data[[v]], na.rm = TRUE)
    if (is.na(gvar) || gvar == 0) next
    for (w in vars) {
      if (v == w) next
      if (sum(is.na(data[[w]])) == 0L) next
      m_miss <- mean(data[[v]][is.na(data[[w]])], na.rm = TRUE)
      m_obs  <- mean(data[[v]][!is.na(data[[w]])], na.rm = TRUE)
      n1 <- sum(is.na(data[[w]]) & !is.na(data[[v]]))
      n2 <- sum(!is.na(data[[w]]) & !is.na(data[[v]]))
      if (n1 < 2L || n2 < 2L) next
      se <- sqrt(gvar * (1 / n1 + 1 / n2))
      if (se > 0) {
        z <- (m_miss - m_obs) / se
        chi2_total <- chi2_total + z^2
        df_total   <- df_total + 1L
      }
    }
  }

  mcar_result <- NULL
  if (df_total > 0L) {
    pval <- stats::pchisq(chi2_total, df = df_total, lower.tail = FALSE)
    message(sprintf("  chi2(%d) = %.3f   p = %.4f   %s",
                    df_total, chi2_total, pval,
                    if (pval < 0.05) "REJECT H0: NOT MCAR" else "Fail to reject: possibly MCAR"))
    mcar_result <- list(chi2 = chi2_total, df = df_total, pval = pval)
  } else {
    message("  Test could not be computed (insufficient variation).")
  }

  ## MAR logistic test
  message(line)
  message("  Test 2: MAR Logistic Regression Test")
  message(sprintf("  %-16s %10s %10s %10s  %s",
                  "Variable", "chi2", "p-value", "Pseudo-R2", "Conclusion"))
  message(line)

  mar_results <- list()
  for (v in vars) {
    nmv <- sum(is.na(data[[v]]))
    if (nmv == 0L) { message(sprintf("  %-16s %10s", substr(v,1,16), "No missing")); next }
    mis_ind <- as.integer(is.na(data[[v]]))
    preds   <- setdiff(vars, v)
    if (length(preds) == 0L) next
    pred_df <- data[, preds, drop = FALSE]
    complete_rows <- stats::complete.cases(pred_df)
    if (sum(complete_rows) < length(preds) + 2L) next
    fit <- tryCatch(
      stats::glm(mis_ind[complete_rows] ~ as.matrix(pred_df[complete_rows, ]),
                 family = stats::binomial()),
      error = function(e) NULL
    )
    if (is.null(fit)) { message(sprintf("  %-16s  FAILED", substr(v,1,16))); next }
    null_dev <- fit$null.deviance; res_dev <- fit$deviance
    chi2  <- null_dev - res_dev
    df_m  <- fit$df.null - fit$df.residual
    pval  <- stats::pchisq(chi2, df = df_m, lower.tail = FALSE)
    pr2   <- 1 - res_dev / null_dev
    conc  <- if (!is.na(pval) && pval < 0.05) "MAR" else "MCAR"
    message(sprintf("  %-16s %10.3f %10.4f %10.4f  %s",
                    substr(v, 1L, 16L), chi2, pval, pr2, conc))
    mar_results[[v]] <- list(chi2 = chi2, df = df_m, pval = pval, conclusion = conc)
  }
  message(line)
  list(mcar = mcar_result, mar = mar_results)
}

#' @keywords internal
.xmp_impute <- function(data, target, pid, tid, method, knn_k) {
  x       <- data[[target]]
  out     <- x
  panels  <- sort(unique(pid))

  method <- tolower(trimws(method))

  for (p in panels) {
    idx  <- which(pid == p)
    vals <- x[idx]
    ord  <- order(tid[idx])
    vals <- vals[ord]

    vals_imp <- switch(method,
      mean    = { m <- mean(vals, na.rm = TRUE); ifelse(is.na(vals), m, vals) },
      median  = { m <- stats::median(vals, na.rm = TRUE); ifelse(is.na(vals), m, vals) },
      locf    = .xmp_locf(vals),
      nocb    = rev(.xmp_locf(rev(vals))),
      linear  = .xmp_linear_interp(vals),
      spline  = .xmp_spline_interp(vals),
      pmm     = .xmp_pmm(vals),
      hotdeck = .xmp_hotdeck(vals),
      knn     = .xmp_knn(data, target, pid, tid, p, knn_k),
      rf      = .xmp_rf(data, target, pid, tid, p),
      em      = .xmp_em(vals),
      {
        message("  xtmispanel: unknown method '", method, "'; using mean.")
        m <- mean(vals, na.rm = TRUE); ifelse(is.na(vals), m, vals)
      }
    )
    out[idx[ord]] <- vals_imp
  }
  out
}

## --- Imputation method implementations ---

#' @keywords internal
.xmp_locf <- function(x) {
  for (i in seq_along(x)) {
    if (is.na(x[i]) && i > 1L) x[i] <- x[i - 1L]
  }
  x
}

#' @keywords internal
.xmp_linear_interp <- function(x) {
  n   <- length(x)
  idx <- which(!is.na(x))
  if (length(idx) < 2L) return(x)
  stats::approx(idx, x[idx], xout = seq_len(n), rule = 2L)$y
}

#' @keywords internal
.xmp_spline_interp <- function(x) {
  n   <- length(x)
  idx <- which(!is.na(x))
  if (length(idx) < 3L) return(.xmp_linear_interp(x))
  tryCatch(
    stats::spline(idx, x[idx], xout = seq_len(n), method = "natural")$y,
    error = function(e) .xmp_linear_interp(x)
  )
}

#' @keywords internal
.xmp_pmm <- function(x) {
  obs_vals <- x[!is.na(x)]
  if (length(obs_vals) == 0L) return(x)
  out <- x
  for (i in which(is.na(x))) {
    i_ctx <- mean(obs_vals)  # simplified: use global mean as donor context
    dists <- abs(obs_vals - i_ctx)
    out[i] <- obs_vals[which.min(dists)]
  }
  out
}

#' @keywords internal
.xmp_hotdeck <- function(x) {
  obs_vals <- x[!is.na(x)]
  if (length(obs_vals) == 0L) return(x)
  out <- x
  miss_idx <- which(is.na(x))
  out[miss_idx] <- sample(obs_vals, length(miss_idx), replace = TRUE)
  out
}

#' @keywords internal
.xmp_knn <- function(data, target, pid, tid, panel, k) {
  idx      <- which(pid == panel)
  vals     <- data[[target]][idx]
  ord      <- order(tid[idx])
  vals     <- vals[ord]
  miss_pos <- which(is.na(vals))
  if (length(miss_pos) == 0L) return(vals)

  # Use other panels' values at same time positions as donors
  other_idx <- which(pid != panel)
  donors    <- data[[target]][other_idx]
  donors    <- donors[!is.na(donors)]
  if (length(donors) == 0L) donors <- stats::na.omit(vals)
  if (length(donors) == 0L) return(vals)

  obs_vals <- vals[!is.na(vals)]
  for (i in miss_pos) {
    # time-position context: average of window
    window <- vals[max(1L, i - 2L):min(length(vals), i + 2L)]
    ctx    <- mean(window, na.rm = TRUE)
    if (is.na(ctx)) ctx <- mean(donors)
    dists  <- abs(donors - ctx)
    nn     <- utils::head(order(dists), k)
    vals[i] <- mean(donors[nn])
  }
  vals
}

#' @keywords internal
.xmp_rf <- function(data, target, pid, tid, panel) {
  ## Simplified random-forest-like imputation via bagged linear models
  idx      <- which(pid == panel)
  vals     <- data[[target]][idx]
  ord      <- order(tid[idx])
  vals     <- vals[ord]
  miss_pos <- which(is.na(vals))
  if (length(miss_pos) == 0L) return(vals)

  obs_idx  <- which(!is.na(vals))
  if (length(obs_idx) < 3L) return(.xmp_linear_interp(vals))

  t_idx <- seq_along(vals)
  preds <- numeric(length(miss_pos))
  B <- 20L
  for (b in seq_len(B)) {
    samp <- sample(obs_idx, replace = TRUE)
    fit  <- tryCatch(stats::lm(vals[samp] ~ t_idx[samp]),
                     error = function(e) NULL)
    if (!is.null(fit)) {
      preds <- preds + stats::predict(fit, newdata = data.frame(`t_idx[samp]` = t_idx[miss_pos]))
    }
  }
  vals[miss_pos] <- preds / B
  vals
}

#' @keywords internal
.xmp_em <- function(x) {
  ## Single-variable EM: iterate between mean estimate and imputation
  obs  <- x[!is.na(x)]
  if (length(obs) == 0L) return(x)
  mu   <- mean(obs); sig2 <- stats::var(obs)
  if (is.na(sig2) || sig2 == 0) sig2 <- 1
  out  <- x
  for (iter in seq_len(50L)) {
    out[is.na(x)] <- mu
    mu_new  <- mean(out)
    sig2_new <- stats::var(out)
    if (abs(mu_new - mu) < 1e-8) break
    mu <- mu_new; sig2 <- sig2_new
  }
  out
}

#' @keywords internal
.xmp_compare <- function(orig, imp, orig_name, imp_name, method) {
  line <- strrep("-", 72)
  message(line)
  message("  Imputation Results: method = ", method)
  message(line)
  nmiss  <- sum(is.na(orig))
  nfill  <- sum(!is.na(imp)) - sum(!is.na(orig))
  nstill <- sum(is.na(imp))
  message("  Missing before   : ", nmiss)
  message("  Values filled    : ", nfill)
  message("  Still missing    : ", nstill)

  o_mean <- mean(orig, na.rm = TRUE); i_mean <- mean(imp, na.rm = TRUE)
  o_sd   <- stats::sd(orig, na.rm = TRUE); i_sd <- stats::sd(imp, na.rm = TRUE)
  d_mean <- if (!is.na(o_mean) && o_mean != 0) (i_mean - o_mean) / abs(o_mean) * 100 else NA_real_

  message(sprintf("  %-12s %10s %10s %10s", "Stat", "Original", "Imputed", "Delta%"))
  message(sprintf("  %-12s %10.4f %10.4f %9.2f%%", "Mean", o_mean, i_mean, d_mean))
  message(sprintf("  %-12s %10.4f %10.4f", "SD", o_sd, i_sd))

  corr <- tryCatch(stats::cor(orig, imp, use = "complete.obs"), error = function(e) NA_real_)
  message(sprintf("  Correlation (obs pairs): %.4f", corr))
  message(line)

  list(n_miss = nmiss, n_filled = nfill, n_remain = nstill,
       orig_mean = o_mean, imp_mean = i_mean, d_mean_pct = d_mean,
       orig_sd = o_sd, imp_sd = i_sd, correlation = corr)
}

#' @keywords internal
.xmp_sensitivity <- function(data, target, pid, tid, knn_k) {
  methods <- c("mean", "median", "locf", "nocb", "linear", "spline",
               "pmm", "hotdeck", "knn", "rf", "em")
  orig      <- data[[target]]
  orig_mean <- mean(orig, na.rm = TRUE)

  line <- strrep("-", 72)
  message(line)
  message("  MODULE 4: Sensitivity Analysis -- Comparing Imputation Methods")
  message(line)
  message(sprintf("  %-12s %10s %10s %10s %10s %10s",
                  "Method", "Mean", "SD", "Min", "Max", "Delta Mean%"))
  message(line)

  results <- list()
  best_method <- methods[1L]; best_d <- Inf

  for (m in methods) {
    imp <- tryCatch(
      .xmp_impute(data, target, pid, tid, method = m, knn_k = knn_k),
      error = function(e) NULL
    )
    if (is.null(imp)) {
      message(sprintf("  %-12s  FAILED", m)); next
    }
    imean <- mean(imp, na.rm = TRUE); isd <- stats::sd(imp, na.rm = TRUE)
    imin  <- min(imp, na.rm = TRUE);  imax <- max(imp, na.rm = TRUE)
    d     <- if (!is.na(orig_mean) && orig_mean != 0)
      (imean - orig_mean) / abs(orig_mean) * 100 else 0
    message(sprintf("  %-12s %10.4f %10.4f %10.4f %10.4f %9.2f%%",
                    m, imean, isd, imin, imax, d))
    results[[m]] <- list(mean = imean, sd = isd, d_mean_pct = d)
    if (abs(d) < best_d) { best_d <- abs(d); best_method <- m }
  }
  message(line)
  message("  Best method (lowest distribution distortion): ", best_method)
  message("  Mean change: ", round(best_d, 2), "%")
  message(line)

  list(results = results, best_method = best_method, best_d_pct = best_d)
}

