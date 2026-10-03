#' Between/Within SD Ratio for Time-Invariant Variables
#'
#' @md
#' @importFrom stats sd
#' @description
#' Computes the between-panel and within-panel standard deviations for
#' specified variables, along with their ratio. This diagnostic helps assess
#' whether FEVD/FEF methods may improve upon standard FE estimation.
#'
#' @param data A data frame containing the panel data.
#' @param variables A character vector of variable names to analyze.
#' @param id Character string naming the panel identifier variable.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{variable}{Variable name}
#'   \item{sd_between}{Between-panel standard deviation}
#'   \item{sd_within}{Within-panel standard deviation}
#'   \item{bw_ratio}{Ratio of between to within SD}
#' }
#'
#' @details
#' For truly time-invariant variables, the within-panel SD should be zero
#' (or near-zero due to numerical precision), giving an infinite ratio.
#' Plumper and Troeger (2007, p. 134) define the b/w ratio as the between SD
#' of a variable divided by its within SD, which is what this function
#' reports.
#'
#' Plumper and Troeger (2007, Section 6) compare the root mean squared error
#' of the FE and FEVD estimators of a rarely changing variable in a Monte
#' Carlo design with N = 30 units and T = 20 periods (their footnote 13
#' reports that the findings for rarely changing variables were confirmed
#' for N in 15 to 100 and T in 20 to 100). The b/w ratio above
#' which FEVD has the lower RMSE depends on the (unobservable) correlation
#' between the rarely changing variable and the unit effects: about 0.2 when
#' the correlation is 0, about 1.7 at a correlation of 0.3, about 2.8 at 0.5
#' and close to 3.8 at 0.8 (their Fig. 4, pp. 135-136). The often quoted
#' threshold of 1.7 is therefore a result for corr(z, u) = 0.3 from their
#' Fig. 4 (N = 30, T = 20), and Plumper and Troeger themselves write that
#' they "cannot offer a simple rule of thumb", adding (p. 136) that "the
#' odds are that at a b/w ratio of at least 2.8, the variable is better
#' included into the stage 2 estimation". This function prints these
#' thresholds for orientation only.
#'
#' Two caveats from Plumper and Troeger apply. First, the FEVD (and FEF)
#' coefficient on a time-invariant or rarely changing variable is biased
#' whenever that variable is correlated with the unobserved unit effects
#' (p. 129); the bias is the usual omitted variable bias and does not vanish
#' with a high b/w ratio. Second, this correlation cannot be observed or
#' tested, because the unit effects are unobservable (p. 135); the trade-off
#' between the bias of FEVD and the inefficiency of FE for rarely changing
#' variables therefore rests on an untestable assumption.
#'
#' @examples
#' # Create example data
#' set.seed(42)
#' N <- 50
#' T <- 5
#' id <- rep(1:N, each = T)
#' z_invariant <- rep(rnorm(N), each = T)  # Truly time-invariant
#' z_slow <- rep(rnorm(N), each = T) + rnorm(N * T, sd = 0.1)  # Slowly varying
#' x_varying <- rnorm(N * T)  # Time-varying
#'
#' data <- data.frame(id = id, z_inv = z_invariant,
#'                    z_slow = z_slow, x = x_varying)
#'
#' bw_ratio(data, c("z_inv", "z_slow", "x"), id = "id")
#'
#' @references
#' Plumper, T. and Troeger, V. E. (2007). Efficient Estimation of Time-Invariant
#' and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
#' Effects. \emph{Political Analysis}, 15(2), 124-139.
#' \doi{10.1093/pan/mpm002}
#'
#' @export
bw_ratio <- function(data, variables, id) {
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  missing_vars <- setdiff(c(variables, id), names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  panel_id <- data[[id]]
  
  results <- data.frame(
    variable = character(),
    sd_between = numeric(),
    sd_within = numeric(),
    bw_ratio = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (v in variables) {
    x <- data[[v]]
    
    # Between SD: SD of panel means
    panel_means <- tapply(x, panel_id, mean, na.rm = TRUE)
    sd_between <- sd(panel_means, na.rm = TRUE)
    
    # Within SD: SD of deviations from panel means
    mean_expanded <- panel_means[as.character(panel_id)]
    deviations <- x - mean_expanded
    sd_within <- sd(deviations, na.rm = TRUE)
    
    # Ratio
    if (sd_within > 1e-10) {
      ratio <- sd_between / sd_within
    } else {
      ratio <- Inf
    }
    
    results <- rbind(results, data.frame(
      variable = v,
      sd_between = sd_between,
      sd_within = sd_within,
      bw_ratio = ratio,
      stringsAsFactors = FALSE
    ))
  }
  
  # Print nicely via message() so output can be suppressed
  lines <- c(
    "",
    "Between/Within SD Ratios",
    paste(rep("-", 60), collapse = ""),
    sprintf("%-15s %12s %12s %12s",
            "Variable", "Between SD", "Within SD", "B/W Ratio"),
    paste(rep("-", 60), collapse = "")
  )
  for (i in seq_len(nrow(results))) {
    if (is.finite(results$bw_ratio[i])) {
      lines <- c(lines, sprintf("%-15s %12.4f %12.4f %12.2f",
                  results$variable[i],
                  results$sd_between[i],
                  results$sd_within[i],
                  results$bw_ratio[i]))
    } else {
      lines <- c(lines, sprintf("%-15s %12.4f %12.4f %12s",
                  results$variable[i],
                  results$sd_between[i],
                  results$sd_within[i],
                  "Inf"))
    }
  }
  lines <- c(
    lines,
    paste(rep("-", 60), collapse = ""),
    "Note: Plumper and Troeger (2007, Fig. 4; N = 30, T = 20) find that FEVD",
    "      has lower RMSE than FE for a rarely changing variable when its",
    "      B/W SD ratio exceeds about 0.2 (corr(z, u) = 0), 1.7 (0.3),",
    "      2.8 (0.5) or 3.8 (0.8). corr(z, u) is not observable or testable",
    "      and the FEVD/FEF coefficient is biased whenever it is non-zero;",
    "      PT offer no simple rule of thumb.",
    ""
  )
  message(paste(lines, collapse = "\n"))
  
  invisible(results)
}
