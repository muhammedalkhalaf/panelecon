\name{xtmispanel}
\alias{xtmispanel}
\title{Missing Data Detection and Imputation for Panel Data}
\usage{
xtmispanel(
  data,
  vars = NULL,
  index,
  detect = TRUE,
  test = FALSE,
  impute = NULL,
  target = NULL,
  new_var = NULL,
  sensitivity = FALSE,
  knn_k = 5L
)
}
\arguments{
\item{data}{A \code{data.frame} in long format.}

\item{vars}{Character vector of variable names to analyse.  If \code{NULL}
(default), all numeric columns except the index are used.}

\item{index}{Character vector of length 2: \code{c("panel_id", "time_id")}.}

\item{detect}{Logical. Run Module 1 (detection tables, default \code{TRUE}).}

\item{test}{Logical. Run Module 2 (MCAR/MAR mechanism tests,
default \code{FALSE}).}

\item{impute}{Character or \code{NULL}. If a method name is given, run
Module 3 (imputation).  Supported methods:
\code{"mean"}, \code{"median"}, \code{"locf"}, \code{"nocb"},
\code{"linear"}, \code{"spline"}, \code{"pmm"}, \code{"hotdeck"},
\code{"knn"}, \code{"rf"}, \code{"em"}.}

\item{target}{Character. Name of the variable to impute (required when
\code{impute} is not \code{NULL}).}

\item{new_var}{Character. Name of the output imputed variable
(default \code{"\{target\}_imp"}).}

\item{sensitivity}{Logical. Run Module 4 (sensitivity analysis across
all imputation methods, default \code{FALSE}).}

\item{knn_k}{Integer. Number of neighbours for KNN imputation (default 5).}
}
\value{
A list (invisibly) with components:
\describe{
  \item{\code{detect}}{Summary statistics per variable/panel/period.}
  \item{\code{test}}{MCAR and MAR test results.}
  \item{\code{imputed}}{The \code{data} frame augmented with the imputed
    column (when imputation is requested).}
  \item{\code{impute_stats}}{Summary comparing original vs imputed.}
  \item{\code{sensitivity}}{Sensitivity analysis results.}
}
}
\description{
Detects, diagnoses, and imputes missing values in panel (longitudinal)
data sets.  The function can produce summary tables (Module 1), test the
missingness mechanism (Module 2), impute a target variable (Module 3),
and run a cross-method sensitivity analysis (Module 4).
}
\examples{
set.seed(1)
df <- data.frame(
  id   = rep(1:4, each = 8),
  time = rep(1:8, times = 4),
  y    = c(rnorm(32))
)
# introduce some NAs
df$y[c(3, 11, 20)] <- NA
res <- xtmispanel(df, vars = "y", index = c("id", "time"), detect = TRUE)
}
\references{
Little, R. J. A. (1988). A test of missing completely at random for
multivariate data with missing values.
\emph{Journal of the American Statistical Association}, 83(404), 1198-1202.
\doi{10.1080/01621459.1988.10478714}
}
