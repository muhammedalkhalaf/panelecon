\name{xtpretest}
\alias{xtpretest}
\title{Comprehensive Panel Data Pre-Testing Suite}
\usage{
xtpretest(
  data,
  formula,
  index,
  tests = "ALL",
  level = 0.05
)
}
\arguments{
\item{data}{A \code{data.frame} in long format containing all variables.}

\item{formula}{A two-sided formula of the form \code{y ~ x1 + x2 + ...}
specifying the dependent and independent variables.}

\item{index}{Character vector of length 2: \code{c("panel_id", "time_id")}.}

\item{tests}{Character vector. Which modules to run. Possible values:
\code{"summary"}, \code{"hsiao"}, \code{"robust"}, \code{"heterogeneity"},
\code{"csd"}.  Default \code{"ALL"} runs everything.}

\item{level}{Numeric. Significance level for decisions (default 0.05).}
}
\value{
A list (invisibly) with components:
\describe{
  \item{\code{summary}}{Panel summary statistics.}
  \item{\code{hsiao}}{Hsiao homogeneity F-test results.}
  \item{\code{robust}}{Robust HC1 F-test results.}
  \item{\code{swamy}}{Swamy heterogeneity test results.}
  \item{\code{csd}}{Cross-sectional dependence test results.}
  \item{\code{recommendation}}{Character. Suggested estimator.}
}
}
\description{
Performs a full battery of panel data pre-tests: Hsiao (2014) homogeneity
F-tests, robust (HC1) versions, Swamy (1970) parameter heterogeneity test,
cross-sectional dependence (Pesaran 2004), and panel summary statistics.
}
\examples{
set.seed(10)
n <- 5; t <- 10
df <- data.frame(
  id   = rep(1:n, each = t),
  time = rep(1:t, times = n),
  y    = rnorm(n * t),
  x1   = rnorm(n * t)
)
res <- xtpretest(df, y ~ x1, index = c("id", "time"),
                 tests = c("hsiao", "csd"))
}
\references{
Hsiao, C. (2014). \emph{Analysis of Panel Data} (3rd ed.).
Cambridge University Press.
\doi{10.1017/CBO9781139839327}

Swamy, P. A. V. B. (1970). Efficient inference in a random coefficient
regression model. \emph{Econometrica}, 38(2), 311-323.
\doi{10.2307/1909405}

Pesaran, M. H. (2004). General diagnostic tests for cross section
dependence in panels. Cambridge Working Paper in Economics, No. 0435.
\doi{10.2139/ssrn.572504}
}
