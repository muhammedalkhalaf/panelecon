\name{xtqsh}
\alias{xtqsh}
\title{Quantile Regression Slope Homogeneity Test for Panel Data}
\description{
  Tests the null hypothesis of slope homogeneity in panel quantile regressions.
  Implements the S-hat and D-hat statistics of Galvao et al. (2017).
}
\usage{
xtqsh(formula, data, index, tau, bw = "hallsheather", marginal = FALSE)
}
\arguments{
  \item{formula}{A formula of the form \code{y ~ x1 + x2 + ...}.}
  \item{data}{A data frame containing the panel data in long format.}
  \item{index}{Character vector of length 2: \code{c("id_var", "time_var")}.}
  \item{tau}{Numeric vector of quantile levels, each strictly between 0 and 1.}
  \item{bw}{Bandwidth method: \code{"hallsheather"} (default) or \code{"bofinger"}.}
  \item{marginal}{Logical. If \code{TRUE}, compute per-variable marginal tests. Default \code{FALSE}.}
}
\value{
  An object of class \code{"xtqsh"} containing test statistics and p-values.
}
\references{
Galvao, A.F., Juhl, T., Montes-Rojas, G. and Olmo, J. (2017).
Testing Slope Homogeneity in Quantile Regression Panel Data.
\emph{Journal of Financial Econometrics}, 16(2), 211-243.
}
\examples{
\donttest{
set.seed(42)
n <- 10; tt <- 20
dat <- data.frame(
  id   = rep(1:n, each = tt),
  time = rep(1:tt, times = n),
  y    = rnorm(n * tt),
  x1   = rnorm(n * tt)
)
res <- xtqsh(y ~ x1, data = dat, index = c("id", "time"), tau = 0.5)
print(res)
}
}
\keyword{htest}
