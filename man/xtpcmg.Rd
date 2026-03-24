\name{xtpcmg}
\alias{xtpcmg}
\title{Panel Cointegrating Polynomial Regressions via FM-OLS}
\description{
  Estimates a polynomial cointegrating regression in a panel using either
  Group-Mean FM-OLS or Pooled FM-OLS.
}
\usage{
xtpcmg(data, y, x, panel_id, time_id,
        model = c("mg", "pmg"), q = 2L,
        controls = NULL, trend = 1L, kernel = "ba",
        bw = "And91", effects = "oneway", corr_rob = FALSE)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{y}{Character. Name of the dependent variable.}
  \item{x}{Character. Name of the polynomial (I(1)) regressor.}
  \item{panel_id}{Character. Name of the panel identifier variable.}
  \item{time_id}{Character. Name of the time variable.}
  \item{model}{Character. \code{"mg"} for Group-Mean FM-OLS (default) or \code{"pmg"} for Pooled FM-OLS.}
  \item{q}{Integer. Polynomial degree: \code{2} (default) or \code{3}.}
  \item{controls}{Character vector of additional I(1) control variable names.}
  \item{trend}{Integer. Trend type: \code{1} (demeaning only) or \code{2} (demeaning + detrending).}
  \item{kernel}{Character. HAC kernel: \code{"ba"}, \code{"pa"}, \code{"qs"}, \code{"tr"}, \code{"bo"}.}
  \item{bw}{Character or numeric. Bandwidth for HAC estimation.}
  \item{effects}{Character. For Pooled FM-OLS: \code{"oneway"} or \code{"twoways"}.}
  \item{corr_rob}{Logical. For Group-Mean FM-OLS: cross-sectionally robust VCV.}
}
\value{
  An object of class \code{"xtpcmg"} containing FM-OLS coefficient estimates,
  standard errors, and turning point estimates.
}
\references{
Wagner, M. and Reichold, K. (2023). Panel cointegrating polynomial regressions.
\emph{Econometric Theory}, 39(3), 483-536.
}
\examples{
\donttest{
dat <- grunfeld_cmg()
res <- xtpcmg(dat, y = "invest", x = "mvalue",
              panel_id = "firm", time_id = "year",
              model = "mg", q = 2L)
print(res)
}
}
\keyword{models}
