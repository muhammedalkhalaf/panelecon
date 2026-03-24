\name{xtfifevd}
\alias{xtfifevd}
\title{Fixed Effects Filtered and Vector Decomposition for Panel Data}
\description{
  Estimates panel data models with time-invariant regressors using the
  Fixed Effects Vector Decomposition (FEVD), Fixed Effects Filtered (FEF),
  or FEF-IV methods. Provides consistent standard errors following
  Pesaran and Zhou (2016).
}
\usage{
xtfifevd(formula, data, index, zinvariants, method = c("fevd", "fef", "fefiv"),
         instruments = NULL, robust = FALSE)
}
\arguments{
  \item{formula}{A formula of the form \code{y ~ x1 + x2}.}
  \item{data}{A data frame in long format.}
  \item{index}{Character vector of length 2: \code{c("id_var", "time_var")}.}
  \item{zinvariants}{Character vector naming the time-invariant regressors.}
  \item{method}{Estimation method: \code{"fevd"}, \code{"fef"}, or \code{"fefiv"}.}
  \item{instruments}{Character vector of instrument names (required for \code{"fefiv"}).}
  \item{robust}{Logical. If \code{TRUE}, use HC robust standard errors. Default \code{FALSE}.}
}
\value{
  An object of class \code{"xtfifevd"} containing coefficient estimates,
  standard errors, and diagnostic statistics.
}
\references{
Pesaran, M.H. and Zhou, Q. (2016). Estimation of Time-Invariant Effects in
Static Panel Data Models. \emph{Econometric Reviews}, 37(10), 1137-1171.

Plumper, T. and Troeger, V.E. (2007). Efficient Estimation of Time-Invariant
and Rarely Changing Variables in Finite Sample Panel Analyses with Unit Fixed
Effects. \emph{Political Analysis}, 15(2), 124-139.
}
\examples{
\donttest{
set.seed(42)
n <- 10; tt <- 15
uid  <- rep(1:n, each = tt)
tval <- rep(1:tt, times = n)
z_i  <- rep(rnorm(n), each = tt)
x_it <- rnorm(n * tt)
y    <- 0.5 * x_it + 0.8 * z_i + rnorm(n * tt, sd = 0.5)
dat  <- data.frame(id = uid, time = tval, y = y, x = x_it, z = z_i)
res <- xtfifevd(y ~ x, data = dat, index = c("id", "time"),
                zinvariants = "z", method = "fevd")
print(res)
}
}
\keyword{models}
