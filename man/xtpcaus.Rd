\name{xtpcaus}
\alias{xtpcaus}
\title{Panel Granger Causality Tests}
\description{
  Tests whether \code{x} Granger-causes \code{y} in a balanced panel using
  either the Panel Fourier Toda-Yamamoto (PFTY) test or the Panel Quantile
  Causality (PQC) test.
}
\usage{
xtpcaus(data, y, x, panel_id, time_id,
        test = c("pfty", "pqc"),
        pmax = 4L, dmax = 1L, nboot = 499L, kmax = 3L,
        ic = c("aic", "bic"),
        quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9), seed = -1L)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{y}{Character. Name of the dependent (caused) variable.}
  \item{x}{Character. Name of the independent (causing) variable.}
  \item{panel_id}{Character. Name of the panel identifier variable.}
  \item{time_id}{Character. Name of the time variable.}
  \item{test}{Character. \code{"pfty"} (default) or \code{"pqc"}.}
  \item{pmax}{Integer. Maximum lag order. Default 4.}
  \item{dmax}{Integer. Maximum integration order for Toda-Yamamoto. Default 1.}
  \item{nboot}{Integer. Bootstrap replications. Default 499.}
  \item{kmax}{Integer. Maximum Fourier frequency (PFTY only). Default 3.}
  \item{ic}{Character. Information criterion: \code{"aic"} or \code{"bic"}.}
  \item{quantiles}{Numeric vector. Quantile grid for PQC test.}
  \item{seed}{Integer. Random seed (-1 = no seed).}
}
\value{
  An object of class \code{"xtpcaus"} containing test statistics and p-values.
}
\references{
Dumitrescu, E.I. and Hurlin, C. (2012). Testing for Granger non-causality in
heterogeneous panels. \emph{Economic Modelling}, 29(4), 1450-1460.
}
\examples{
\donttest{
dat <- grunfeld_panel()
res <- xtpcaus(dat, y = "invest", x = "mvalue",
               panel_id = "firm", time_id = "year",
               test = "pfty", pmax = 2L, nboot = 99L)
print(res)
}
}
\keyword{htest}
