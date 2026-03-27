\name{xtpqroot}
\alias{xtpqroot}
\title{Panel Quantile Unit Root Tests}
\description{
  Tests for a panel unit root using either the CIPS(tau) quantile test
  (Yang, Wei & Cai 2022) or the tFR Fourier-smooth-break test
  (Corakci & Omay 2023).
}
\usage{
xtpqroot(data, var, panel_id, time_id,
         test = c("cipstau", "tfr"),
         model = c("intercept", "trend", "trendshift"),
         quantiles = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
         maxlag = -1L, reps = 500L, bootreps = 500L)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{var}{Character. Name of the variable to test.}
  \item{panel_id}{Character. Name of the panel identifier variable.}
  \item{time_id}{Character. Name of the time variable.}
  \item{test}{Character. \code{"cipstau"} (default) or \code{"tfr"}.}
  \item{model}{Character. Deterministic specification: \code{"intercept"} or \code{"trend"}.}
  \item{quantiles}{Numeric vector. Quantiles for CIPS(tau) test.}
  \item{maxlag}{Integer. Maximum lag (-1 for automatic selection).}
  \item{reps}{Integer. Monte Carlo replications for CIPS(tau) p-values.}
  \item{bootreps}{Integer. Bootstrap replications for tFR p-values.}
}
\value{
  An object of class \code{"xtpqroot"}.
}
\references{
Yang, Z., Wei, Z. and Cai, Y. (2022). Quantile unit root inference for
panel data with common shocks. \emph{Economics Letters}, 219, 110809.

Corakci, A. and Omay, T. (2023). Is there convergence in renewable energy
deployment? \emph{Renewable Energy}, 205, 648-662.
}
\examples{
\donttest{
dat <- grunfeld_pqroot()
res <- xtpqroot(dat, var = "invest", panel_id = "firm",
                time_id = "year", test = "cipstau",
                quantiles = c(0.1, 0.5, 0.9), reps = 100L)
print(res)
}
}
\keyword{htest}
