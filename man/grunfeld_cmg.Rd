\name{grunfeld_cmg}
\alias{grunfeld_cmg}
\title{Grunfeld Investment Panel Dataset for FM-OLS}
\description{
  Returns a balanced panel dataset for Grunfeld (1958) firms.
}
\usage{
grunfeld_cmg()
}
\value{
  A data frame with columns: \code{firm}, \code{year}, \code{invest}, \code{mvalue}, \code{capital}.
}
\references{
Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment}.
PhD thesis, University of Chicago.
}
\examples{
dat <- grunfeld_cmg()
head(dat)
}
\keyword{datasets}
