\name{grunfeld_pqroot}
\alias{grunfeld_pqroot}
\title{Grunfeld Investment Data for Panel Unit Root Tests}
\description{
  Returns a balanced panel dataset containing investment and capital stock
  data for Grunfeld (1958) firms, used to illustrate panel unit root tests.
}
\usage{
grunfeld_pqroot()
}
\value{
  A data frame with columns: \code{firm}, \code{year}, \code{invest}, \code{capital}.
}
\references{
Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment}.
PhD thesis, University of Chicago.
}
\examples{
dat <- grunfeld_pqroot()
head(dat)
}
\keyword{datasets}
