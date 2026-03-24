\name{grunfeld_panel}
\alias{grunfeld_panel}
\title{Grunfeld Investment Panel Dataset}
\description{
  Returns a balanced panel dataset for Grunfeld (1958) firms, used to
  illustrate panel causality tests.
}
\usage{
grunfeld_panel()
}
\value{
  A data frame with columns: \code{firm}, \code{year}, \code{invest}, \code{mvalue}.
}
\references{
Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment}.
PhD thesis, University of Chicago.
}
\examples{
dat <- grunfeld_panel()
head(dat)
}
\keyword{datasets}
