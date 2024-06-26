% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multispec.R
\name{+.tsgarch.spec}
\alias{+.tsgarch.spec}
\title{Combine univariate GARCH specifications into a multi-specification object}
\usage{
\method{+}{tsgarch.spec}(x, y)
}
\arguments{
\item{x}{an object of class \dQuote{tsgarch.spec}}

\item{y}{an object of class \dQuote{tsgarch.spec}}
}
\value{
an object of class \dQuote{tsgarch.multispec}
}
\description{
Combine univariate GARCH specifications into a multi-specification object
}
\details{
A simple method for combining multiple specifications into an object
which can then be estimated using parallel resources. Note that the returned
object is effectively a validated list of specification objects with no names.
Names can be assigned post-construction (see example).
}
\examples{
library(xts)
x <- xts(rnorm(1000), as.Date(1:1000, origin = "1970-01-01"))
y <- xts(rnorm(1000), as.Date(1:1000, origin = "1970-01-01"))
z <- xts(rnorm(1000), as.Date(1:1000, origin = "1970-01-01"))
mspec <- garch_modelspec(x, model = "egarch") +
garch_modelspec(y, model = "cgarch") +
garch_modelspec(z, model = "aparch")
names(mspec) <- c("x", "y", "z")
sapply(mspec, function(x) x$model$model)
}
