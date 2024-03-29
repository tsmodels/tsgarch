% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{estimate.tsgarch.spec}
\alias{estimate.tsgarch.spec}
\alias{estimate}
\title{Estimates an GARCH model given a specification object using maximum likelihood and autodiff}
\usage{
\method{estimate}{tsgarch.spec}(
  object,
  solver = "nloptr",
  control = NULL,
  stationarity_constraint = 0.999,
  ...
)
}
\arguments{
\item{object}{an object of class tsgarch.spec.}

\item{solver}{only \dQuote{nloptr} is currently supported (see \code{\link[nloptr]{nloptr}}).}

\item{control}{solver control parameters.}

\item{stationarity_constraint}{the bound on the inequality constraint for ensuring
the stationary of the GARCH process (see details).}

\item{...}{not currently used.}
}
\value{
An object of class \dQuote{tsgarch.estimate}.
}
\description{
Estimates an GARCH model given a specification object using maximum likelihood and autodiff
}
\details{
The underlying code is written using the TMB framework which uses
automatic differentiation and hence allows the generation of analytic
derivatives.
Stationarity is usually based on the condition that the persistence of the model
is less than 1. The argument \dQuote{stationarity_constraint} allows to fine tune
this. For example, setting it to a very high value will effectively render
this constraint inactive. The default of 0.999 has been found to be a reasonable
bound since values close to one may lead to problems.
Since the nloptr solver make use of analytic Jacobians for the inequality constraint,
these are either provided in closed form or calculated as part of the automatic
differentiation algorithms implemented in the package.
The estimation makes 2 passes to the solver. The first pass uses no parameter
scaling, whilst in the second pass the parameters (as well as bounds) are scaled
making use of the estimated hessian from the first pass in order to generate
a hopefully more robust solution.
}
\author{
Alexios Galanos
}
