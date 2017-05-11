\name{reweight}

\alias{reweight}
\alias{reweight.linearModel}
\alias{reweight.glinearModel}

\title{Re-fitting Models with New Weights}

\description{
 Generic function for re-fitting a model object using the same
 observations but different weights.
}

\usage{
reweight(object, weights, \dots)
}

\arguments{
 \item{object}{a fitted model object.}
 \item{weights}{a vector of weights.}
 \item{\dots}{arguments passed to methods.}
}

\details{
  The method is not unsimilar in spirit to \code{\link[stats]{update}}, but
  much more narrowly focused. It should return an updated fitted model 
  derived from re-fitting the model on the same observations but using
  different weights.
}

\value{The re-weighted fitted model object.}

\seealso{\code{\link{update}}}

\examples{
  ## fit cars regression
  mf <- dpp(linearModel, dist ~ speed, data = cars)
  fm <- fit(linearModel, mf)
  fm
  
  ## re-fit, excluding the last 4 observations
  ww <- c(rep(1, 46), rep(0, 4))
  reweight(fm, ww)
}

\keyword{regression}
