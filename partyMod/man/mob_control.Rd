\name{mob_control}
\alias{mob_control}

\title{Control Parameters for Model-based Partitioning}

\description{
  Various parameters that control aspects the fitting algorithm
  for recursively partitioned \code{\link{mob}} models.
}

\usage{
mob_control(alpha = 0.05, bonferroni = TRUE, minsplit = 20, trim = 0.1,
  objfun = deviance, breakties = FALSE, parm = NULL, verbose = FALSE)
}

\arguments{
  \item{alpha}{numeric significance level. A node is splitted when
    the (possibly Bonferroni-corrected) \eqn{p} value for any parameter
    stability test in that node falls below \code{alpha}.}
  \item{bonferroni}{logical. Should \eqn{p} values be Bonferroni corrected?}
  \item{minsplit}{integer. The minimum number of observations (sum of the
    weights) in a node.}
  \item{trim}{numeric. This specifies the trimming in the parameter instability
    test for the numerical variables. If smaller than 1, it is interpreted
    as the fraction relative to the current node size.}
  \item{objfun}{function. A function for extracting the minimized value of
    the objective function from a fitted model in a node.}
  \item{breakties}{logical. Should ties in numeric variables be broken
    randomly for computing the associated parameter instability test?}
  \item{parm}{numeric or character. Number or name of model parameters
    included in the parameter instability tests (by default all parameters
    are included).}
  \item{verbose}{logical. Should information about the fitting process
    of \code{\link{mob}} (such as test statistics, \eqn{p} values, selected
    splitting variables and split points) be printed to the screen?}
}

\details{
  See \code{\link{mob}} for more details and references.
}

\seealso{\code{\link{mob}}}

\value{
  A list of class \code{mob_control} containing the control parameters.
}

\keyword{misc}
