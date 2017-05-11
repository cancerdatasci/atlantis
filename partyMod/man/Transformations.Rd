\name{Transformations}
\alias{ptrafo}
\alias{ff_trafo}
\title{ Function for Data Transformations }
\description{
  Transformations of Response or Input Variables
}
\usage{
ptrafo(data, numeric_trafo = id_trafo, factor_trafo = ff_trafo, 
    ordered_trafo = of_trafo, surv_trafo = logrank_trafo, 
    var_trafo = NULL)
ff_trafo(x)
}
\arguments{
  \item{data}{an object of class \code{data.frame}.}
  \item{numeric_trafo}{a function to by applied to \code{numeric} 
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns.}
  \item{ordered_trafo}{a function to by applied to \code{ordered}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns (usually some scores).}
  \item{factor_trafo}{a function to by applied to \code{factor}
    elements of \code{data} returning a matrix with \code{nrow(data)}
    rows and an arbitrary number of columns (usually a dummy or contrast 
    matrix).}
  \item{surv_trafo}{a function to by applied to 
    elements of class \code{Surv} of \code{data} returning a 
    matrix with \code{nrow(data)} rows and an arbitrary number of columns.}
  \item{var_trafo}{an optional named list of functions to be applied to the
                   corresponding variables in \code{data}.}
  \item{x}{ a factor}
}
\details{

  \code{trafo} applies its arguments to the elements of \code{data}
  according to the classes of the elements. See \code{\link[coin]{Transformations}}
  for more documentation and examples.

  In the presence of missing values, one needs to make sure that all
  user-supplied functions deal with that. 

}
\value{

  A named matrix with \code{nrow(data)} rows and 
  arbitrary number of columns.

}
\examples{

  ### rank a variable
  ptrafo(data.frame(y = 1:20), 
         numeric_trafo = function(x) rank(x, na.last = "keep"))

  ### dummy coding of a factor
  ptrafo(data.frame(y = gl(3, 9)))

}
\keyword{manip}
