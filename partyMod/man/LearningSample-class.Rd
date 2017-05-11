\name{LearningSample Class}
\docType{class}
\alias{LearningSample-class}

\title{Class "LearningSample"}
\description{  
    Objects of this class represent data for fitting tree-based models.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LearningSample", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{responses}:}{Object of class \code{"VariableFrame"} with the
                             response variables. }
    \item{\code{inputs}:}{Object of class \code{"VariableFrame"} with the
                          input variables.}
    \item{\code{weights}:}{Object of class \code{"numeric"}, a vector of
                           case counts or weights. }
    \item{\code{nobs}:}{Object of class \code{"integer"}, the number of
                        observations. }
    \item{\code{ninputs}:}{Object of class \code{"integer"}, the number of
                           input variables.}
  }
}
\section{Methods}{
No methods defined with class "LearningSample" in the signature.
}
\keyword{classes}
