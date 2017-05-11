\name{ForestControl-class}
\docType{class}
\alias{ForestControl-class}

\title{Class "ForestControl" }
\description{
  Objects of this class represent the hyper parameter setting for forest
  growing.
}
\section{Objects from the Class}{
Objects can be created by \code{\link{cforest_control}}.
}
\section{Slots}{
	 \describe{
    \item{\code{ntree}:}{number of trees in the forest.}
    \item{\code{replace}:}{sampling with or without replacement.}
    \item{\code{fraction}:}{fraction of observations to sample without replacement.}
    \item{\code{trace}:}{logical indicating if a progress bar shall be printed.}
    \item{\code{varctrl}:}{Object of class \code{"VariableControl"}}
    \item{\code{splitctrl}:}{Object of class \code{"SplitControl"}}
    \item{\code{gtctrl}:}{Object of class \code{"GlobalTestControl"}}
    \item{\code{tgctrl}:}{Object of class \code{"TreeGrowControl"}}
  }
}
\section{Extends}{
Class \code{"TreeControl"}, directly.
}
\section{Methods}{
No methods defined with class "ForestControl" in the signature.
}
\keyword{classes}
