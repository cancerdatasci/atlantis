\name{RandomForest-class}
\docType{class}
\alias{RandomForest-class}
\alias{treeresponse,RandomForest-method}
\alias{weights,RandomForest-method}
\alias{where,RandomForest-method}
\alias{show,RandomForest-method}

\title{Class "RandomForest"}
\description{A class for representing random forest ensembles. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RandomForest", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ensemble}:}{Object of class \code{"list"}, each element
                           being an object of class \code{"\linkS4class{BinaryTree}"}.}
    \item{\code{data}:}{ an object of class \code{"\linkS4class{ModelEnv}"}.}
    \item{\code{initweights}:}{ a vector of initial weights.}
    \item{\code{weights}:}{ a list of weights defining the sub-samples.}
    \item{\code{where}:}{ a matrix of integers vectors of length n (number of
                          observations in the learning sample) giving the
                          number of the terminal node the corresponding
                          observations is element of (in each tree).}
    \item{\code{data}:}{ an object of class \code{"\linkS4class{ModelEnv}"}.}
    \item{\code{responses}:}{ an object of class \code{"VariableFrame"}
                              storing the values of the response variable(s). }
    \item{\code{cond_distr_response}:}{ a function computing the conditional
                                        distribution of the response. }
    \item{\code{predict_response}:}{ a function for computing predictions. }
    \item{\code{prediction_weights}:}{ a function for extracting weights from
                                     terminal nodes. }
    \item{\code{get_where}:}{ a function for determining the number
        of terminal nodes observations fall into. }
    \item{\code{update}:}{ a function for updating weights.}
  }
}
\section{Methods}{
  \describe{
    \item{treeresponse}{\code{signature(object = "RandomForest")}: ... }
    \item{weights}{\code{signature(object = "RandomForest")}: ... }
    \item{where}{\code{signature(object = "RandomForest")}: ... }
	 }
}
\examples{

    set.seed(290875)

    ### honest (i.e., out-of-bag) cross-classification of 
    ### true vs. predicted classes
    data("mammoexp", package = "TH.data")
    table(mammoexp$ME, predict(cforest(ME ~ ., data = mammoexp, 
                               control = cforest_unbiased(ntree = 50)), 
                               OOB = TRUE))
}
\keyword{classes}
