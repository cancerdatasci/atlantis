\name{BinaryTree Class}
\docType{class}
\alias{BinaryTree-class}
\alias{weights}
\alias{weights-methods}
\alias{weights,BinaryTree-method}
\alias{show,BinaryTree-method}
\alias{where}
\alias{where-methods}
\alias{where,BinaryTree-method}
\alias{response}
\alias{response-methods}
\alias{response,BinaryTree-method}
\alias{nodes}
\alias{nodes-methods}
\alias{nodes,BinaryTree,integer-method}
\alias{nodes,BinaryTree,numeric-method}
\alias{treeresponse}
\alias{treeresponse-methods}
\alias{treeresponse,BinaryTree-method}
\title{Class "BinaryTree"}
\description{A class for representing binary trees.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BinaryTree", ...)}.
The most important slot is \code{tree}, a (recursive) list with elements
\describe{
  \item{nodeID}{ an integer giving the number of the node, starting with
                \code{1} in the root node.}
  \item{weights}{ the case weights (of the learning sample) corresponding to
                this node.}
  \item{criterion}{ a list with test statistics and p-values for each partial
                  hypothesis.}
  \item{terminal}{ a logical specifying if this is a terminal node.}
  \item{psplit}{ primary split: a list with elements \code{variableID} (the
               number of the input variable splitted), \code{ordered} (a
               logical whether the input variable is ordered),
               \code{splitpoint} (the cutpoint or set of levels to the left),
               \code{splitstatistics} saves the process of standardized 
               two-sample statistics the split point estimation is based on.
               The logical \code{toleft} determines if observations
               go left or right down the tree. For nominal splits, the slot
               \code{table} is a vector being greater zero if the
               corresponding level is available in the corresponding node.}
  \item{ssplits}{ a list of surrogate splits, each with the same elements as
                 \code{psplit}.}
  \item{prediction}{ the prediction of the node: the mean for numeric
                     responses and the conditional class probabilities for 
                     nominal or ordered respones. For censored responses,
                     this is the mean of the logrank scores and useless as
                     such.}
  \item{left}{ a list representing the left daughter node. }
  \item{right}{ a list representing the right daugther node.}
}

Please note that this data structure may be subject to change in future
releases of the package.

}
\section{Slots}{
  \describe{
    \item{\code{data}:}{ an object of class \code{"\linkS4class{ModelEnv}"}.}
    \item{\code{responses}:}{ an object of class \code{"VariableFrame"}
                              storing the values of the response variable(s). }
    \item{\code{cond_distr_response}:}{ a function computing the conditional
                                        distribution of the response. } 
    \item{\code{predict_response}:}{ a function for computing predictions. }
    \item{\code{tree}:}{ a recursive list representing the tree. See above. }
    \item{\code{where}:}{ an integer vector of length n (number of
                          observations in the learning sample) giving the
                          number of the terminal node the corresponding
                          observations is element of. }
    \item{\code{prediction_weights}:}{ a function for extracting weights from
                                     terminal nodes. }
    \item{\code{get_where}:}{ a function for determining the number
        of terminal nodes observations fall into. }
    \item{\code{update}:}{ a function for updating weights.}
  }
}
\section{Extends}{
Class \code{"BinaryTreePartition"}, directly.
}
\section{Methods}{
\describe{
    \item{\code{response(object, ...)}:}{extract the response variables the
       tree was fitted to.}
    \item{\code{treeresponse(object, newdata = NULL, ...)}:}{compute
      statistics for the conditional distribution of the response as
      modelled by the tree. For regression problems, this is just the mean.
      For nominal or ordered responses, estimated conditional class
      probabilities are returned. Kaplan-Meier curves are computed for
      censored responses. Note that a list with one element for each
      observation is returned.}
    \item{\code{Predict(object, newdata = NULL, ...)}:}{ compute predictions.}
    \item{\code{weights(object, newdata = NULL, ...)}:}{ extract the weight
      vector from terminal nodes each element of the learning sample is
      element of (\code{newdata = NULL}) and for new observations, 
      respectively.}
    \item{\code{where(object, newdata = NULL, ...)}:}{ extract the number of 
      the terminal nodes each element of the learning sample is
      element of (\code{newdata = NULL}) and for new observations, 
      respectively.}
    \item{\code{nodes(object, where, ...)}:}{ extract the nodes with
      given number (\code{where}).}
    \item{\code{plot(x, ...)}:}{ a plot method for \code{BinaryTree}
      objects, see \code{\link{plot.BinaryTree}}.}
    \item{\code{print(x, ...)}:}{ a print method for \code{BinaryTree}
      objects.}
}
}
\examples{

  set.seed(290875)

  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq,   
                 controls = ctree_control(maxsurrogate = 3))

  ### distribution of responses in the terminal nodes
  plot(airq$Ozone ~ as.factor(where(airct)))

  ### get all terminal nodes from the tree
  nodes(airct, unique(where(airct)))

  ### extract weights and compute predictions
  pmean <- sapply(weights(airct), function(w) weighted.mean(airq$Ozone, w))

  ### the same as
  drop(Predict(airct))

  ### or
  unlist(treeresponse(airct))

  ### don't use the mean but the median as prediction in each terminal node
  pmedian <- sapply(weights(airct), function(w) 
                 median(airq$Ozone[rep(1:nrow(airq), w)]))

  plot(airq$Ozone, pmean, col = "red")
  points(airq$Ozone, pmedian, col = "blue")
}
\keyword{classes}
