\name{cforest}
\alias{cforest}
\alias{proximity}
\title{ Random Forest }
\description{
    An implementation of the random forest and bagging ensemble algorithms
    utilizing conditional inference trees as base learners.
}
\usage{
cforest(formula, data = list(), subset = NULL, weights = NULL, 
        controls = cforest_unbiased(),
        xtrafo = ptrafo, ytrafo = ptrafo, scores = NULL)
proximity(object, newdata = NULL)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. Note 
                  that symbols like \code{:} and \code{-} will not work
                  and the tree will make use of all variables listed on the
                  rhs of \code{formula}.}                  
  \item{data}{ an data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. Non-negative integer valued weights are
                  allowed as well as non-negative real weights.
                  Observations are sampled (with or without replacement)
                  according to probabilities \code{weights / sum(weights)}.
                  The fraction of observations to be sampled (without replacement)
                  is computed based on the sum of the weights if all weights
                  are integer-valued and based on the number of weights greater zero
                  else. Alternatively, \code{weights} can be a double matrix defining
                  case weights for all \code{ncol(weights)} trees in the forest directly.
                  This requires more storage but gives the user more control.}
  \item{controls}{an object of class \code{\link{ForestControl-class}}, which can be
                  obtained using \code{\link{cforest_control}} (and its
		  convenience interfaces \code{cforest_unbiased} and \code{cforest_classical}).}
  \item{xtrafo}{ a function to be applied to all input variables.
                By default, the \code{\link{ptrafo}} function is applied.}
  \item{ytrafo}{ a function to be applied to all response variables.
                By default, the \code{\link{ptrafo}} function is applied.}
  \item{scores}{ an optional named list of scores to be attached to ordered
               factors.}
  \item{object}{ an object as returned by \code{cforest}.}
  \item{newdata}{ an optional data frame containing test data.}
}
\details{

  This implementation of the random forest (and bagging) algorithm differs
  from the reference implementation in \code{\link[randomForest]{randomForest}}
  with respect to the base learners used and the aggregation scheme applied.
  
  Conditional inference trees, see \code{\link{ctree}}, are fitted to each
  of the \code{ntree} (defined via \code{\link{cforest_control}}) 
  bootstrap samples of the learning sample. Most of the hyper parameters in 
  \code{\link{cforest_control}} regulate the construction of the conditional inference trees.
  Therefore you MUST NOT change anything you don't understand completely.

  Hyper parameters you might want to change in \code{\link{cforest_control}} are: 

  1. The number of randomly preselected variables \code{mtry}, which is fixed 
  to the value 5 by default here for technical reasons, while in 
  \code{\link[randomForest]{randomForest}} the default values for classification and regression
  vary with the number of input variables. 

  2. The number of trees \code{ntree}. Use more trees if you have more variables.

  3. The depth of the trees, regulated by \code{mincriterion}. Usually unstopped and unpruned
  trees are used in random forests. To grow large trees, set \code{mincriterion} to a small value.

  The aggregation scheme works by averaging observation weights extracted
  from each of the \code{ntree} trees and NOT by averaging predictions directly
  as in \code{\link[randomForest]{randomForest}}.
  See Hothorn et al. (2004) for a description. 

  Predictions can be computed using \code{\link{predict}}. For observations
  with zero weights, predictions are computed from the fitted tree 
  when \code{newdata = NULL}. While \code{\link{predict}} returns predictions
  of the same type as the response in the data set by default (i.e., predicted class labels for factors), 
  \code{\link{treeresponse}} returns the statistics of the conditional distribution of the response
  (i.e., predicted class probabilities for factors). The same is done by \code{predict(..., type = "prob")}.
  Note that for multivariate responses \code{predict} does not convert predictions to the type 
  of the response, i.e., \code{type = "prob"} is used.

  Ensembles of conditional inference trees have not yet been extensively
  tested, so this routine is meant for the expert user only and its current
  state is rather experimental. However, there are some things available 
  in \code{\link{cforest}} that can't be done with \code{\link[randomForest]{randomForest}}, 
  for example fitting forests to censored response variables (see Hothorn et al., 2006a) or to
  multivariate and ordered responses.
  
  Moreover, when predictors vary in their scale of measurement of number 
  of categories, variable selection and computation of variable importance is biased 
  in favor of variables with many potential cutpoints in \code{\link[randomForest]{randomForest}}, 
  while in \code{\link{cforest}} unbiased trees and an adequate resampling scheme 
  are used by default. See Hothorn et al. (2006b) and Strobl et al. (2007) 
  as well as Strobl et al. (2009). 

  The \code{proximity} matrix is an \eqn{n \times n} matrix \eqn{P} with \eqn{P_{ij}}
  equal to the fraction of trees where observations \eqn{i} and \eqn{j} 
  are element of the same terminal node (when both \eqn{i} and \eqn{j}
  had non-zero weights in the same bootstrap sample).

}
\value{
  An object of class \code{\link{RandomForest-class}}.
}
\references{ 

    Leo Breiman (2001). Random Forests. \emph{Machine Learning}, 45(1), 5--32.

    Torsten Hothorn, Berthold Lausen, Axel Benner and Martin Radespiel-Troeger
    (2004). Bagging Survival Trees. \emph{Statistics in Medicine}, \bold{23}(1), 77--91.

    Torsten Hothorn, Peter Buhlmann, Sandrine Dudoit, Annette Molinaro 
    and Mark J. van der Laan (2006a). Survival Ensembles. \emph{Biostatistics}, 
    \bold{7}(3), 355--373.

    Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006b). Unbiased
    Recursive Partitioning: A Conditional Inference Framework.
    \emph{Journal of Computational and Graphical Statistics}, \bold{15}(3),
    651--674.  Preprint available from 
    \url{http://statmath.wu-wien.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf}

    Carolin Strobl, Anne-Laure Boulesteix, Achim Zeileis and Torsten Hothorn (2007).
    Bias in Random Forest Variable Importance Measures: Illustrations, Sources and 
    a Solution. \emph{BMC Bioinformatics}, \bold{8}, 25. 
    \url{http://www.biomedcentral.com/1471-2105/8/25}

    Carolin Strobl, James Malley and Gerhard Tutz (2009).
    An Introduction to Recursive Partitioning: Rationale, Application, and Characteristics of
    Classification and Regression Trees, Bagging, and Random forests.
    \emph{Psychological Methods}, \bold{14}(4), 323--348.

}
\examples{

    set.seed(290875)

    ### honest (i.e., out-of-bag) cross-classification of
    ### true vs. predicted classes
    data("mammoexp", package = "TH.data")
    table(mammoexp$ME, predict(cforest(ME ~ ., data = mammoexp, 
                               control = cforest_unbiased(ntree = 50)),
                               OOB = TRUE))

    ### fit forest to censored response
    if (require("TH.data") && require("survival")) {

        data("GBSG2", package = "TH.data")
        bst <- cforest(Surv(time, cens) ~ ., data = GBSG2, 
                   control = cforest_unbiased(ntree = 50))

        ### estimate conditional Kaplan-Meier curves
        treeresponse(bst, newdata = GBSG2[1:2,], OOB = TRUE)

        ### if you can't resist to look at individual trees ...
        party:::prettytree(bst@ensemble[[1]], names(bst@data@get("input")))
    }

    ### proximity, see ?randomForest
    iris.cf <- cforest(Species ~ ., data = iris, 
                       control = cforest_unbiased(mtry = 2))
    iris.mds <- cmdscale(1 - proximity(iris.cf), eig = TRUE)
    op <- par(pty="s")
    pairs(cbind(iris[,1:4], iris.mds$points), cex = 0.6, gap = 0, 
          col = c("red", "green", "blue")[as.numeric(iris$Species)],
          main = "Iris Data: Predictors and MDS of Proximity Based on cforest")
    par(op)

}
\keyword{tree}
