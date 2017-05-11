\name{Control Forest Hyper Parameters}
\alias{cforest_control}
\alias{cforest_classical}
\alias{cforest_unbiased}
\title{ Control for Conditional Tree Forests }
\description{

  Various parameters that control aspects of the `cforest' fit via
  its `control' argument.

}
\usage{
cforest_unbiased(\dots)
cforest_classical(\dots)
cforest_control(teststat = "max",
                testtype = "Teststatistic",
                mincriterion = qnorm(0.9),
                savesplitstats = FALSE,
                ntree = 500, mtry = 5, replace = TRUE,
                fraction = 0.632, trace = FALSE, \dots)
}
\arguments{
  \item{teststat}{ a character specifying the type of the test statistic
                       to be applied. }
  \item{testtype}{ a character specifying how to compute the distribution of
                   the test statistic. }
  \item{mincriterion}{ the value of the test statistic (for \code{testtype == "Teststatistic"}),
                       or 1 - p-value (for other values of \code{testtype}) that
                       must be exceeded in order to implement a split. }
  \item{mtry}{ number of input variables randomly sampled as candidates 
               at each node for random forest like algorithms. Bagging, as special case 
               of a random forest without random input variable sampling, can 
               be performed by setting \code{mtry} either equal to \code{NULL} or 
               manually equal to the number of input variables.}
  \item{savesplitstats}{ a logical determining whether the process of standardized
                         two-sample statistics for split point estimate
                         is saved for each primary split.}
  \item{ntree}{ number of trees to grow in a forest.}
  \item{replace}{ a logical indicating whether sampling of observations is 
                 done with or without replacement.}
  \item{fraction}{ fraction of number of observations to draw without 
                   replacement (only relevant if \code{replace = FALSE}).}
  \item{trace}{ a logical indicating if a progress bar shall be printed
                while the forest grows.}
  \item{\dots}{ additional arguments to be passed to 
                \code{\link{ctree_control}}.}
}
\details{

  All three functions return an object of class \code{\link{ForestControl-class}}
  defining hyper parameters to be specified via the \code{control} argument
  of \code{\link{cforest}}.

  The arguments \code{teststat}, \code{testtype} and \code{mincriterion}
  determine how the global null hypothesis of independence between all input
  variables and the response is tested (see \code{\link{ctree}}). The 
  argument \code{nresample} is the number of Monte-Carlo replications to be
  used when \code{testtype = "MonteCarlo"}.

  A split is established when the sum of the weights in both daugther nodes
  is larger than \code{minsplit}, this avoids pathological splits at the
  borders. When \code{stump = TRUE}, a tree with at most two terminal nodes
  is computed.

  The \code{mtry} argument regulates a random selection of \code{mtry} input 
  variables in each node. Note that here \code{mtry} is fixed to the value 5 by 
  default for merely technical reasons, while in \code{\link[randomForest]{randomForest}} 
  the default values for classification and regression vary with the number of input 
  variables. Make sure that \code{mtry} is defined properly before using \code{cforest}.

  It might be informative to look at scatterplots of input variables against
  the standardized two-sample split statistics, those are available when
  \code{savesplitstats = TRUE}. Each node is then associated with a vector
  whose length is determined by the number of observations in the learning
  sample and thus much more memory is required.

  The number of trees \code{ntree} can be increased for large numbers of input variables.

  Function \code{cforest_unbiased} returns the settings suggested 
  for the construction of unbiased random forests (\code{teststat = "quad", testtype = "Univ", 
    replace = FALSE}) by Strobl et al. (2007)
  and is the default since version 0.9-90.
  Hyper parameter settings mimicing the behaviour of
  \code{\link[randomForest]{randomForest}} are available in
  \code{cforest_classical} which have been used as default up to
  version 0.9-14. 

  Please note that \code{\link{cforest}}, in contrast to 
  \code{\link[randomForest]{randomForest}}, doesn't grow trees of
  maximal depth. To grow large trees, set \code{mincriterion = 0}.

}
\value{
  An object of class \code{\link{ForestControl-class}}.
}
\references{

    Carolin Strobl, Anne-Laure Boulesteix, Achim Zeileis and Torsten Hothorn (2007).
    Bias in Random Forest Variable Importance Measures: Illustrations, Sources and  
    a Solution. \emph{BMC Bioinformatics}, \bold{8}, 25. 
    \url{http://www.BioMedCentral.com/1471-2105/8/25/}
}
\keyword{misc}
