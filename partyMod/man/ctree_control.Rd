\name{Control ctree Hyper Parameters}
\alias{ctree_control}
\title{ Control for Conditional Inference Trees }
\description{

  Various parameters that control aspects of the `ctree' fit.

}
\usage{
ctree_control(teststat = c("quad", "max"), 
              testtype = c("Bonferroni", "MonteCarlo", 
                           "Univariate", "Teststatistic"), 
              mincriterion = 0.95, minsplit = 20, minbucket = 7, 
              stump = FALSE, nresample = 9999, maxsurrogate = 0, 
              mtry = 0, savesplitstats = TRUE, maxdepth = 0)
}
\arguments{
  \item{teststat}{ a character specifying the type of the test statistic
                       to be applied. }
  \item{testtype}{ a character specifying how to compute the distribution of
                   the test statistic. }
  \item{mincriterion}{ the value of the test statistic (for \code{testtype == "Teststatistic"}),
                       or 1 - p-value (for other values of \code{testtype}) that
                       must be exceeded in order to implement a split. }
  \item{minsplit}{ the minimum sum of weights in a node in order to be considered
                   for splitting. }
  \item{minbucket}{ the minimum sum of weights in a terminal node. }
  \item{stump}{ a logical determining whether a stump (a tree with three
                nodes only) is to be computed. }
  \item{nresample}{ number of Monte-Carlo replications to use when the
                    distribution of the test statistic is simulated.}
  \item{maxsurrogate}{ number of surrogate splits to evaluate. Note the
                       currently only surrogate splits in ordered
                       covariables are implemented. }
  \item{mtry}{ number of input variables randomly sampled as candidates 
               at each node for random forest like algorithms. The default
               \code{mtry = 0} means that no random selection takes place.}
  \item{savesplitstats}{ a logical determining if the process of standardized
                         two-sample statistics for split point estimate
                         is saved for each primary split.}
  \item{maxdepth}{ maximum depth of the tree. The default \code{maxdepth = 0}
                   means that no restrictions are applied to tree sizes.}
}
\details{

  The arguments \code{teststat}, \code{testtype} and \code{mincriterion}
  determine how the global null hypothesis of independence between all input
  variables and the response is tested (see \code{\link{ctree}}). The 
  argument \code{nresample} is the number of Monte-Carlo replications to be
  used when \code{testtype = "MonteCarlo"}.

  A split is established when the sum of the weights in both daugther nodes
  is larger than \code{minsplit}, this avoids pathological splits at the
  borders. When \code{stump = TRUE}, a tree with at most two terminal nodes
  is computed.

  The argument \code{mtry > 0} means that a random forest like `variable
  selection', i.e., a random selection of \code{mtry} input variables, is
  performed in each node.

  It might be informative to look at scatterplots of input variables against
  the standardized two-sample split statistics, those are available when
  \code{savesplitstats = TRUE}. Each node is then associated with a vector
  whose length is determined by the number of observations in the learning
  sample and thus much more memory is required.

}
\value{
  An object of class \code{\link{TreeControl}}.
}
\keyword{misc}
