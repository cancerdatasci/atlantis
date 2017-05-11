\name{Conditional Inference Trees}
\alias{ctree}
\alias{conditionalTree}
\title{ Conditional Inference Trees }
\description{
  Recursive partitioning for continuous, censored, ordered, nominal and
  multivariate response variables in a conditional inference framework. 
}
\usage{
ctree(formula, data, subset = NULL, weights = NULL, 
      controls = ctree_control(), xtrafo = ptrafo, ytrafo = ptrafo, 
      scores = NULL)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. Note
                  that symbols like \code{:} and \code{-} will not work
                  and the tree will make use of all variables listed on the
                  rhs of \code{formula}.}
  \item{data}{ a data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. Only non-negative integer valued weights are
                  allowed.}
  \item{controls}{an object of class \code{\link{TreeControl}}, which can be
                  obtained using \code{\link{ctree_control}}.}
  \item{xtrafo}{ a function to be applied to all input variables.
                By default, the \code{\link{ptrafo}} function is applied.}
  \item{ytrafo}{ a function to be applied to all response variables. 
                By default, the \code{\link{ptrafo}} function is applied.}
  \item{scores}{ an optional named list of scores to be attached to ordered
               factors.}
}
\details{

  Conditional inference trees estimate a regression relationship by binary recursive
  partitioning in a conditional inference framework. Roughly, the algorithm
  works as follows: 1) Test the global null hypothesis of independence between
  any of the input variables and the response (which may be multivariate as well). 
  Stop if this hypothesis cannot be rejected. Otherwise select the input
  variable with strongest association to the resonse. This
  association is measured by a p-value corresponding to a test for the
  partial null hypothesis of a single input variable and the response.
  2) Implement a binary split in the selected input variable. 
  3) Recursively repeate steps 1) and 2). 

  The implementation utilizes a unified framework for conditional inference,
  or permutation tests, developed by Strasser and Weber (1999). The stop
  criterion in step 1) is either based on multiplicity adjusted p-values 
  (\code{testtype == "Bonferroni"}
  or \code{testtype == "MonteCarlo"} in \code{\link{ctree_control}}),
  on the univariate p-values (\code{testtype == "Univariate"}),
  or on values of the test statistic
  (\code{testtype == "Teststatistic"}). In both cases, the
  criterion is maximized, i.e., 1 - p-value is used. A split is implemented 
  when the criterion exceeds the value given by \code{mincriterion} as
  specified in \code{\link{ctree_control}}. For example, when 
  \code{mincriterion = 0.95}, the p-value must be smaller than
  $0.05$ in order to split this node. This statistical approach ensures that
  the right sized tree is grown and no form of pruning or cross-validation
  or whatsoever is needed. The selection of the input variable to split in
  is based on the univariate p-values avoiding a variable selection bias
  towards input variables with many possible cutpoints.

  Multiplicity-adjusted Monte-Carlo p-values are computed 
  following a "min-p" approach. The univariate 
  p-values based on the limiting distribution (chi-square
  or normal) are computed for each of the random 
  permutations of the data. This means that one should
  use a quadratic test statistic when factors are in
  play (because the evaluation of the corresponding
  multivariate normal distribution is time-consuming).

  By default, the scores for each ordinal factor \code{x} are
  \code{1:length(x)}, this may be changed using \code{scores = list(x =
  c(1,5,6))}, for example.

  Predictions can be computed using \code{\link{predict}} or
  \code{\link{treeresponse}}.  The first function accepts arguments
  \code{type = c("response", "node", "prob")} where \code{type = "response"}
  returns predicted means, predicted classes or median predicted survival
  times, \code{type = "node"} returns terminal node IDs (identical to
  \code{\link{where}}) and \code{type = "prob"} gives more information about
  the conditional distribution of the response, i.e., class probabilities or
  predicted Kaplan-Meier curves and is identical to
  \code{\link{treeresponse}}.  For observations with zero weights,
  predictions are computed from the fitted tree when \code{newdata = NULL}.

  For a general description of the methodology see Hothorn, Hornik and
  Zeileis (2006) and Hothorn, Hornik, van de Wiel and Zeileis (2006). 
  Introductions for novices can be found in Strobl et al. (2009) and
  at \url{http://github.com/christophM/overview-ctrees.git}.

}
\value{
  An object of class \code{\link{BinaryTree-class}}.
}
\references{ 

   Helmut Strasser and Christian Weber (1999). On the asymptotic theory of permutation
   statistics. \emph{Mathematical Methods of Statistics}, \bold{8}, 220--250.

   Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel and Achim Zeileis (2006).
   A Lego System for Conditional Inference. \emph{The American Statistician},
   \bold{60}(3), 257--263.

   Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased Recursive
   Partitioning: A Conditional Inference Framework. \emph{Journal of
   Computational and Graphical Statistics}, \bold{15}(3), 651--674. 
   Preprint available
   from \url{http://statmath.wu-wien.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf}

   Carolin Strobl, James Malley and Gerhard Tutz (2009).
   An Introduction to Recursive Partitioning: Rationale, Application, and Characteristics of 
   Classification and Regression Trees, Bagging, and Random forests.
   \emph{Psychological Methods}, \bold{14}(4), 323--348. 

}
\examples{

    set.seed(290875)

    ### regression
    airq <- subset(airquality, !is.na(Ozone))
    airct <- ctree(Ozone ~ ., data = airq, 
                   controls = ctree_control(maxsurrogate = 3))
    airct
    plot(airct)
    mean((airq$Ozone - predict(airct))^2)
    ### extract terminal node ID, two ways
    all.equal(predict(airct, type = "node"), where(airct))

    ### classification
    irisct <- ctree(Species ~ .,data = iris)
    irisct
    plot(irisct)
    table(predict(irisct), iris$Species)

    ### estimated class probabilities, a list
    tr <- treeresponse(irisct, newdata = iris[1:10,])

    ### ordinal regression
    data("mammoexp", package = "TH.data")
    mammoct <- ctree(ME ~ ., data = mammoexp) 
    plot(mammoct)

    ### estimated class probabilities
    treeresponse(mammoct, newdata = mammoexp[1:10,])

    ### survival analysis
    if (require("TH.data") && require("survival")) {
        data("GBSG2", package = "TH.data")
        GBSG2ct <- ctree(Surv(time, cens) ~ .,data = GBSG2)
        plot(GBSG2ct)
        treeresponse(GBSG2ct, newdata = GBSG2[1:2,])        
    }

    ### if you are interested in the internals:
    ### generate doxygen documentation
    \dontrun{

        ### download src package into temp dir
        tmpdir <- tempdir()
        tgz <- download.packages("party", destdir = tmpdir)[2]
        ### extract
        untar(tgz, exdir = tmpdir)
        wd <- setwd(file.path(tmpdir, "party"))
        ### run doxygen (assuming it is there)
        system("doxygen inst/doxygen.cfg")
        setwd(wd)
        ### have fun
        browseURL(file.path(tmpdir, "party", "inst", 
                            "documentation", "html", "index.html")) 
    }
}
\keyword{tree}
