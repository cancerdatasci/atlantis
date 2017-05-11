\name{varimp}
\alias{varimp}
\alias{varimpAUC}
\title{ Variable Importance }
\description{
    Standard and conditional variable importance for `cforest', following the permutation
    principle of the `mean decrease in accuracy' importance in `randomForest'.
}
\usage{
varimp(object, mincriterion = 0, conditional = FALSE, 
       threshold = 0.2, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
varimpAUC(object, mincriterion = 0, conditional = FALSE, 
       threshold = 0.2, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
}
\arguments{
  \item{object}{ an object as returned by \code{cforest}.}
  \item{mincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to include a split in the 
                       computation of the importance. The default \code{mincriterion = 0}
                       guarantees that all splits are included.}
  \item{conditional}{ a logical determining whether unconditional or conditional 
                      computation of the importance is performed. }
  \item{threshold}{ the value of the test statistic or 1 - p-value of the association 
                    between the variable of interest and a covariate that must be 
                    exceeded inorder to include the covariate in the conditioning 
                    scheme for the variable of interest (only relevant if 
                    \code{conditional = TRUE}). }
  \item{nperm}{ the number of permutations performed.}
  \item{OOB}{ a logical determining whether the importance is computed from the out-of-bag 
              sample or the learning sample (not suggested).}
  \item{pre1.0_0}{ Prior to party version 1.0-0, the actual data values
                   were permuted according to the original permutation
                   importance suggested by Breiman (2001). Now the assignments
                   to child nodes of splits in the variable of interest
                   are permuted as described by Hapfelmeier et al. (2012),
                   which allows for missing values in the explanatory
                   variables and is more efficient wrt memory consumption and 
                   computing time. This method does not apply to conditional
                   variable importances.}
}
\details{

  Function \code{varimp} can be used to compute variable importance measures
  similar to those computed by \code{\link[randomForest]{importance}}. Besides the
  standard version, a conditional version is available, that adjusts for correlations between
  predictor variables. 
  
  If \code{conditional = TRUE}, the importance of each variable is computed by permuting 
  within a grid defined by the covariates that are associated  (with 1 - p-value 
  greater than \code{threshold}) to the variable of interest.
  The resulting variable importance score is conditional in the sense of beta coefficients in   
  regression models, but represents the effect of a variable in both main effects and interactions.
  See Strobl et al. (2008) for details.

  Note, however, that all random forest results are subject to random variation. Thus, before
  interpreting the importance ranking, check whether the same ranking is achieved with a
  different random seed -- or otherwise increase the number of trees \code{ntree} in 
  \code{\link{ctree_control}}.

  Note that in the presence of missings in the predictor variables the procedure
  described in Hapfelmeier et al. (2012) is performed.

  Function \code{varimpAUC} implements AUC-based variables importances as
  described by Janitza et al. (2012).  Here, the area under the curve
  instead of the accuracy is used to calculate the importance of each variable. 
  This AUC-based variable importance measure is more robust towards class imbalance.

  For right-censored responses, \code{varimp} uses the integrated Brier score as a 
  risk measure for computing variable importances. This feature is extremely slow and
  experimental; use at your own risk.

}
\value{
  A vector of `mean decrease in accuracy' importance scores.
}
\references{ 

    Leo Breiman (2001). Random Forests. \emph{Machine Learning}, 45(1), 5--32.

    Alexander Hapfelmeier, Torsten Hothorn, Kurt Ulm, and Carolin Strobl (2012).
    A New Variable Importance Measure for Random Forests with Missing Data.
    \emph{Statistics and Computing}, \url{http://dx.doi.org/10.1007/s11222-012-9349-1}

    Torsten Hothorn, Kurt Hornik, and Achim Zeileis (2006b). Unbiased
    Recursive Partitioning: A Conditional Inference Framework.
    \emph{Journal of Computational and Graphical Statistics}, \bold{15} (3),
    651-674.  Preprint available from 
    \url{http://statmath.wu-wien.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf}

    Silke Janitza, Carolin Strobl and Anne-Laure Boulesteix (2013). An AUC-based Permutation 
    Variable Importance Measure for Random Forests. BMC Bioinformatics.2013, \bold{14} 119.
    \url{http://www.biomedcentral.com/1471-2105/14/119}

    Carolin Strobl, Anne-Laure Boulesteix, Thomas Kneib, Thomas Augustin, and Achim Zeileis (2008).
    Conditional Variable Importance for Random Forests. \emph{BMC Bioinformatics}, \bold{9}, 307. 
    \url{http://www.biomedcentral.com/1471-2105/9/307}
}
\examples{
    
   set.seed(290875)
   readingSkills.cf <- cforest(score ~ ., data = readingSkills, 
       control = cforest_unbiased(mtry = 2, ntree = 50))

   # standard importance
   varimp(readingSkills.cf)
   # the same modulo random variation
   varimp(readingSkills.cf, pre1.0_0 = TRUE)

   # conditional importance, may take a while...
   varimp(readingSkills.cf, conditional = TRUE)

   \dontrun{
   data("GBSG2", package = "TH.data")
   ### add a random covariate for sanity check
   set.seed(29)
   GBSG2$rand <- runif(nrow(GBSG2))
   object <- cforest(Surv(time, cens) ~ ., data = GBSG2, 
                     control = cforest_unbiased(ntree = 20)) 
   vi <- varimp(object)
   ### compare variable importances and absolute z-statistics
   layout(matrix(1:2))
   barplot(vi)
   barplot(abs(summary(coxph(Surv(time, cens) ~ ., data = GBSG2))$coeff[,"z"]))
   ### looks more or less the same
   }
}
\keyword{tree}
