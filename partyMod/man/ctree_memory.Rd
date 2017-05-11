\name{Memory Allocation}
\alias{ctree_memory}
\title{ Memory Allocation }
\description{

    This function sets up the memory needed for tree growing. It might be 
    convenient to allocate memory only once but build multiple trees.
}
\usage{
ctree_memory(object, MPinv = FALSE)
}
\arguments{
  \item{object}{an object of class \code{LearningSample}.}
  \item{MPinv}{a logical indicating whether memory for the Moore-Penrose
               inverse of covariance matrices should be allocated. }
}
\details{

  This function is normally not to be called by users. However, for
  performance reasons it might be nice to allocate memory and re-fit trees
  using the same memory for the computations. Below is an example.

}
\value{
  An object of class \code{TreeFitMemory}.
}
\examples{

    set.seed(290875)

    ### setup learning sample
    airq <- subset(airquality, !is.na(Ozone))
    ls <- dpp(conditionalTree, Ozone ~ ., data = airq)

    ### setup memory and controls 
    mem <- ctree_memory(ls)
    ct <- ctree_control(teststat = "max")

    ### fit 50 trees on bootstrap samples
    bs <- rmultinom(50, nrow(airq), rep(1, nrow(airq))/nrow(airq))
    storage.mode(bs) <- "double"
    cfit <- conditionalTree@fit
    ens <- apply(bs, 2, function(w) cfit(ls, ct, weights = w, 
                                         fitmem = mem))

}
\keyword{misc}
