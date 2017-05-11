\name{plot.mob}
\alias{plot.mob}

\title{ Visualization of MOB Trees }

\description{
  \code{plot} method for \code{mob} objects with
  extended facilities for plugging in panel functions.
}

\usage{
\method{plot}{mob}(x, terminal_panel = node_bivplot, tnex = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{mob}.}
  \item{terminal_panel}{a panel function or panel-generating function of
    class \code{"grapcon_generator"}. See \code{\link{plot.BinaryTree}} for
    more details.}
  \item{tnex}{a numeric value giving the terminal node extension in relation
    to the inner nodes.}
  \item{\dots}{ further arguments passed to \code{\link{plot.BinaryTree}}.}
}

\details{
  This \code{plot} method for \code{mob} objects simply calls the
  \code{\link{plot.BinaryTree}} method, setting a different \code{terminal_panel}
  function by default (\code{\link{node_bivplot}}) and \code{tnex} value.
}


\seealso{\code{\link{node_bivplot}}, \code{\link{node_scatterplot}},
  \code{\link{plot.BinaryTree}}, \code{\link{mob}}}

\examples{

set.seed(290875)

if(require("mlbench")) {

## recursive partitioning of a linear regression model
## load data
data("BostonHousing", package = "mlbench")
## and transform variables appropriately (for a linear regression)
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2
## as well as partitioning variables (for fluctuation testing)
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, 
                             labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)

## partition the linear regression model medv ~ lstat + rm
## with respect to all remaining variables:
fm <- mob(medv ~ lstat + rm | zn + indus + chas + nox + age + dis + 
                              rad + tax + crim + b + ptratio,
  control = mob_control(minsplit = 40), data = BostonHousing, 
  model = linearModel)

## visualize medv ~ lstat and medv ~ rm
plot(fm)

## visualize only one of the two regressors
plot(fm, tp_args = list(which = "lstat"), tnex = 2)
plot(fm, tp_args = list(which = 2), tnex = 2)

## omit fitted mean lines
plot(fm, tp_args = list(fitmean = FALSE))

## mixed numerical and categorical regressors 
fm2 <- mob(medv ~ lstat + rm + chas | zn + indus + nox + age + 
                                      dis + rad,
  control = mob_control(minsplit = 100), data = BostonHousing, 
  model = linearModel)
plot(fm2)

## recursive partitioning of a logistic regression model
data("PimaIndiansDiabetes", package = "mlbench")
fmPID <- mob(diabetes ~ glucose | pregnant + pressure + triceps + 
                                  insulin + mass + pedigree + age,
  data = PimaIndiansDiabetes, model = glinearModel, 
  family = binomial())
## default plot: spinograms with breaks from five point summary
plot(fmPID)
## use the breaks from hist() instead
plot(fmPID, tp_args = list(fivenum = FALSE))
## user-defined breaks
plot(fmPID, tp_args = list(breaks = 0:4 * 50))
## CD plots instead of spinograms
plot(fmPID, tp_args = list(cdplot = TRUE))
## different smoothing bandwidth
plot(fmPID, tp_args = list(cdplot = TRUE, bw = 15))

}
}
\keyword{hplot}
