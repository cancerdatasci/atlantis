\name{Panel Generating Functions}
\encoding{latin1}

\alias{node_inner}
\alias{node_terminal}
\alias{edge_simple}
\alias{node_surv}
\alias{node_barplot}
\alias{node_boxplot}
\alias{node_hist}
\alias{node_density}
\alias{node_scatterplot}
\alias{node_bivplot}

\title{ Panel-Generators for Visualization of Party Trees }
\description{
    The plot method for \code{BinaryTree} and \code{mob} objects are rather
    flexible and can be extended by panel functions. Some pre-defined
    panel-generating functions of class \code{grapcon_generator}
    for the most important cases are documented here.
}
\usage{
node_inner(ctreeobj, digits = 3, abbreviate = FALSE, 
  fill = "white", pval = TRUE, id = TRUE)
node_terminal(ctreeobj, digits = 3, abbreviate = FALSE, 
  fill = c("lightgray", "white"), id = TRUE)
edge_simple(treeobj, digits = 3, abbreviate = FALSE)
node_surv(ctreeobj, ylines = 2, id = TRUE, \dots)
node_barplot(ctreeobj, col = "black", fill = NULL, beside = NULL,
  ymax = NULL, ylines = NULL, widths = 1, gap = NULL, 
  reverse = NULL, id = TRUE)
node_boxplot(ctreeobj, col = "black", fill = "lightgray", 
  width = 0.5,  yscale = NULL, ylines = 3, cex = 0.5, id = TRUE)
node_hist(ctreeobj, col = "black", fill = "lightgray", 
  freq = FALSE, horizontal = TRUE, xscale = NULL, ymax = NULL, 
  ylines = 3, id = TRUE, \dots)
node_density(ctreeobj, col = "black", rug = TRUE, 
  horizontal = TRUE, xscale = NULL, yscale = NULL, ylines = 3, 
  id = TRUE)
node_scatterplot(mobobj, which = NULL, col = "black", 
  linecol = "red", cex = 0.5, pch = NULL, jitter = FALSE, 
  xscale = NULL, yscale = NULL, ylines = 1.5, id = TRUE, 
  labels = FALSE)
node_bivplot(mobobj, which = NULL, id = TRUE, pop = TRUE,
  pointcol = "black", pointcex = 0.5,
  boxcol = "black", boxwidth = 0.5, boxfill = "lightgray",
  fitmean = TRUE, linecol = "red",
  cdplot = FALSE, fivenum = TRUE, breaks = NULL,
  ylines = NULL, xlab = FALSE, ylab = FALSE, margins = rep(1.5, 4), \dots)
}
\arguments{
  \item{ctreeobj}{ an object of class \code{BinaryTree}.}
  \item{treeobj}{ an object of class \code{BinaryTree} or \code{mob}.}
  \item{mobobj}{ an object of class \code{mob}.}
  \item{digits}{ integer, used for formating numbers. }
  \item{abbreviate}{ logical indicating whether strings should be 
    abbreviated. }
  \item{col, pointcol}{ a color for points and lines. }
  \item{fill}{ a color to filling rectangles. }
  \item{pval}{ logical. Should p values be plotted?}
  \item{id}{ logical. Should node IDs be plotted?}
  \item{ylines}{ number of lines for spaces in y-direction. }
  \item{widths}{ widths in barplots. }
  \item{width, boxwidth}{ width in boxplots. }
  \item{gap}{ gap between bars in a barplot (\code{node_barplot}). }
  \item{yscale}{ limits in y-direction}
  \item{xscale}{ limits in x-direction}
  \item{ymax}{ upper limit in y-direction}
  \item{beside}{ logical indicating if barplots should be side by side or stacked. }
  \item{reverse}{logical indicating whether the order of levels should be reversed
    for barplots.}
  \item{horizontal}{ logical indicating if the plots should be horizontal. }
  \item{freq}{logical; if \code{TRUE}, the histogram graphic is a representation
    of frequencies. If \code{FALSE}, probabilities are plotted.}
  \item{rug}{logical indicating if a rug representation should be added. }
  \item{which}{ numeric or character vector indicating which of the regressor
    variables should be plotted (default = all).}
  \item{linecol}{ color for fitted model lines.}
  \item{cex, pointcex}{character extension of points in scatter plots.}
  \item{pch}{plotting character of points in scatter plots.}
  \item{jitter}{logical. Should the points be jittered in y-direction?}
  \item{labels}{logical. Should axis labels be plotted?}
  \item{pop}{logical. Should the panel viewports be popped?}
  \item{boxcol}{color for box plot borders.}
  \item{boxfill}{fill color for box plots.}
  \item{fitmean}{logical. Should lines for the predicted means from the model
    be added?}
  \item{cdplot}{logical. Should CD plots (or spinograms) be used for visualizing
    the dependence of a categorical on a numeric variable?}
  \item{fivenum}{logical. When using spinograms, should the five point summary
    of the explanatory variable be used for determining the breaks?}
  \item{breaks}{a (list of) numeric vector(s) of breaks for the spinograms. If set to \code{NULL}
    (the default), the \code{breaks} are chosen according to the \code{fivenum}
    argument.}
  \item{xlab, ylab}{ character with x- and y-axis labels. Can also be logical: if \code{FALSE}
    axis labels are surpressed, if \code{TRUE} they are taken from the underlying data.
    Can be a vector of labels for \code{xlab}. }
  \item{margins}{margins of the viewports.}
  \item{\dots}{ additional arguments passed to callies.}
}

\details{
  The \code{plot} methods for \code{BinaryTree} and \code{mob} objects provide an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. The panel functions to be used
  should depend only on the node being visualzied, however, for setting up
  an appropriate panel function, information from the whole tree is typically
  required. Hence, \pkg{party} adopts the framework of \code{grapcon_generator}
  (graphical appearance control) from the \pkg{vcd} package (Meyer, Zeileis and
  Hornik, 2005) and provides several panel-generating functions. For convenience,
  the panel-generating functions \code{node_inner} and \code{edge_simple} 
  return panel functions to draw inner nodes and left and right edges. 
  For drawing terminal nodes, the functions returned by the other panel 
  functions can be used. The panel generating function \code{node_terminal} 
  is a terse text-based representation of terminal nodes.

  Graphical representations of terminal nodes are available and depend on
  the kind of model and the measurement scale of the variables modelled.

  For univariate regressions (typically fitted by \code{ctree}),
  \code{node_surv} returns a functions that plots Kaplan-Meier curves in each 
  terminal node; \code{node_barplot}, \code{node_boxplot}, \code{node_hist} and
  \code{node_density} can be used to plot bar plots, box plots, histograms and 
  estimated densities into the terminal nodes.

  For multivariate regressions (typically fitted by \code{mob}),
  \code{node_bivplot} returns a panel function that creates bivariate plots
  of the response against all regressors in the model. Depending on the scale
  of the variables involved, scatter plots, box plots, spinograms (or CD plots)
  and spine plots are created. For the latter two \code{\link[vcd]{spine}} and
  \code{\link[vcd]{cd_plot}} from the \pkg{vcd} package are re-used.
}

\references{
  David Meyer, Achim Zeileis, and Kurt Hornik (2006).
  The Strucplot Framework: Visualizing Multi-Way Contingency Tables with vcd.
  \emph{Journal of Statistical Software}, \bold{17}(3).
  \url{http://www.jstatsoft.org/v17/i03/}
}

\examples{

  set.seed(290875)

  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq)

  ## default: boxplots
  plot(airct)
  
  ## change colors
  plot(airct, tp_args = list(col = "blue", fill = hsv(2/3, 0.5, 1)))
  ## equivalent to
  plot(airct, terminal_panel = node_boxplot(airct, col = "blue", 
                                            fill = hsv(2/3, 0.5, 1)))

  ### very simple; the mean is given in each terminal node
  plot(airct, type = "simple")

  ### density estimates
  plot(airct, terminal_panel = node_density)
    
  ### histograms 
  plot(airct, terminal_panel = node_hist(airct, ymax = 0.06, 
                                         xscale = c(0, 250)))
}
\keyword{hplot}
