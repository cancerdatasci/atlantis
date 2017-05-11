\name{Plot BinaryTree}
\encoding{latin1}

\alias{plot.BinaryTree}

\title{ Visualization of Binary Regression Trees }

\description{
  \code{plot} method for \code{BinaryTree} objects with
  extended facilities for plugging in panel functions.
}

\usage{
\method{plot}{BinaryTree}(x, main = NULL, type = c("extended", "simple"),
     terminal_panel = NULL, tp_args = list(),
     inner_panel = node_inner, ip_args = list(),
     edge_panel = edge_simple, ep_args = list(),
     drop_terminal = (type[1] == "extended"), 
     tnex = (type[1] == "extended") + 1, newpage = TRUE,  
     pop = TRUE, \dots)
}

\arguments{
  \item{x}{ an object of class \code{BinaryTree}.}
  \item{main}{ an optional title for the plot.}
  \item{type}{ a character specifying the complexity of the plot:
               \code{extended} tries to visualize the distribution of the
               response variable in each terminal node whereas \code{simple} 
               only gives some summary information.}
  \item{terminal_panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the terminal nodes.
                         Alternatively, a panel generating function of class
			 \code{"grapcon_generator"} that is called with arguments
			 \code{x} and \code{tp_args} to set up a panel function.
			 By default, an appropriate panel function is chosen 
			 depending on the scale of the dependent variable.}
  \item{tp_args}{ a list of arguments passed to \code{terminal_panel} if this
                  is a \code{"grapcon_generator"} object.}
  \item{inner_panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the inner nodes.
                         Alternatively, a panel generating function of class
			 \code{"grapcon_generator"} that is called with arguments
			 \code{x} and \code{ip_args} to set up a panel function.}
  \item{ip_args}{ a list of arguments passed to \code{inner_panel} if this
                  is a \code{"grapcon_generator"} object.}
  \item{edge_panel}{ an optional panel function of the form 
                    \code{function(split, ordered = FALSE, left = TRUE)}
                    plotting the edges.
                    Alternatively, a panel generating function of class
		    \code{"grapcon_generator"} that is called with arguments
		    \code{x} and \code{ip_args} to set up a panel function.}		       
  \item{ep_args}{ a list of arguments passed to \code{edge_panel} if this
                  is a \code{"grapcon_generator"} object.}
  \item{drop_terminal}{ a logical indicating whether all terminal nodes
                        should be plotted at the bottom.}
  \item{tnex}{a numeric value giving the terminal node extension in relation
    to the inner nodes.}
  \item{newpage}{ a logical indicating whether \code{grid.newpage()} should be called. }
  \item{pop}{ a logical whether the viewport tree should be popped before
              return. }
  \item{\dots}{ additional arguments passed to callies.}
}

\details{
  This \code{plot} method for \code{BinaryTree} objects provides an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. Panel functions for plotting
  inner nodes, edges and terminal nodes are available for the most important
  cases and can serve as the basis for user-supplied extensions, see
  \code{\link{node_inner}} and \code{vignette("party")}.

  More details on the ideas and concepts of panel-generating functions and
  \code{"grapcon_generator"} objects in general can be found in Meyer, Zeileis
  and Hornik (2005).
}

\references{
  David Meyer, Achim Zeileis, and Kurt Hornik (2006).
  The Strucplot Framework: Visualizing Multi-Way Contingency Tables with vcd.
  \emph{Journal of Statistical Software}, \bold{17}(3).
  \url{http://www.jstatsoft.org/v17/i03/}
}


\seealso{\code{\link{node_inner}}, \code{\link{node_terminal}}, \code{\link{edge_simple}}, 
  \code{\link{node_surv}}, \code{\link{node_barplot}}, \code{\link{node_boxplot}}, 
  \code{\link{node_hist}}, \code{\link{node_density}}}

\examples{

  set.seed(290875)

  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq)

  ### regression: boxplots in each node
  plot(airct, terminal_panel = node_boxplot, drop_terminal = TRUE)

  if(require("TH.data")) {
  ## classification: barplots in each node
  data("GlaucomaM", package = "TH.data")
  glauct <- ctree(Class ~ ., data = GlaucomaM)
  plot(glauct)
  plot(glauct, inner_panel = node_barplot,
    edge_panel = function(ctreeobj, ...) { function(...) invisible() },
    tnex = 1)

  ## survival: Kaplan-Meier curves in each node
  data("GBSG2", package = "TH.data")
  library("survival")
  gbsg2ct <- ctree(Surv(time, cens) ~ ., data = GBSG2)
  plot(gbsg2ct)
  plot(gbsg2ct, type = "simple")  
  }

}
\keyword{hplot}
