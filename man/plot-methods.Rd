\name{plot-methods}
\docType{methods}

\alias{plot}
\alias{plot-methods}
\alias{plot,RtreemixModel,missing-method}



\title{Method for visualizing mutagenetic trees mixture models}
\description{
   Function for visualizing the tree components comprising a mutagenetic trees mixture model.  
   The user can also specify the \code{fontSize} used for the text labels of the nodes and the edges of the plotted trees. 
   Additionally, one can use the parameter \code{k} to plot a certain tree component from the 
   mixture model. 
}
\details{
    The value of \code{k} (that gives the tree component to be plotted) can take \code{"integer"} values from one to the number of tree components.  
}
 
\usage{
\S4method{plot}{RtreemixModel,missing}(x, y, ...)
}

\arguments{
  \item{x}{An \code{RtreemixModel} object giving the mixture model that should be visualized.}
  \item{y}{Not specified.}
  \item{...}{
    \code{fontSize} is the size of the text labels of the nodes and the edges of the tree components. The default value is 8.
    \code{k} is a \code{numeric} giving the specific tree component from the given mixture model that should be plotted. Its value 
    can be from one to the number of tree components in the given model.
  }
}
\value{
  The method returns a plot of the mixture model model.
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixData-class}}, \code{\link{RtreemixModel-class}},
  \code{\link{generate-methods}}
}

\examples{
## Generate a random RtreemixModel object.
rand.mod <- generate(K = 2, no.events = 7, noise.tree = TRUE, prob = c(0.2, 0.8))

## Visualize it.
plot(rand.mod)

## Increase the font size of the text labels in the plot.
plot(rand.mod, fontSize = 10)

## Plot the second component of the mixture model rand.mod
plot(rand.mod, k = 2)

}

\keyword{methods}
