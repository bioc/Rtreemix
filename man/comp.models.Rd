\name{Models}
\alias{comp.models}
\alias{comp.models.levels}

\title{Functions for comparing the tree topologies of two mutagenetic trees mixture models}
\description{
  These functions implement a similarity measure
  for comparing the topologies of the trees of two mixture models
  \code{mixture1} and \code{mixture2}. \code{comp.models} chaaracterizes
  the similarity of the models based on sum of the number of different
  edges of matched tree components (similarity
  pairs). \code{comp.models.levels} quantifies the similarity of two
  mixture models by adding to the edge ddifference of each similarity
  pair in the previously described sum the L1 distance of the level vectors of the
  trees comprising the pair. A level vector can be associated to each
  tree component and denotes the depth of each of the genetic
  events in the tree. 
  It is necessary that the two models have the same number of tree
  components build on the same number of genetic events. It is assumed
  that the mixtures have at least two tree components.   
}
\usage{
comp.models(mixture1, mixture2)
comp.models.levels(mixture1, mixture2)
}

\arguments{
  \item{mixture1}{An \code{RtreemixModel} object specifying the first
    component for the similarity calculation.}
  \item{mixture2}{An \code{RtreemixModel} object specifying the second
    component for the similarity calculation. The number of tree
    components equals the one of \code{mixture1}.}
}

\details{
  The value returned by the function \code{comp.models} is between 0 (no
  similarity) and 1 (identical models).
}

\value{
  The functions return a numeric value that quantifies the similarity
  of the tree topologies of two mixture models.
}

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixModel-class}}, \code{\link{comp.trees}},
  \code{\link{fit-methods}}, \code{\link{stability.sim}}
}


\examples{
## Generate two random RtreemixModel objects each with 3 components.
rand.mod1 <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob =
c(0.2, 0.8))
rand.mod2 <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob =
c(0.2, 0.8))

## Compare the topologies of the tree components of the two randomly
## generated models
comp.models(rand.mod1, rand.mod2)
comp.models.levels(rand.mod1, rand.mod2)
}

\keyword{misc}
