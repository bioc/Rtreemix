\name{comp.trees}

\alias{comp.trees}
\alias{comp.trees.levels}

\title{Functions for quantifying the diversity of the nontrivial trees in
a mutagenetic trees mixture model}
\description{
 These functions implement a similarity measure for comparing the topologies
 of the nontrivial tree components of a specified mixture
 model, and thereby quantifying their diversity.
 All possible pairs of nontrivial components are considered when
 computing the similarity. \code{comp.trees} uses the sum of the number of
 different edges of all pairs for caracterizing
 the difference of the trees in the model. \code{comp.trees.levels} uses
 the sum of the number of different edges of all pairs and the
 corresponding L1 distances of their level vectors. The model must have at
 least two nontrivial components.
}
\usage{
comp.trees(model)
comp.trees.levels(model)
}

\arguments{
  \item{model}{An \code{RtreemixModel} object.}
}

\value{
  The functions return a numeric value that quantifies the similarity
  (or diversity) of the nontrivial tree topologies of a given mixture models.
}

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixModel-class}}, \code{\link{comp.models}},
  \code{\link{fit-methods}}, \code{\link{stability.sim}}
}

\examples{
## Generate two random RtreemixModel objects each with 3 components.
mix1 <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
mix2 <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
## Inspect the diversity of the nontrivial tree components in a given model
## using the number of distinct edges and the levels of the events in
## the treesas dissimilarity measure.
comp.trees.levels(model = mix1)
comp.trees.levels(model = mix2)
}


\keyword{misc}

