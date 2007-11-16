\name{generate-methods}
\docType{methods}

\alias{generate}
\alias{generate-methods}
\alias{generate,numeric,numeric-method}

\title{Method for generating a random mutagenetic trees mixture model}
\description{
  Function for generating a random mutagenetic mixture model. Each tree component
  from the model is drawn uniformly at random from the tree topology
  space by using the Pr\"ufer encoding of trees.  
  The number of tree components and the number of genetic events
  have to be specified.
}
\usage{
\S4method{generate}{numeric,numeric}(K, no.events, ...)
}
\arguments{
  \item{K}{An \code{integer} larger than 0 specifying the number of
    branchings in the mixture model.}
  \item{no.events}{An \code{integer} larger than 0 specifying the number of
    genetic events in the mixture model.}
  \item{...}{
    \code{noise.tree} is a \code{logical} indicating the presence of a noise
    (star) component in the random mixture model. The default value is
    \code{TRUE}.
    \code{equal.edgeweights} is a \code{logical} specifying whether to use
    equal edge weights in the noise component. The default value is
    \code{TRUE}.
    \code{prob} is a \code{numeric} vector of length 2 specifying the
    boundaries for the edge weights of the randomly generated trees. The
    first component of the vector (the lower boundary) must be smaller
    than the second component (the upper boundary). The default value
    is (0.0, 1.0).
    \code{seed} is a positive \code{integer} specifying the random generator
    seed. The default value is (-1) and then the time is used as a
    random generator.
  }
}
\value{
  The method returns an \code{RtreemixModel} object that represents the
  randomly generated K-trees mixture model.
}
\references{Beweis eines Satzes \"uber Permutationen, H. Pr\"ufer; Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.; Model Selection for Mixtures of
  Mutagenetic Trees, Yin et al. }

\author{Jasmina Bogojeska}

\seealso{
  \code{\link{RtreemixModel-class}}
}

\examples{
## Generate a random RtreemixModel object with 3 components and 9 genetic events.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
show(rand.mod)
}

\keyword{methods}
\keyword{models}
