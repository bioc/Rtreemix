\name{kullback.leibler}
\alias{kullback.leibler}

\title{Kullback-Leibler divergence}
\description{
  A function for calculating the Kullback-Leibler divergence between
  two discrete probability distributions. The vectors specifying the
  probability distributions must have the same length. 
}
\usage{
kullback.leibler(p, q)
}

\arguments{
  \item{p}{A \code{numeric} vector specifying the the first probability
    distribution. It has to have the same length as \code{q}.}
  \item{q}{A \code{numeric} vector specifying the the second probability
    distribution.}
}

\value{
 The function returns the Kullback-Leibler divergence between the two
 specified descrete probability distributions.
}

\author{Jasmina Bogojeska}
\section{Warning }{
  The function does not check whether the values in the vectors
  specifying the discrete probability distributions sum up to one.
} 
\seealso{
  \code{\link{L1.dist}}, \code{\link{L2.norm}},  \code{\link{stability.sim}} 
}

\examples{
## Define two discrete probability distributions with equal lengths.
p <- c(0.1, 0.2, 0.3, 0.4)
q <- c(0.2, 0.5, 0.1, 0.2)

## Calculate the Kullback-Leibler divergence
## between the probability distributions p and q
kullback.leibler(p, q)
}

\keyword{misc}

