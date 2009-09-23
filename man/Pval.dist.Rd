\name{Pval.dist}
\alias{Pval.dist}

\title{p-value of a given similarity value}
\description{
  This function calculates the p-value of a given similarity value,
  i.e. the probability for obtaining the same or a smaller value than
  the given one in a vector of random similarity values. The p-value is
  used to determine whether the given similarity value is significant.
}
\usage{
Pval.dist(dist.val, random.vals)
}

\arguments{
  \item{dist.val}{ A \code{numeric} value quanifying the similarity
    for which a p-value should be calculated. }
  \item{random.vals}{A \code{numeric} vector of random similarities used
  for calculating the p-value.}
}


\value{
 It returns a \code{numeric} value between 0 and 1 that specifies the
 p-value of the given \code{dist.val}.
}

\author{ Jasmina Bogojeska }

\seealso{
 \code{\link{L1.dist}}, \code{\link{kullback.leibler}},
    \code{\link{comp.models}}, \code{\link{stability.sim}}}

\examples{
## The function is currently defined as
function(dist.val, random.vals) {
  return((sum(random.vals <= dist.val) + 1) /(length(random.vals) + 1))
  }

## Define the similarity value and a vector of random similarities
sim.val <- 0.2
rand.vals <- c(0.1, 0.24, 0.28, 0.35, 0.15, 0.5, 0.14, 0.6, 0.8, 0.3)

## Calculate the p-value of sim.val using the vector of random
## similarities
Pval.dist(dist.val = sim.val, random.vals = rand.vals)
 
}

\keyword{misc}
