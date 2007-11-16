\name{distances}

\alias{L1.dist}
\alias{euclidian.dist}
\alias{cosin.dist}
\alias{rank.cor.dist}

\title{Different distances between two given vectors}
\description{
 These functions are used for calculating different distances between
 two given vectors. Thus, \code{L1.dist} calculates the L1 distance,
 \code{cosin.dist} calculates the cosine distance, \code{euclidian.dist}
 computes the Euclidian distance, and \code{rank.cor.dist} computes
 the rank correlation distance. The vectors have to have same length.
 When using \code{rank.cor.dist} the vectors have to have length larger
 than 4.
}
\usage{
L1.dist(p, q)
cosin.dist(p, q)
euclidian.dist(x, y)
rank.cor.dist(x, y)
}

\arguments{
  \item{p}{ A \code{numeric} vector specifying the first component for
    the distance calculation. It has to have the same length as \code{q}.}
  \item{q}{ A \code{numeric} vector specifying the second component for
    the distance calculation. }
  \item{x}{ Same as \code{p}. }
  \item{y}{ Same as \code{q}.}

}

\value{
  The functions return the distance between the two given vectors. 
}

\author{Jasmina Bogojeska}

\seealso{\code{\link{kullback.leibler}}, \code{\link{L2.norm}}, \code{\link{stability.sim}}}
\examples{
## Define two numeric vectors with equal lengths (> 4).
x <- c(1, 2, 3, 4, 5)
y <- c(5, 6, 7, 8, 9)

## Calculate the L1 distance between the vectors x and y
L1.dist(x, y)

## Calculate the Euclidian distance between the vectors x and y
euclidian.dist(x, y)

## Calculate the cosine distance between the vectors x and y
cosin.dist(x, y)

## Calculate the rank-correlation distance between the vectors x and y
rank.cor.dist(x, y)

}

\keyword{misc}
