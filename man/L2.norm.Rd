\name{L2.norm}
\alias{L2.norm}

\title{L2 norm of a given vector}
\description{
  A function for calculating the L2 norm of a given \code{numeric vector}.
}
\usage{
L2.norm(x)
}

\arguments{
  \item{x}{ A \code{numeric} vector.}
}


\value{
  The function returns the L2 norm of the given vector \code{x}.
}

\author{Jasmina Bogojeska }

\seealso{\code{\link{L1.dist}}}
\examples{
## Define a numeric vector
y <- c(1, 2, 3, 4)

## Calculate the L2 norm of the vector y
L2.norm(y)
}

\keyword{misc}
