\name{get.tree.levels}
\alias{get.tree.levels}

\title{Function for constructing level vectors}
\description{
 Function that assignes to each node the level at which that node is
 in a specific tree (\code{tree.num}) of the mutagenetic trees mixture model \code{mixture}.
 The \code{start.val} is the number assigned to the events pruned from
 the tree. This usually is the maximum depth of the tree with which the
 tree specified with \code{tree.num} will be compared.
}
\usage{
get.tree.levels(mixture, tree.num, start.val)
}

\arguments{
  \item{mixture}{An object of the class \code{RtreemixModel}.}
  \item{tree.num}{A \code{numeric} specifying the tree component from
    \code{mixture} used for creating the level vector.}
  \item{start.val}{A \code{numeric} specifying the number assigned to
    the pruned events.}
}

\value{
 The function returns a named \code{numeric} vector. Its length equals the
 number of genetic events in \code{mixture} minus one (for the initial
 null event which is always on level 0). The vector names correspond to
 the names of the genetic events and each vector component gives the
 level at which the respective event is in the \code{num.tree} tree of \code{mixture}.
}


\author{Jasmina Bogojeska}


\seealso{
  \code{\link{comp.models}}, \code{\link{comp.trees}}, \code{\link{stability.sim}}, 
  \code{\link{RtreemixModel-class}}, \code{\link{fit-methods}}
}




\examples{
## Generate two random RtreemixModel objects each with 3 components.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob =
c(0.2, 0.8))

## Get the tree levels of the 2nd component of the model rand.mod.
get.tree.levels(mixture = rand.mod, tree.num = 2, start.val = 10)
}

\keyword{misc}
