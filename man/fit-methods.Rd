\name{fit-methods}
\docType{methods}

\alias{fit}
\alias{fit-methods}
\alias{fit,RtreemixData,numeric-method}



\title{Method for fitting mutagenetic trees mixture model to a given dataset}
\description{
   Function for fitting a mutagenetic trees mixture model to a given dataset
   \code{data}. The dataset and the number of trees \code{K} have to be specified.
   The function estimates K-oncogenetic trees mixture model from the
   specified data by using an EM-like learning algorithm. The first tree
   component of the model has a star topology and is referred to as the
   noise component.   
}
\details{
   When K = 1 and noise = FALSE a single mutagenetic tree is fit to the data.
   When K = 1 and noise = TRUE a star mutagenetic tree is fit to the data.
   If K > 1 the first mutagenetic tree is always the star,
   i.e. the case K > 1 and noise = FALSE is not possible.  
}
 
\usage{
\S4method{fit}{RtreemixData,numeric}(data, K, ...)
}

\arguments{
  \item{data}{An \code{RtreemixData} object giving the dataset used for
    learning the trees mixture model.}
  \item{K}{An \code{integer} larger than 0 specifying the number of
    branchings in the mixture model.}
  \item{...}{
    \code{no.start.sol} is an \code{integer} larger than 0 specifying the number of starting solutions for the k-means
    algorithm. The default value is 100.
    \code{eps} is a \code{numeric} giving the minimum conditional probability to include edge. The
    default value is 0.01.
    \code{weighing} is a \code{logical} specifying whether to use special
    weights log(Pr(v)) for the edges (root, v). The default value is \code{FALSE}.
    \code{equal.edgeweights} is a \code{logical} specifying whether to use
    equal edge weights in the noise component. The default value is
    \code{TRUE}.
    \code{seed} is a positive \code{integer} specifying the random generator
    seed. The default value is (-1) and then the time is used as a
    random generator.
    \code{noise} is a \code{logical} indicating the presence of a noise
    (star) component in the fitted mixture model. It is mostly relevant
    for models with a single tree component, since it is assumed that mixture models with
    at least two components always have the noise as a first component.
  }
}
\value{
  The method returns an \code{RtreemixModel} object that represents the
  K-trees mixture model learned from the given dataset.
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixData-class}}, \code{\link{RtreemixModel-class}},
  \code{\link{generate-methods}}, \code{\link{bootstrap-methods}},
  \code{\link{confIntGPS-methods}}
}

\examples{
## Create an RtreemixData object from a randomly generated RtreemixModel object.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
data <- sim(model = rand.mod, no.draws = 300)
show(data)

## Create an RtreemixModel object by fitting model to the given data.
mod <- fit(data = data, K = 3, equal.edgeweights = TRUE, noise = TRUE)
show(mod)
## See the number of tree components in the mixture model.
numTrees(mod)
## See the weights of the branchings from the fitted mixture model.
Weights(mod)
## See a specific tree component k.
getTree(object = mod, k = 2)
}

\keyword{methods}

