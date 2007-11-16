\name{likelihoods-methods}

\docType{methods}

\alias{likelihoods}
\alias{likelihoods-methods}
\alias{likelihoods,RtreemixModel,RtreemixData-method}

\title{Method for predicting the likelihoods of a set of samples with
  respect to a mutagenetic trees mixture model}
\description{
  This function predicts the (log, weighted) likelihoods of the samples in a given
  dataset according to a given mutagenetic trees mixture model. The
  dataset and the model have to be specified. 
}
\usage{
\S4method{likelihoods}{RtreemixModel,RtreemixData}(model, data)
}
\arguments{
  \item{model}{An \code{RtreemixModel} object specifying the
    probabilistic framework in which the likelihoods of the genetic
    patterns are computed.}
  \item{data}{An \code{RtreemixData} object giving the samples for
    which the likelihoods are to be calculated.}
}
\value{
  This method returns an \code{RtreemixStats} object that containes the
  weghted- and log-likelihoods of the samples in the given dataset with
  respect to the given mutagenetic trees mixture model. 
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixData-class}}, \code{\link{RtreemixModel-class}},
  \code{\link{fit-methods}}, \code{\link{distribution-methods}}
}

\examples{
## Create an RtreemixData object from a randomly generated RtreemixModel object.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
data <- sim(model = rand.mod, no.draws = 300)
show(data)

## Compute the likelihoods of the samples in data with respect to the model rand.mod
mod.stat <- likelihoods(model = rand.mod, data = data)
show(mod.stat)
}

\keyword{methods}
