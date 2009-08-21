\name{distribution-methods}
\docType{methods}


\alias{distribution}
\alias{distribution-methods}
\alias{distribution,RtreemixModel,missing,missing,missing-method}
\alias{distribution,RtreemixModel,character,numeric,numeric-method}


\title{Method for generating the (scaled) probablility distribution
  induced with a given mutagenetic trees mixture model}
\description{
 These functions generate the probability distribution induced with a
 given mutagenetic trees mixture model \code{model} on the space of all possible
 patterns of genetic events. The \code{model} has to be specified. The
 sampling mode and the parameters for the sampling times of the observed
 input and output probabilities are optional. The number of genetic
 events in the \code{model} cannot exceed 30.
}
\usage{
distribution(model, ...)

\method{distribution}{model}(model)

\method{distribution}{model, sampling.mode, sampling.param, output.param}(model, 
                                    sampling.mode, sampling.param, output.param)
}
\arguments{
  \item{model}{An \code{RtreemixModel} object that encodes a
    probability distribution on the set of all possible patterns.}
  \item{sampling.mode}{A \code{character} that specifies the
    sampling mode ("constant" or "exponential") for the observed input
    and output probabilities.}
  \item{sampling.param}{A \code{numeric} that specifies the
    sampling parameter for the observed input probabilities
    corresponding to the sampling mode given by \code{sampling.mode}.}
  \item{output.param}{A \code{numeric} that specifies the
    sampling parameter for the observed output probabilities
    corresponding to the sampling mode given by \code{sampling.mode}.}  
}

\value{
  The function returns a \code{dataframe} of all possible patterns with their
  corresponding probabilities derived from the specified trees mixture
  model. When the sampling mode and the sampling parameters (input and
  output) are specified their values are printed out.
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}

\author{Jasmina Bogojeska }

\seealso{
  \code{\link{RtreemixModel-class}}, \code{\link{fit-methods}}
}

\examples{
## Generate a random RtreemixModel object with 3 components.
mod <- generate(K = 3, no.events = 8, noise.tree = TRUE, prob = c(0.2, 0.8))
show(mod)

## See the probability distribution encoded by the model on the set of all possible patterns.
distr <- distribution(model = mod)
distr

## Get the probabilities.
distr$probability

## See the probability distribution encoded by the model on the set of all possible patterns
## calculated for given sampling mode, and corresponding input and output parameters.
distr1 <- distribution(model = mod, sampling.mode = "exponential", sampling.param = 1, output.param = 1)
distr1
}

\keyword{methods}
\keyword{distribution}
