\name{sim-methods}
\docType{methods}
\alias{sim}
\alias{sim-methods}
\alias{sim,RtreemixModel,missing,missing-method}
\alias{sim,RtreemixModel,character,numeric-method}

\title{Method for simulating data from a mutagenetic trees mixture model}
\description{
  This function draws a certain number of patterns from a specified
  mutagenetic trees mixture model. Thus, the mixture model has to be specified.
  When besides the mixture model also the sampling mode and its
  respective sampling parameter are specified, this function simulates
  patterns together with their waiting and sampling times from the respective model.
}
\usage{
\S4method{sim}{RtreemixModel,missing,missing}(model, sampling.mode,
sampling.param, ...)
\S4method{sim}{RtreemixModel,character,numeric}(model, sampling.mode,
                            sampling.param, ...)
}

\arguments{
  \item{model}{An object of the class \code{RtreemixModel} specifying
    the mutagenetic trees mixture model used for drawing the patterns,
    or for simulating patterns with their sampling and waiting times.}
  \item{sampling.mode}{A \code{character} that specifies the
  sampling mode ("constant" or "exponential") used in the time simulations.}
  \item{sampling.param}{A \code{numeric} that specifies the
    sampling parameter corresponding to the sampling mode given by
    \code{sampling.mode}.}
  \item{...}{
    \code{no.draws} is an \code{integer} larger than zero specifying the
    number of patterns that should be drawn from the given mixture model.
    \code{no.sim} is an \code{integer} larger than 0 giving the number of
    iterations for the waiting time simulations. Its default value is
    10.
    \code{seed} is a positive \code{integer} specifying the random generator
    seed. Its default value is (-1) and then the time is used as a
    random generator.
  }
}

\value{
  The function returns an \code{RtreemixData} object in the case when
  one wants to draw a certain number of patterns from a given mixture
  model, i.e. when only the mutagenetic trees mixture model and the
  number of patterns to be drawn are specified. When besides the model
  also the sampling mode and the sampling parameter are given, the
  function returns an object from the \code{RtreemixSim} class where the
  simulated patterns together with their sampling and waiting times are stored.
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixSim-class}}, \code{\link{RtreemixModel-class}},
  \code{\link{RtreemixData-class}}
}

\examples{
## Create a random RtreemixModel object with 3 branchings and 9 genetic events.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))

## Draw 300 samples from the randomly generated model rand.mod 
data <- sim(model = rand.mod, no.draws = 300)
show(data)

## Create an RtreemixSim object by simulating patterns with their sampling and waiting times from a given mixture model.
sim.data <- sim(model = rand.mod, sampling.mode = "exponential", sampling.param = 1, no.sim = 100)
show(sim.data)
WaitingTimes(sim.data)
SamplingTimes(sim.data)
}

\keyword{methods}
\keyword{misc}
