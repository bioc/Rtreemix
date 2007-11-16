\name{RtreemixSim-class}
\docType{class}

\alias{SamplingMode,RtreemixSim}
\alias{SamplingParam,RtreemixSim}
\alias{SamplingTimes}
\alias{SimPatterns}
\alias{WaitingTimes}
\alias{getModel}
\alias{noDraws,RtreemixSim}
\alias{noDraws}

\alias{RtreemixSim-class}
\alias{SamplingMode,RtreemixSim-method}
\alias{SamplingParam,RtreemixSim-method}
\alias{SamplingTimes,RtreemixSim-method}
\alias{SimPatterns,RtreemixSim-method}
\alias{WaitingTimes,RtreemixSim-method}
\alias{getModel,RtreemixSim-method}
\alias{noDraws,RtreemixSim-method}
\alias{initialize,RtreemixSim-method}
\alias{print,RtreemixSim-method}
\alias{show,RtreemixSim-method}


\title{Class "RtreemixSim"}
\description{
  This class contains data simulated from the \code{RtreemixModel} it
  extends together with their sampling and waiting times. It also
  includes the sampling mode and the sampling parameter used for the time simulation.     
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RtreemixSim",
    Model, SimPatterns, SamplingMode, SamplingParam, WaitingTimes,
    SamplingTimes)}.
  The \code{RtreemixSim} class specifies patterns (\code{RtreemixData})
  simulated from the parent \code{RtreemixModel} together with their
  waiting and sampling times resulting from the waiting time simulation
  along the branchings in the parent model.

  The \code{Model} is an \code{RtreemixModel} object used in the data
  and time simulation process. In other words, this model is used for
  simulating patterns with their sampling and waiting times.

  The \code{SimPatterns} is an \code{RtreemixData} object that contains
  the patterns simulated from the given \code{Model}.

  The \code{SamplingMode} is a \code{character} that specifies the
  sampling mode ("constant" or "exponential") used in the time simulations.

  The \code{SamplingParam} is a \code{numeric} that specifies the
  sampling parameter corresponding to the sampling mode given by
  \code{SamplingMode}.

  The \code{WaitingTimes} is a numeric \code{vector} that specifies the
  waiting times for the simulated patterns. Its length equals the number
  of patterns in \code{SimPatterns}.

  The \code{SamplingTimes} is a numeric \code{vector} that specifies the
  sampling times for the simulated patterns. Its length equals the number
  of patterns in \code{SimPatterns}.  
}
\section{Slots}{
	 \describe{
    \item{\code{SimPatterns}:}{Object of class \code{"RtreemixData"}.}
    \item{\code{SamplingMode}:}{Object of class \code{"character"}. It
      can have one of the two possible values: "constant" or "exponential".}
    \item{\code{SamplingParam}:}{Object of class \code{"numeric"}.}
    \item{\code{WaitingTimes}:}{Object of class \code{"numeric"}. The
      length of \code{WaitingTimes} must be equal to the number
      of patterns in \code{SimPatterns}.}
    \item{\code{SamplingTimes}:}{Object of class \code{"numeric"}. The
      length of \code{SamplingTimes} must be equal to the number
      of patterns in \code{SimPatterns}.}
   }
}
\section{Extends}{
Class \code{"RtreemixModel"}, directly.
Class \code{"RtreemixData"}, by class "RtreemixModel", distance 2.
}
\section{Methods}{
  \describe{
    \item{SamplingMode}{\code{signature(object = "RtreemixSim")}: A
      method for obtaining the sampling mode ("constant" or
      "exponential") used for the time simulations.}
    \item{SamplingParam}{\code{signature(object = "RtreemixSim")}: A
      method for obtaining the sampling parameter corresponding to the
      specified \code{SamplingMode}.}
    \item{SamplingTimes}{\code{signature(object = "RtreemixSim")}: A
      method used for obtaining the sampling times of the patterns
      in \code{SimPatterns}.}
    \item{SimPatterns}{\code{signature(object = "RtreemixSim")}: A
      method used for obtaining the patterns simulated from the parent model.}
    \item{WaitingTimes}{\code{signature(object = "RtreemixSim")}: A
      method used for obtaining the waiting times of the patterns
      in \code{SimPatterns}.}
    \item{getModel}{\code{signature(object = "RtreemixSim")}: A method
      for obtaining the mixture model used in the simulations.}
    \item{noDraws}{\code{signature(object = "RtreemixSim")}: A method
      for obtaining the number of simulated patterns, i.e. the size of \code{SimPatterns}.}
	 }
}
       
\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}

\author{Jasmina Bogojeska}
\seealso{
  \code{\link{RtreemixGPS-class}}, \code{\link{RtreemixData-class}},
  \code{\link{RtreemixModel-class}},
  \code{\link{fit-methods}}, \code{\link{sim-methods}}
}
\examples{
## Generate a random RtreemixModel object with 3 components and 9 genetic events.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
show(rand.mod)

## Create an RtreemixSim object by simulating patterns with their sampling and waiting times from a given mixture model.
sim.data <- sim(model = rand.mod, sampling.mode = "exponential", sampling.param = 1, no.sim = 200)
show(sim.data)

## See the slots from the RtreemixSim object.
SimPatterns(sim.data)
SamplingMode(sim.data)
SamplingParam(sim.data)
WaitingTimes(sim.data)
SamplingTimes(sim.data)
## See model.
getModel(sim.data)
## See number of simulated patterns.
noDraws(sim.data)
}
\keyword{classes}
