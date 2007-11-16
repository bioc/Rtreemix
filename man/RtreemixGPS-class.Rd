\name{RtreemixGPS-class}

\docType{class}

\alias{GPS}
\alias{Model}
\alias{SamplingMode}
\alias{SamplingParam}
\alias{gpsCI}
\alias{getData}



\alias{RtreemixGPS-class}
\alias{GPS,RtreemixGPS-method}
\alias{Model,RtreemixGPS-method}
\alias{SamplingMode,RtreemixGPS-method}
\alias{SamplingParam,RtreemixGPS-method}
\alias{getData,RtreemixGPS-method}
\alias{gpsCI,RtreemixGPS-method}
\alias{initialize,RtreemixGPS-method}
\alias{print,RtreemixGPS-method}
\alias{show,RtreemixGPS-method}



\title{Class "RtreemixGPS"}
\description{
  A class for describing the genetic progression scores (GPS) of a given
  set of patterns resulting from a waiting time simulation along the
  edges of the tree components of a given mutagenetic trees mixture
  model. It also contains GPS confidence intervals derived with the
  bootstrap method. 
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RtreemixGPS",
    Data, Model, SamplingMode, SamplingParam, GPS, gpsCI)}.
  The \code{RtreemixGPS} class contains the GPS values each assigned to 
  the corresponding pattern from the dataset given by \code{Data} (the
  parent class). The GPS values are derived in a waiting time simulation
  for a specified sampling mode and its corresponding sampling
  parameter. Moreover, this class specifies the confidence intervals for
  the GPS values derived with the bootstrap method.

  The \code{Data} is an \code{RtreemixData} object that specifies the
  patterns for which the GPS values are calculated.

  The \code{Model} is an \code{RtreemixModel} object that specifies the
  mutagenetic trees mixture model used for deriving the GPS values.

  The \code{SamplingMode} is a \code{character} that specifies the
  sampling mode ("constant" or "exponential") used in the waiting time
  simulations.

  The \code{SamplingParam} is a \code{numeric} that specifies the
  sampling parameter corresponding to the sampling mode given by
  \code{SamplingMode}.

  The \code{GPS} is a numeric \code{vector} that specifies the
  GPS value of each pattern in the given dataset \code{Data}.
  Its length equals the number of patterns in \code{Data}.

  The \code{gpsCI} is a numeric \code{matrix} that specifies the
  confidence intervals for each GPS value in the vector \code{GPS}. The
  number of rows equals the number of patients in \code{Data} and the
  number of columns equals 2. The first column gives the lower bound and
  the second column gives the upper bound of each confidence interval.  
}
\section{Slots}{
	 \describe{
    \item{\code{Model}:}{Object of class \code{"RtreemixModel"}.}
    \item{\code{SamplingMode}:}{Object of class \code{"character"}. It
      can have one of the two possible values: "constant" or "exponential".}
    \item{\code{SamplingParam}:}{Object of class \code{"numeric"}.}
    \item{\code{GPS}:}{Object of class \code{"numeric"}. The
      length of \code{GPS} must be equal to the number
      of patterns in the parent \code{RtreemixData} object.}
    \item{\code{gpsCI}:}{Object of class \code{"matrix"}. It number of
      columns has to be 2 and the number of rows has to be equal to the
      length of \code{GPS}.}
  }
}
\section{Extends}{
Class \code{"RtreemixData"}, directly.
}
\section{Methods}{
  \describe{
    \item{GPS}{\code{signature(object = "RtreemixGPS")}: A method for
      obtaining the GPS values corresponding to the patterns in the
      parent \code{RtreemixData} object.}
    \item{Model}{\code{signature(object = "RtreemixGPS")}: A method for
      obtaining the model used for deriving the GPS values.}
    \item{SamplingMode}{\code{signature(object = "RtreemixGPS")}: A
      method for obtaining the sampling mode ("constant" or
      "exponential") used for the waiting time simulations.}
    \item{SamplingParam}{\code{signature(object = "RtreemixGPS")}: A
      method for obtaining the sampling parameter corresponding to the
      specified \code{SamplingMode}.}
    \item{getData}{\code{signature(object = "RtreemixGPS")}: A method
      for obtaining the set of patterns for which the GPS values are calculated. }
    \item{gpsCI}{\code{signature(object = "RtreemixGPS")}: A method for
      obtaining the GPS confidence intervals.}
	 }
}
\references{Estimating cancer survival and clinical outcome based on
  genetic tumor progression scores, J. Rahnenf\"urer et al. }
\author{Jasmina Bogojeska}

\note{
    The GPS examples are time consuming. They are commented out because of the time restrictions of the check of the package.
    For trying out the code please copy it and uncomment it.
} 

\seealso{
  \code{\link{RtreemixData-class}}, 
  \code{\link{RtreemixModel-class}},
  \code{\link{gps-methods}},
  \code{\link{fit-methods}}, \code{\link{confIntGPS-methods}}
}
\examples{
## Generate a random RtreemixModel object with 3 components and 9 genetic events.
#mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
#show(mod)
## Generate an artificial dataset from the model mod.
#data <- sim(model = mod, no.draws = 300)
#show(data)

## Create an RtreemixGPS object by calculating the GPS for all possible patterns.
#modGPS.all <- gps(model = mod, no.sim = 1000)
#show(modGPS.all)
## Create an RtreemixGPS object by calculating the GPS for the data based on the model mod.
#modGPS <- gps(model = mod, data = data, no.sim = 1000)
#show(modGPS)

## See the slots from the RtreemixGPS object. 
#Model(modGPS)
#SamplingMode(modGPS)
#SamplingParam(modGPS)
#GPS(modGPS)
## See data.
#getData(modGPS)

## Create an RtreemixGPS object by calculating GPS values for a given dataset
## and their 95\% confidence intervals using the bootstrap method.
#modGPS2 <- confIntGPS(data = data, K = 2, B = 10)
#show(modGPS2)

## See the GPS values for the object modGPS2 and their confidence intervals.
#GPS(modGPS2)
#gpsCI(modGPS2)
}
\keyword{classes}
