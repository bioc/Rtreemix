\name{gps-methods}
\docType{methods}

\alias{gps}
\alias{gps-methods}
\alias{gps,RtreemixModel,RtreemixData-method}
\alias{gps,RtreemixModel,matrix-method}
\alias{gps,RtreemixModel,missing-method}




\title{Methods for predicting the GPS of given dataset by
  using a given mutagenetic trees mixture model}
\description{
  These functions compute the genetic progression score (GPS) of each
  sample in the given \code{data} by performing a waiting time
  simulation along the branchings of the mixture model \code{model}. The
  model has to be specified. If a dataset is missing a GPS for all
  possible patterns is calculated. The number of events of the samples
  in \code{data} equals the number of genetic events in the \code{model}.    
}

\usage{
\S4method{gps}{RtreemixModel,RtreemixData}(model, data, ...)
\S4method{gps}{RtreemixModel,matrix}(model, data, ...)
\S4method{gps}{RtreemixModel,missing}(model, data, ...)
}

\section{Methods}{\describe{
    \item{model = "RtreemixModel", data = "RtreemixData", ...}{A method for calculating
      the GPS values of the data given as \code{RtreemixData} object.}
    \item{model = "RtreemixModel", data = "matrix", ...}{A method for calculating
      the GPS values of the data given as 0-1 \code{matrix}.}
    \item{model = "RtreemixModel", data = "missing", ...}{A method for calculating
      the GPS values of the set of all possible patterns.}    
}}

\arguments{
  \item{model}{An object of the class \code{RtreemixModel} specifying
    the mutagenetic trees mixture model used for deriving the GPS values.}  
  \item{data}{An \code{RtreemixData} object or a 0-1 \code{matrix}
    containing the samples (patterns of genetic events) for which the GPS values
    are to be calculated. The length of each of them has to be equal
    to the number of genetic events in the \code{model}.}
  \item{...}{ 
    \code{sampling.mode} is a \code{character} that specifies the
    sampling mode ("constant" or "exponential") used in the waiting time
    simulations. Its default value is "exponential".
    \code{sampling.param} is a \code{numeric} that specifies the
    sampling parameter corresponding to the sampling mode given by
    \code{sampling.mode}. Its default value is 1.
    \code{no.sim} is an \code{integer} larger than 0 giving the number of
    iterations for the waiting time simulations. Its default value is 10.
    \code{seed} is a positive \code{integer} specifying the random generator
    seed. Its default value is (-1) and then the time is used as a
    random generator.
  }
}

\value{
  The function returns an object from the \code{RtreemixGPS} class that
  containes the calculated GPS values, the model used for the
  computation, the data, and so on (see
  \code{\link{RtreemixGPS-class}}). The GPS values are represented as a
  \code{numeric} vector with length equal to the number of samples in \code{data}.
}

\references{Estimating cancer survival and clinical outcome based on
  genetic tumor progression scores, J. Rahnenf\"urer et al. }

\note{
    The GPS examples are time consuming. They are commented out because of the time restrictions of the check of the package.
    For trying out the code please copy it and uncomment it.
} 

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixGPS-class}}, \code{\link{RtreemixData-class}},
  \code{\link{RtreemixModel-class}},
  \code{\link{fit-methods}}, \code{\link{confIntGPS-methods}}
}

\examples{
## Create an RtreemixData object from a randomly generated RtreemixModel object.
#rand.mod <- generate(K = 2, no.events = 7, noise.tree = TRUE, prob = c(0.2, 0.8))
#data <- sim(model = rand.mod, no.draws = 400)

## Create an RtreemixModel object by fitting model to the given data.
#mod <- fit(data = data, K = 2, equal.edgeweights = TRUE, noise = TRUE)
#show(mod)

## Create an RtreemixGPS object by calculating the GPS for all possible patterns.
#modGPS.all <- gps(model = mod, no.sim = 1000) ## time consuming copmutations
#show(modGPS.all)

## See the GPS values for all possible data.
#GPS(modGPS.all) ## time consuming copmutations

## Create an RtreemixGPS object by calculating the GPS for the data based on the model mod.
#modGPS <- gps(model = mod, data = data, no.sim = 1000)
#show(modGPS) ## time consuming copmutations

## See the GPS values for data.
#GPS(modGPS) ## time consuming copmutations
}


\keyword{methods}
\keyword{survival}
