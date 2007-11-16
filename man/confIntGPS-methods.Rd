\name{confIntGPS-methods}
\docType{methods}

\alias{confIntGPS}
\alias{confIntGPS-methods}
\alias{confIntGPS,RtreemixData,numeric-method}


\title{Method for calculating GPS values and their 95\% bootstrap
  confidence intervals}

\description{
  The method first calculates the genetic progression score (GPS) for the
  patterns in a given dataset \code{data} based on a fitted mutagenetic trees
  mixture model with \code{K} components. The \code{data} and \code{K}
  have to be specified. Then, it derives a 95\% confidence intervals for
  the GPS values with bootstrap analysis.   
}

\usage{
\S4method{confIntGPS}{RtreemixData,numeric}(data, K, ...)
}

\arguments{
  \item{data}{An \code{RtreemixData} object containing the samples
    (patterns of genetic events) for which the GPS values and their
    bootstrap confidence intervals are to be calculated.}
  \item{K}{An \code{integer} larger than 0 specifying the number of
    branchings in the mixture model.}
  \item{...}{
    \code{sampling.mode} is a \code{character} that specifies the
      sampling mode ("constant" or "exponential") used in the waiting time
      simulations. Its default value is "exponential".
    \code{sampling.param} is a \code{numeric} that specifies the
      sampling parameter corresponding to the sampling mode given by
      \code{sampling.mode}. Its default value is 1.
    \code{no.sim} is an \code{integer} larger than 0 giving the number of
      iterations for the waiting time simulation. Its default values is
      10000.
    \code{B} is an \code{integer} larger than 0 specifying the number of
      bootstrap samples used in the bootstrap analysis. Its default value
      is 1000.
    \code{equal.star} is a \code{logical} specifying whether to use
      equal edge weights in the noise component. The default value is
      \code{TRUE}.
  }
}

\value{
  The function returns an object from the \code{RtreemixGPS} class that
  containes the calculated GPS values, their 95\% confidence intervals,
  the model used for the computation, the data, and so on (see
  \code{\link{RtreemixGPS-class}}). The GPS values are represented as a
  \code{numeric} vector with length equal to the number of samples in
  \code{data}. Their corresponding confidence intervals are given in a
  matrix with two columns.
}

\note{
    The GPS examples are time consuming. They are commented out because of the time restrictions of the check of the package.
    For trying out the code please copy it and uncomment it.
} 

\author{Jasmina Bogojeska }

\seealso{
  \code{\link{RtreemixGPS-class}},  \code{\link{gps-methods}},
  \code{\link{RtreemixData-class}}, \code{\link{RtreemixModel-class}},
  \code{\link{fit-methods}}
}

\examples{
## Create an RtreemixData object from a randomly generated RtreemixModel object.
#rand.mod <- generate(K = 2, no.events = 7, noise.tree = TRUE, prob = c(0.2, 0.8))
#data <- sim(model = rand.mod, no.draws = 400)

## Create an RtreemixGPS object by calculating GPS values for a given dataset
## and their 95\% confidence intervals using the bootstrap method.
#modGPS2 <- confIntGPS(data = data, K = 2, B = 100) ## time consuming computation
#show(modGPS2)

## See the GPS values for the object modGPS2 and their confidence intervals.
#GPS(modGPS2)
#gpsCI(modGPS2)

## See data.
#getData(modGPS2)
}

\keyword{methods}
