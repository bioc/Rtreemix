\name{bootstrap-methods}
\docType{methods}


\alias{bootstrap}
\alias{bootstrap-methods}
\alias{bootstrap,RtreemixData,numeric-method}

\title{Method for fitting a mutagenetic trees mixture model and analyzing its variance}
\description{
This method fits an \code{RtreemixModel} to a given dataset and then
analyzes its variance with the bootstrap method. The \code{data} and
the number of trees \code{K} have to be specified. 
}
\usage{
\S4method{bootstrap}{RtreemixData,numeric}(data, K, ...)
}

\arguments{
  \item{data}{An \code{RtreemixData} object giving the dataset used for
    learning the trees mixture model. }
  \item{K}{An \code{integer} larger than 0 specifying the number of
    branchings in the mixture model.}
  \item{...}{  
    \code{no.start.sol} is an \code{integer} larger than 0 specifying the number of starting solutions for the k-means
    algorithm. The default value is 100.
    \code{eps} is a \code{numeric} giving the minimum conditional probability to include edge. The
    default value is 0.
    \code{weighing} is a \code{logical} specifying whether to use special
      weights log(Pr(v)) for the edges (root, v). The default value is \code{FALSE}.
    \code{equal.edgeweights} is a \code{logical} specifying whether to use
      equal edge weights in the noise component. The default value is
      \code{TRUE}.
    \code{seed} is a positive \code{integer} specifying the random generator
      seed. The default value is (-1) and then the time is used as a
      random generator.
    \code{B} is an \code{integer} larger than 0 specifying the number of
      bootstrap samples. Its default value is 1000.
    \code{conf.interval} is a \code{numeric} specifying the Confidence level
      for the intervals. Its default value is 0.05.
  }
}

\value{
  The function returns an object from the class
  \code{RtreemixModel}. This is the mixture model learned on the given
  \code{data}. Besides the edge weights it also contains their confidence
  intervals resulting from the bootstrap analysis. Confidence intervals
  for the mixture parameters are also comupted and available.  
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}

\note{
    The bootstrap examples are time consuming. They are commented out because of the time restrictions of the check of the package.
    For trying out the code please copy it and uncomment it.
}    

\author{Jasmina Bogojeska }

\seealso{\code{\link{RtreemixData-class}},
  \code{\link{RtreemixModel-class}}, \code{\link{fit-methods}} 
}

\examples{
## Create an RtreemixData object from a randomly generated RtreemixModel object.
#rand.mod <- generate(K = 2, no.events = 7, noise.tree = TRUE, prob = c(0.2, 0.8))
#data <- sim(model = rand.mod, no.draws = 300)

## Create a RtreemixModel and analyze its variance with the bootstrap method.
#mod.boot <- bootstrap(data = data, K = 2, equal.edgeweights = TRUE, B = 10) ## time consuming computation

## See the confidence intervals for the mixture parameters (the weights).
#WeightsCI(mod.boot)
## See the confidence intervals of the conditional probabilities assigned to the edges.
#edgeData(getTree(mod.boot, 2), attr = "ci")
}

\keyword{methods}
\keyword{datagen}
