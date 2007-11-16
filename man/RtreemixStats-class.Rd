\name{RtreemixStats-class}
\docType{class}

\alias{LogLikelihoods}
\alias{Model,RtreemixStats}
\alias{WLikelihoods}
\alias{getData,RtreemixStats}
\alias{getResp}

\alias{RtreemixStats-class}
\alias{LogLikelihoods,RtreemixStats-method}
\alias{Model,RtreemixStats-method}
\alias{WLikelihoods,RtreemixStats-method}
\alias{getData,RtreemixStats-method}
\alias{getResp,RtreemixStats-method}
\alias{initialize,RtreemixStats-method}
\alias{print,RtreemixStats-method}
\alias{show,RtreemixStats-method}




\title{Class "RtreemixStats"}
\description{
  The \code{RtreemixStats} class contains the (weighted, log)
  likelihoods for a given dataset (specified by the parent class)
  derived from the probability distribution induced by an underlying
  mutagenetic trees mixture model.    
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RtreemixStats",
    Data, Model, LogLikelihoods, WLikelihoods)}.
  The class \code{RtreemixStats} extends the \code{RtreemixData} class
  and specifies (log, weighted) likelihoods for these data derived from
  a given \code{RtreemixModel}. The number of the genetic events in the
  patterns from the given dataset (\code{Data}) has to be equal to the number
  of genetic events in the branchings from the mixture model
  given by the slot \code{Model}. When having the weighted likelihoods,
  one can easily derive the responsibilities of the model components in
  \code{Model} for generating the patterns in the specified dataset
  (\code{Data}).

  The \code{Data} is an \code{RtreemixData} object that specifies the
  patterns for which the likelihoods are calculated.

  The \code{Model} is an \code{RtreemixModel} object that specifies the
  mutagenetic trees mixture model used for deriving the likelihoods of
  the given data. 

  The \code{LogLikelihoods} is a numeric \code{vector} that contains the
  log-likelihoods of the patterns in \code{Data}. Its length equals the
  sample size, i.e. the number of patients in \code{Data}.

  The \code{WLikelihoods} is a numeric \code{matrix} that specifies the
  weighted likelihoods of each pattern in the given dataset
  \code{Data}. The number of columns in \code{WLikelihoods} equals the
  number of tree components in \code{Model} and the number of rows
  equals the number of patients in \code{Data}.   
}
\section{Slots}{
	 \describe{
    \item{\code{Model}:}{Object of class \code{"RtreemixModel"}.}
    \item{\code{LogLikelihoods}:}{Object of class \code{"numeric"}. The
      length of \code{LogLikelihoods} must be equal to the number of
      patients of the dataset specified with the parent
      \code{"RtreemixData"} class.}
    \item{\code{WLikelihoods}:}{Object of class \code{"matrix"}. The
      number of rows must be equal to the sample size of the dataset
      specified with the parent \code{"RtreemixData"} class. The number
      of columns must be identical with the number of tree components in
      the mixture model \code{Model}.}
  }
}
\section{Extends}{
Class \code{"RtreemixData"}, directly.
}
\section{Methods}{
  \describe{
    \item{LogLikelihoods}{\code{signature(object = "RtreemixStats")}: A
      method for obtaining the log-likelihoods of the patterns in
      the dataset specified with the parent \code{"RtreemixData"} class.}
    \item{Model}{\code{signature(object = "RtreemixStats")}: A method
      for obtaining the mutagenetic trees mixture model used for
      deriving the likelihoods.}
    \item{WLikelihoods}{\code{signature(object = "RtreemixStats")}: A
      method for obtaining the weighted likelihoods of the patterns in
      the dataset specified with the parent \code{"RtreemixData"} class.}
    \item{getData}{\code{signature(object = "RtreemixStats")}: A method
      for obtaining the dataset specified with the
      parent \code{"RtreemixData"} class.}
    \item{getResp}{\code{signature(object = "RtreemixStats")}: A method for
    computing the matrix of responsibilities for the trees to generate
    each of the samples in the parent dataset from their weighted
    likelihoods \code{WLikelihoods}.}
  }
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}       

\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixData-class}},
  \code{\link{RtreemixModel-class}},
  \code{\link{fit-methods}}, \code{\link{likelihoods-methods}}
}


\examples{
## Generate a random RtreemixModel object with 3 components and 9 genetic events.
mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
show(mod)

## Draw a data sample from the model mod.
data <- sim(model = mod, no.draws = 400)

## Create an RtreemixStats object.
mod.stat <- likelihoods(model = mod, data = data)
show(mod.stat)

## See the slots from the RtreemixStats object.
Model(mod.stat)
LogLikelihoods(mod.stat)
WLikelihoods(mod.stat)
## See data.
getData(mod.stat)
## Calculate the responsibilities from the weighted likelihoods.
getResp(mod.stat)
}
\keyword{classes}
