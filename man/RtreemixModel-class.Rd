\name{RtreemixModel-class}
\docType{class}

\alias{CompleteMat}
\alias{Resp}
\alias{Star}
\alias{Trees}
\alias{Weights}
\alias{WeightsCI}
\alias{getTree}
\alias{getData,RtreemixModel}
\alias{numTrees}
\alias{Weights<-}


\alias{RtreemixModel-class}
\alias{CompleteMat,RtreemixModel-method}
\alias{Resp,RtreemixModel-method}
\alias{Star,RtreemixModel-method}
\alias{Trees,RtreemixModel-method}
\alias{Weights,RtreemixModel-method}
\alias{Weights<-,RtreemixModel-method}
\alias{WeightsCI,RtreemixModel-method}
\alias{getData,RtreemixModel-method}
\alias{getTree,RtreemixModel,numeric-method}
\alias{numTrees,RtreemixModel-method}
\alias{initialize,RtreemixModel-method}
\alias{print,RtreemixModel-method}
\alias{show,RtreemixModel-method}
\alias{plot,RtreemixModel-method, missing-method}

\title{Class "RtreemixModel"}
\description{
  This class contains all the data needed for characterizing the mutagenetic
  trees mixture model (mixture parameters, mixture components, ...).
  The tree components of the model are given as a list of directed \code{graphNEL} objects.   
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RtreemixModel",
    ParentData, Weights, WeightsCI, Resp, CompleteMat, Star, Trees)}.
  The \code{RtreemixModel} class extends the \code{RtreemixData} class
  and specifies the mutagenetic trees mixture model. If the model is not
  randomly generated the parent class gives the \code{RtreemixData}
  used for learning the mixture model. The directed trees that build up
  the model are represented as a list of directed \code{graphNEL}
  objects, and their weights (the mixture parameters) are given as a
  numeric vector. This class can also contain other useful information
  connected with the mixture model like confidence intervals for the
  mixture parameters and the edge weights (resulting from a bootstrap
  analysis), an indicator for the presence of the star component, etc.
  They are all listed in the text below with brief descriptions.

  The \code{ParentData} is an \code{RtreemixData} object that specifies the
  data used for estimating the mutagenetic trees mixture model. It is
  not specified for random mixture models, since they are not estimated
  from a given dataset but generated randomly.

  The \code{Weights} is a numeric \code{vector} that contains the mixture
  parameters of the model. Its length equals the length of the
  \code{list} of tree components \code{Trees}.

  The \code{WeightsCI} is a named \code{list} with length equal to the
  length of the \code{Weights}. Each list element is a numeric
  \code{vector} of length two specifying the lower and upper bound of
  the confidence interval for the corresponding mixture parametar. The
  confidence intervals are derived using the bootstrap method.

  The \code{Resp} is a numeric \code{matrix} that specifies the responsibility
  of each tree component to generate each of the patterns in the
  \code{ParentData}. The number of rows in \code{Resp} is equal
  to the number of trees in the mixture (the length of the list
  \code{Trees}) and the number of columns equals the number of patients
  in \code{ParentData}. For random mixture models it is an empty matrix,
  since they are not estimated from a given dataset.

  The \code{CompleteMat} is a binary \code{matrix} that specifies the complete
  data in case some measurements for some patients are missing in
  the data used for learning the model (the \code{ParentData}). It has
  the same size as the matrix specifying the data in \code{ParentData}.
  The missing data are estimated in the EM-algorithm used for fitting
  the mixture model. When there are no missing data in
  \code{ParentData}, or the model is randomly generated the \code{CompleteMat} is an
  empty matrix. 

  The \code{Star} is an indicator of the presence of a noise (star) component
  and is mostly relevant for models with a single tree component, since it is assumed that 
  mixture models with at least two components always have the noise
  as a first component. It is of type \code{logical}.

  The \code{Trees} is a \code{list} of directed \code{graphNEL}
  objects, each for every tree component in the mixture model. The
  genetic events are represented as nodes in the graphs. The
  \code{edgeData} of each tree can have two attributes: \code{"weight"}
  and \code{"ci"}. Please see the help page on \code{graph-class} and
  \code{graphNEL-class} in the package \code{graph}. The \code{"weight"} attribute is for edge weight,
  i.e. the conditional probability that the child event of the edge occured given
  that the parent event already occured. The \code{"ci"} attribute is
  for the bootstrap confidence intervals for the edge weight (a numeric vector
  with length two).     
}
\section{Slots}{
     \describe{
    \item{\code{Weights}:}{Object of class \code{"numeric"}. The length
      of the \code{Weights} must be equal to the length of \code{Trees}. }
    \item{\code{WeightsCI}:}{Object of class \code{"list"}. The length
      of the \code{WeightsCI} must be equal to the length of \code{Weights}.  }
    \item{\code{Resp}:}{Object of class \code{"matrix"}. The number of
      rows of \code{Resp} must be identical to the length of
      \code{Trees}, and its number of columns to the number of patients
      in the dataset used for estimating the mixture model (\code{ParentData}).}
    \item{\code{CompleteMat}:}{Object of class \code{"matrix"}. When
      specified (when there are missing data) the size of the
      \code{CompleteMat} must be equal to the size of the matrix used to
      estimate the model.}
    \item{\code{Star}:}{Object of class \code{"logical"}. }
    \item{\code{Trees}:}{Object of class \code{"list"}. The length of
      \code{Trees} equals the length of \code{Weights}.}
  }
}
\section{Extends}{
Class \code{"RtreemixData"}, directly.
}
\section{Methods}{
  \describe{
    \item{CompleteMat}{\code{signature(object = "RtreemixModel")}: A
      method used for obtaining the complete dataset, in case there were
      any missing measurements for some patients in the dataset used to
      learn the mixture model.}
    \item{Resp}{\code{signature(object = "RtreemixModel")}: A method for
    obtaining the matrix of responsibilities for the trees to generate
    each of the samples in the dataset used for learning the model (\code{ParentData}).}    
    \item{Star}{\code{signature(object = "RtreemixModel")}: A method for
    checking the presence of a noise component in the mixture model
    (informative only for models with one tree component).}   
    \item{Trees}{\code{signature(object = "RtreemixModel")}: A method
      for obtaining the tree components of the mixture model as a list
      of directed \code{graphNEL} objects. }
    \item{Weights}{\code{signature(object = "RtreemixModel")}: A method
      for obtaining the mixture parameters (the weights of the trees in
      the model).}
    \item{Weights<-}{\code{signature(object = "RtreemixModel")}: A
      method for replacing the value of the slot \code{Weights} with a
      specified \code{numeric} vector. The components of this vector
      have to sum up to one.}
    \item{WeightsCI}{\code{signature(object = "RtreemixModel")}: A method
      for obtaining the weights of the mixture parameters. }
    \item{getData}{\code{signature(object = "RtreemixModel")}: A method
      for obtaining the \code{ParentData} of the mixture model, i.e. the
      data used for learning the model. }
    \item{getTree}{\code{signature(object = "RtreemixModel", k =
    "numeric")}: A method for obtaining the k-th tree component of the
      mixture model as a directed \code{graphNEL} object.}
    \item{numTrees}{\code{signature(object = "RtreemixModel")}: A method
    for obtaining the number of tree components building up the mixture model.}
  }
}
\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.}


\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixGPS-class}}, \code{\link{RtreemixStats-class}},
  \code{\link{RtreemixData-class}}, \code{\link{RtreemixSim-class}},
  \code{\link{fit-methods}}, \code{\link{bootstrap-methods}},
  \code{\link{generate-methods}}, \code{\link{comp.models}}, \code{\link{comp.trees}}
}

\examples{
## Generate a random RtreemixModel object with 2 components.
rand.mod <- generate(K = 2, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
show(rand.mod)
plot(rand.mod) ## plot the tree components of the model
plot(rand.mod, k = 2) ## plot the second component of the model

## Draw data from a specified mixture model.
draws <- sim(model = rand.mod, no.draws = 200)
show(draws)

## Create an RtreemixModel object by fitting model to the drawn data.
mod <- fit(data = draws, K = 2, equal.edgeweights = TRUE, noise = TRUE)
show(mod)

## See the values of the slots of the RtreemixModel object.
Weights(mod)
Resp(mod)
CompleteMat(mod)
Star(mod)
Trees(mod)
## See data used for learning the model.
getData(mod)
## See the number of tree components in the mixture model.
numTrees(mod)
## See a specific tree component k.
getTree(object = mod, k = 2)
## See the conditional probabilities assigned to edges of the second tree component.
edgeData(getTree(object = mod, k = 2), attr = "weight")
## See the probability distribution encoded by the model on the set of all possible patterns.
distr <- distribution(model = mod)
distr
## Get the probabilities.
distr$probability
## See the probability distribution encoded by the model on the set of all possible patterns
## calculated for given sampling mode, and input and output parameters.
distr1 <- distribution(model = mod, sampling.mode = "exponential", sampling.param = 1, output.param = 1)
distr1

## Create a RtreemixModel and analyze its variance with the bootstrap method.
mod.boot <- bootstrap(data = draws, K = 2, equal.edgeweights = TRUE, B = 100)

## See the confidence intervals for the mixture parameters (the weights).
WeightsCI(mod.boot)
## See the confidence intervals of the conditional probabilities assigned to the edges.
edgeData(getTree(mod.boot, 2), attr = "ci")
}
\keyword{classes}
