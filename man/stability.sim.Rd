\name{stability.sim}
\alias{stability.sim}

\title{Stability analysis of the mutagenetic trees mixture model}
\description{
  The function includes stability analysis on different levels of the mutagenetic trees mixture
  model: GPS values, encoded probability distribution, tree topologies.
  Each analysis contains the values of different similarity measures
  with their corresponding p-values.
}
\usage{
stability.sim(no.trees = 3, no.events = 9, prob = c(0.2, 0.8),
              no.draws = 300, no.rands = 100, no.sim = 1)
}
\arguments{
  \item{no.trees}{An \code{integer} larger than 2 giving the number of tree components
    of the mixture models considered in the stability analysis. The
    default value is 3.}
  \item{no.events}{An \code{integer} larger than 0 giving the number of genetic events of the mixture models
    considered in the stability analysis.}
  \item{prob}{A \code{numeric} vector of length 2 specifying the
    boundaries for the edge weights of the randomly generated trees. The
    first component of the vector (the lower boundary) must be smaller
    than the second component (the upper boundary). The default value is (0.2, 0.8).}
  \item{no.draws}{An \code{integer} larger than 0 giving the size of the
    data sample drawn from the random models used for learning the
    mixture models. The default value is 300.}
  \item{no.rands}{An \code{integer} larger than 0 specifying the number
    of random models used for calculating the p-values. The default value is 100.}
  \item{no.sim}{An \code{integer} larger than 0 specifying the number of
  iterations used for the waiting time simulations (a part of the GPS
  calculation). The default value is 1.}
}
\details{
  The stability analysis is performed by first drawing a true mixture
  model uniformly at random from the model space, and drawing a data
  sample from it. Afterwards, a mutagenetic trees model is fitted to the
  drawn sample. The quality of the features derived from the model is
  then assessed by comparing its quality with the quality of the
  corresponding features of a sufficient number of random mixture models
  sampled uniformly from the model space. A p-value is obtained as a
  percentage of cases in which the true model is closer to a random
  model tnah to the fitted model. 
}
\value{
  \item{comp1 }{Results from the stability analysis of the GPS values
    derived from a fitted mixture model.
    A \code{matrix} with 4 columns and \code{no.sim} rows. The first two columns give the similarity
    values and their corresponding p-values when the Euclidian distance
    is used as a similarity measure for comparing the respective GPS
    vectors. The last two columns depict the same results, but with the
    rank correlation distance used as a similarity measure.
  }
  \item{comp2 }{Results from the stability analysis of the probability
    distributions induced by a fitted mixture model. A \code{matrix} with 6
    columns and \code{no.sim} rows. Each two columns give the values of the comparissons between
    the true and the fitted probability distributions and their
    corresponding p-values, when using the cosine distance, the L1 distance, and the
    Kullback-Leibler divergence as similarity measures.}
  \item{comp3 }{Results from the stability analysis of the topologies
    of the tree components of a fitted mixture model. A \code{matrix} with 2
    columns and \code{no.sim} rows that give the value of the comparisson of the topologies
    between the true and the corresponding fitted model and their
    p-values. The similarity measure underlying the number of different
    edges was used.}
  \item{comp4 }{Similar to \code{comp3}. However, the similarity measure for
    comparing the tree topologies besides the number of distinct edges
    includes the L1 distances of the level vectors of events. See
    \code{\link{get.tree.levels}}.}
  \item{comp5 }{A \code{matrix} where the columns correspond to the
    true GPS vector from each simulation iteration. The matrix has
    \code{no.sim} columns and \code{no.draws} rows.}
  \item{comp6 }{Same as \code{comp5}, but the matrix contains the
    fitted GPS values from each simulation iteration.}
  \item{comp7 }{A \code{list} where each component corresponds to the
    true models generated in each simulation iteration. the length of
    the list is \code{no.sim}.}
  \item{comp8 }{Same as \code{comp7}, but the list contains the fitted models.}
}

\references{Learning multiple evolutionary pathways from cross-sectional
  data, N. Beerenwinkel et al.; Estimating cancer survival and clinical outcome based on
  genetic tumor progression scores, J. Rahnenf\"urer et al.}

\note{
    The stability simulation examples are time consuming. They are commented out because of the time restrictions of the check of the package.
    For trying out the code please copy it and uncomment it.
} 
 
\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixData-class}}, \code{\link{RtreemixModel-class}}, \code{\link{RtreemixGPS-class}},
  \code{\link{RtreemixStats-class}}, \code{\link{fit-methods}},
  \code{\link{gps-methods}}, \code{\link{distribution-methods}},
  \code{\link{generate-methods}}, \code{\link{sim-methods}},
  \code{\link{L1.dist}}, \code{\link{Pval.dist}},
  \code{\link{comp.models}}, \code{\link{comp.trees}},
  \code{\link{get.tree.levels}}, \code{\link{kullback.leibler}}
  
}
\examples{
## Stability analysis - a toy example
#stability.sim(no.trees = 3, no.rands = 5, no.sim = 4, no.draws = 300)
}

\keyword{misc}

