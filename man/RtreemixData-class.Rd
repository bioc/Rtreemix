\name{RtreemixData-class}
\docType{class}
\alias{RtreemixData-class}


\alias{Description}
\alias{Events}
\alias{Patients}
\alias{Sample}
\alias{eventsNum}
\alias{sampleSize}
\alias{Description<-}
\alias{Events<-}
\alias{Patients<-}

\alias{Description,RtreemixData-method}
\alias{Description<-,RtreemixData-method}
\alias{Events,RtreemixData-method}
\alias{Events<-,RtreemixData-method}
\alias{Patients,RtreemixData-method}
\alias{Patients<-,RtreemixData-method}
\alias{Sample,RtreemixData-method}
\alias{eventsNum,RtreemixData-method}
\alias{sampleSize,RtreemixData-method}
\alias{initialize,RtreemixData-method}
\alias{print,RtreemixData-method}
\alias{show,RtreemixData-method}





\title{Class "RtreemixData"}

\description{
  This class is used to represent the results of genetic
  measurements of the occurence of subsets of a given set of genetic
  events in a group of patients. Each observation is a binary vector
  that indicates which events occured in a specific patient. The length
  of the vector equals the size of the set of genetic events that is
  taken into consideration. 
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("RtreemixData",
    Sample, Patients, Events, Description, File)}.
  The \code{RtreemixData} class represents patterns of occurences of subsets of
  a given set of genetic events in a specific group of patients. The
  patterns are given as binary vectors with length equal to the size of
  the set of genetic events. In other words, it provides a
  representation of the dataset used for learning an mutagenetic trees
  mixture model.  
  
  The \code{Sample} is a binary \code{matrix} where each row corresponds to the pattern of
  genetic events observed in one of the given patients. Hence, the
  number of rows gives the number of patients, i.e. the size of the
  dataset. Each column corresponds to one of the genetic events.  
  Missing measurement for the presence or absence of a certain genetic
  event in a given pattern is marked with -1. The initial null event
  (that initially occurs in all patients) is
  not present in the sample, i.e. the first component in each
  observation (which is always equal to 1) is left out. This is done for
  saving space and avoiding the process of checking for correctly
  specified samples. 

  The \code{Patients} is a character \code{vector} that contains the IDs
  of the patients. The length of this vector must be equal to the number
  of rows in the \code{Sample}.

  The \code{Events} is a character \code{vector} that contains the
  labels of the genetic events taken into consideration. Its length
  equals one plus the number of columns in the \code{Sample}. This is
  because of the label of the null event. When the object of class
  \code{RtreemixData} is a parent of a randomly generated
  \code{RtreemixModel} object, the events specify the labels of the
  genetic events present in the random model, although the \code{Sample} slot
  is an empty \code{matrix}. This is because the random mixture models
  are not estimated from a given dataset, but generated randomly for
  some set of genetic events.

  The \code{Description} is a character giving a short description for
  the created object.
  
  The \code{File} specifies the path to a text file with a specific
  format which contains the infromation needed to create an
  \code{RtreemixData} object (the patient IDs, the names of the events,
  the matrix with the observations).    
}

\section{Slots}{
	 \describe{
    \item{\code{Sample}:}{Object of class \code{"matrix"}.}
    \item{\code{Patients}:}{Object of class \code{"character"}. The
      \code{Patients} must be of same length as the number of rows in \code{Sample}.}
    \item{\code{Events}:}{Object of class \code{"character"}. The length
      of \code{Events} must be identical to the number of columns in
      \code{Sample} plus one (for specifying the label of the null event).  }
    \item{\code{Description}:}{Object of class \code{"character"}.}
  }
}

\section{Methods}{
  \describe{
    \item{Description}{\code{signature(object = "RtreemixData")}: A
      method for obtaining the description of the \code{"RtreemixData"} object.}
    \item{Description<-}{\code{signature(object = "RtreemixData")}: A
      method for specifying the \code{Description} of the data object. }
    \item{Events}{\code{signature(object = "RtreemixData")}: A
      method for obtaining the labels of the genetic events. }
    \item{Events<-}{\code{signature(object = "RtreemixData")}: A method
      for replacing the names of the genetic events in the data object. It checks
      to be sure the values have the right length. As a parent data
      of a random \code{RtreemixModel} object the suitable labels of
      events present in the model components can be specified although the
      \code{Sample} slot is an empty \code{matrix}.}
    \item{Patients}{\code{signature(object = "RtreemixData")}: A
      method for obtaining the IDs of the patients. }
    \item{Patients<-}{\code{signature(object = "RtreemixData")}: A method
      for replacing the IDs of the patients in the data object. It checks
      to be sure the values have the right length. }
    \item{Sample}{\code{signature(object = "RtreemixData")}: A
      method for obtaining the matrix of observations. }
    \item{eventsNum}{\code{signature(object = "RtreemixData")}: A
      method for obtaining the number of genetic events. }
    \item{sampleSize}{\code{signature(object = "RtreemixData")}: A
      method for obtaining the size of the sample (the number of
      patients). }
	 }
}
 
\author{Jasmina Bogojeska}
 
\seealso{
  \code{\link{RtreemixGPS-class}}, \code{\link{RtreemixStats-class}},
  \code{\link{RtreemixModel-class}},
  \code{\link{fit-methods}}, \code{\link{bootstrap-methods}}
}


\examples{
## Create an RtreemixData object from a file given in the examples directory of the package.
data1 <- new("RtreemixData", File = paste(system.file(package = "Rtreemix"), "/examples/treemix.pat", sep = ""))
show(data1) ## show the RtreemixData object

## Create an RtreemixData object from a randomly generated RtreemixModel object.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
data2 <- sim(model = rand.mod, no.draws = 300)
show(data2)

## Create an RtreemixData object from a given binary matrix.
bin.mat <- cbind(c(1, 0, 0, 1, 1), c(0, 1, 0, 0, 1), c(1, 1, 0, 1, 0))
data3 <- new("RtreemixData", Sample = bin.mat, Events = c("0", "1", "2",
"3"))
show(data3)
}
\keyword{classes}
