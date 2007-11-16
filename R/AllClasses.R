##################################################################################
###                             CLASSES DEFINITIONS                            ###
##################################################################################


## Class RtreemixData.
setClass("RtreemixData",
         representation = representation(
           ## Binary patterns (matrix with rows=patients,cols=genetic events), without the null event
           Sample = "matrix",
           ## Patient IDs
           Patients = "character",
           ## Names of genetic events (length L)
           Events = "character",
           ## Description of the object
           Description = "character"),
         prototype = prototype(
           Sample = matrix(integer(0), 0, 0),
           Patients = character(),
           Events = character(),
           Description = character(0)))


################################################################################

## Class RtreemixModel. Extends the class RtreemixData.
setClass("RtreemixModel",
         representation = representation(
           ## The weight vector of the model
           Weights = "numeric",
           ## Confidence intervals for the weights (from bootstrap analysis)
           WeightsCI = "list",
           ## The responsibilities
           Resp = "matrix",
           ## The complete sample matrix (if there were some missing data otherwise empty)
           CompleteMat = "matrix",
           ## Indicator of the presence of a star component (mostly relevant for models with a single tree component)
           Star = "logical",
           ## The list of the graphs each for every tree component of the mixture model
           Trees = "list" ),
         prototype = prototype(
           Weights  = numeric(0),
           WeightsCI = list(),
           Resp = matrix(numeric(0), 0, 0),
           CompleteMat = matrix(integer(0), 0, 0),  
           Star = logical(0),
           Trees = list()),
         contains = "RtreemixData")


################################################################################

## Class RtreemixSim. Extends the class RtreemixModel.
setClass("RtreemixSim",
         representation = representation(
           ## Data drawn (or simulated in a waiting time simulation) from an oncogenetic trees mixture model  
           SimPatterns = "RtreemixData",
           ## Sampling mode for the simulations: exponential or constant
           SamplingMode = "character",
           ## Sampling parameter that corresponds to the sampling mode
           SamplingParam = "numeric",
           ## Waiting times of the simulated patterns
           WaitingTimes = "numeric",
           ## Sampling times of the simulated patterns
           SamplingTimes = "numeric"
           ),
         prototype = prototype(
           SimPatterns = new("RtreemixData", Sample = matrix(integer(0), 0, 0)),
           SamplingMode = character(0),
           SamplingParam = numeric(0),
           WaitingTimes  = numeric(0),
           SamplingTimes = numeric(0)  
           ),
         contains = "RtreemixModel")

################################################################################

## Class RtreemixStats. Extends the class RtreemixData.
setClass("RtreemixStats",
         representation = representation(
           ## The underlying model for calculating the (log, weighted) likelihoods.
           Model = "RtreemixModel",
           ## The log-likelihoods of the set of patterns
           LogLikelihoods = "numeric",
           ## The weighted likelihoods for the set of patterns
           WLikelihoods = "matrix"),
         prototype = prototype(
           Model = new("RtreemixModel", Weights = numeric(0), Trees = list()),
           LogLikelihoods = numeric(0), 
           WLikelihoods  = matrix(numeric(0), 0, 0)
         ),
         contains = "RtreemixData")

################################################################################

## Class RtreemixGPS. Extends the class RtreemixModel.
setClass("RtreemixGPS",
         representation = representation(
           ## Underlying model for the GPS is calculation. 
           Model = "RtreemixModel",
           ## Sampling mode for the simulations of the waiting time process: exponential or constant
           SamplingMode = "character",
           ## Sampling parameter that corresponds to the sampling mode
           SamplingParam = "numeric",
           ## GPS vector associated to the corresponding dataset of patterns 
           GPS = "numeric",
           ## Confidence intervals for the GPS values (from bootstrap analysis)
           gpsCI = "matrix"           
           ),
         prototype = prototype(
           Model = new("RtreemixModel", Weights = numeric(0), Trees = list()),
           SamplingMode = character(0),
           SamplingParam = numeric(0),
           GPS  = numeric(0),
           gpsCI = matrix(numeric(0), 0, 0)
           ),
         contains = "RtreemixData")

################################################################################
