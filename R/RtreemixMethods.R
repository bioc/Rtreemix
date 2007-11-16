##################################################################################
## Function for fitting Rtreemix model to a given RtreemixData dataset.
## The data and the number of trees (K) have to be specified.
## The function estimates K-oncogenetic trees mixture model from the
## specified data (with a noise component).
## When K == 1, and noise == FALSE a single mutagenetic tree is fit to the data.
## When K == 1 and noise == TRUE a star mutagenetic tree is fit to the data.
## If K > 1, the first mutagenetic tree is always the star,
## i.e. the case K > 1 and noise == FALSE is not possible.
if(!isGeneric("fit"))
  setGeneric("fit", function(data, K, ...) standardGeneric("fit"))

setMethod("fit",
          signature(data = "RtreemixData", K = "numeric"),
          function(data, K, no.start.sol = 100,
                   eps = 0.01, weighing = FALSE, equal.edgeweights = TRUE,
                   seed = (-1), noise = TRUE) {

            ## Setting the true type
            K <- as.integer(K)
            no.start.sol <- as.integer(no.start.sol)
            equal.edgeweights <- as.integer(equal.edgeweights)
            eps <- as.numeric(eps)
            weighing <- as.integer(weighing)
            seed <- as.integer(seed)
            ## Check if the necessary parameters are provided and have the correct form.                
            if(K < 1)
              stop("The number of trees must be an integer greater than zero!")

            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0 )))
              stop("The seed must be a positive integer!")

            if(no.start.sol < 1)
              stop("The number of starting solutions for the k-means must be an integer greater than zero!")            
  
            ##Add the first column of ones in the sample.
            sample <- cbind(rep(1, sampleSize(data)), Sample(data))

            ## Fitting the corresponding RtreemixModel to the provided data.
            if(noise && (K > 1)) {
              fit.res <- .Call("R_fit", Events(data), sample, K,
                               no.start.sol, equal.edgeweights,
                               eps, weighing, seed)
    
              for(i in 1:K)
                edgeDataDefaults(fit.res$graphs.mixture[[i]], "weight") <- numeric(1)
              
              return(new("RtreemixModel",
                         ParentData = data,
                         Weights =  as.numeric(fit.res$alpha),
                         Resp = as.matrix(fit.res$resp),
                         CompleteMat = as.matrix(fit.res$pat.hat),
                         Star = TRUE,
                         Trees = fit.res$graphs.mixture))
            } else {
              if(K == 1) {
                if (!noise) {
                  ## Fit single RtreemixModel to the data.
                  fit.res <- .Call("R_fit1", Events(data), sample, eps, weighing)
                  
                  edgeDataDefaults(fit.res$graphs.mixture[[1]], "weight") <- numeric(1)
                  
                  return(new("RtreemixModel",
                             ParentData = data,
                             Weights =  as.numeric(fit.res$alpha),
                             Resp = as.matrix(t(fit.res$resp)),
                             Star = FALSE,
                             Trees = fit.res$graphs.mixture))
                } else {
                  ## Fit single star RtreemixModel to the data.
                  fit.res <- .Call("R_fit0", Events(data), sample,
                                   equal.edgeweights, weighing)
                  
                  edgeDataDefaults(fit.res$graphs.mixture[[1]], "weight") <- numeric(1)
                  
                  return(new("RtreemixModel",
                             ParentData = data,
                             Weights =  as.numeric(fit.res$alpha),
                             Resp = as.matrix(t(fit.res$resp)),
                             Star = TRUE,
                             Trees = fit.res$graphs.mixture))
                }
              } else
              stop("Mixture model with more than one tree and without noise component is not implemented!")
            }
          })

##################################################################################
## Function for fitting Rtreemix model to given data and
## analyzing its variance with the bootstrap method.
## The data and number of trees (K) have to be specified.
if(!isGeneric("bootstrap"))
  setGeneric("bootstrap", function(data, K, ...) standardGeneric("bootstrap"))

setMethod("bootstrap",
          signature(data = "RtreemixData", K = "numeric"),
          function(data, K, no.start.sol = 100,
                      eps = 0.0, weighing = FALSE, equal.edgeweights = TRUE,
                      seed = (-1), B = 1000,
                      conf.interval = 0.05) {

            ## Setting the true type
            K <- as.integer(K)
            no.start.sol <- as.integer(no.start.sol)
            equal.edgeweights <- as.integer(equal.edgeweights)
            eps <- as.numeric(eps)
            weighing <- as.integer(weighing)
            seed <- as.integer(seed)
            B <- as.integer(B)
            conf.interval <- as.numeric(conf.interval)
            events <- Events(data)
            ## Check if the necessary parameters are provided and have the correct form.    
            if(K < 1)
              stop("The number of trees must be an integer greater than zero!")
  
            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0 )))
              stop("The seed must be a positive integer!")
            
            if(no.start.sol < 1)
              stop("The number of starting solutions for the k-means must be an integer greater than zero!")

            if(B < 1)
              stop("The number of bootstrap samples integer greater than zero!")            
            
            ## Add the first column of ones in the sample.
            sample <- cbind(replicate(sampleSize(data), 1), Sample(data))

            ## Fit a model and analyze its variance with the bootstrap method.
            fit.boot <- .Call("R_bootstrap", events, sample, K,
                              no.start.sol, equal.edgeweights,
                              eps, weighing, seed,
                              B, conf.interval
                              )

            for(i in 1:K) {
              edgeDataDefaults(fit.boot$graphs.mixture[[i]], "weight") <- numeric(1)
              edgeDataDefaults(fit.boot$graphs.mixture[[i]], "ci") <- numeric(0)
            }
            
            return(new("RtreemixModel",
                       ParentData = data,
                       Weights =  as.numeric(fit.boot$alpha),
                       WeightsCI = fit.boot$alpha.ci,
                       Resp = as.matrix(fit.boot$resp),
                       CompleteMat = as.matrix(fit.boot$pat.hat),
                       Star = TRUE,
                       Trees = fit.boot$graphs.mixture
                       ))  
          })

####################################################################################
## Function for predicting the likelihoods of given data by using a given Rtreemix model.
## The model and the data have to be specified.
if(!isGeneric("likelihoods"))
  setGeneric("likelihoods", function(model, data) standardGeneric("likelihoods"))

setMethod("likelihoods",
          signature(model = "RtreemixModel", data = "RtreemixData"),
          function(model, data) {
            ## Check if the patterns have the correct length.
            if(eventsNum(model) !=  eventsNum(data)) {
              cat("The patterns are expected to have length: ", eventsNum(model), ". ")
              stop("\n")
            } else {
              ## Give the node names and edge weights in a list structure.             
              listG <- list()
              for(i in 1:numTrees(model)) {
                listG[[i]] <- list()
                listG[[i]][[1]] <- nodes(Trees(model)[[i]])
                listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
              }
              
              ## Call to the C function for calculating the likelihoods.              
              res.like <- .Call("R_likelihood", eventsNum(model),
                                listG, Weights(model), cbind(rep(1, sampleSize(data)), Sample(data)))
    
              return(new("RtreemixStats",
                         Data = data,
                         Model = model,
                         LogLikelihoods = res.like[[2]],
                         WLikelihoods = res.like[[1]]))  
            }
          })

###################################################################################
## Function for predicting the GPS of given dataset by using a given Rtreemix model.
## The model has to be specified. If the dataset is missing the GPS is calculated for all possible patterns.
if(!isGeneric("gps"))
  setGeneric("gps", function(model, data, ...) standardGeneric("gps"))

## The dataset is specified as an RtreemixData object.
setMethod("gps",
          signature(model = "RtreemixModel", data = "RtreemixData"),
          function(model, data, sampling.mode = "exponential", sampling.param = 1,
                    no.sim = 10, seed = (-1)) {

            ## Setting the true type
            sampling.param <- as.numeric(sampling.param)
            no.sim <- as.integer(no.sim)
            seed <- as.integer(seed)           
            ## Check if the necessary parameters are provided and have the correct form.    
            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0 )))
              stop("The seed must be a positive integer!")
            
            if (!missing(sampling.mode)) {
              if (!(identical(sampling.mode, "constant") || identical(sampling.mode, "exponential")))
                stop("Please give correct sampling mode (constant or exponential)!")
            }
            
            if (no.sim < 1)
              stop("The number of iterations for the waiting time simulation must be an integer greater than zero!")            
                   
            mat <- cbind(rep(1, sampleSize(data)), Sample(data))
            ## Check if the patterns have the correct length.            
            if(eventsNum(model) !=  ncol(mat)) {
              cat("The patterns are expected to have length: ", eventsNum(model), ". \n")
              stop("\n")
            }       
                       
            ## Give the node names and edge weights in a list structure. 
            listG <- list()
            for(i in 1:numTrees(model)) {
              listG[[i]] <- list()
              listG[[i]][[1]] <- nodes(Trees(model)[[i]])
              listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
            }
            ## Call to the C function for calculating the GPS.  
            res.GPS <- .Call("R_time", eventsNum(model), listG, Weights(model),
                             as.integer(ifelse(identical(sampling.mode, "constant"), 0, 1)),
                             sampling.param, no.sim, seed)
            ## All possible patterns.
            all.data <- new("RtreemixData", Sample = res.GPS[[1]][, -1])

            resp <- WLikelihoods(likelihoods(model = model, data = all.data))
            bind.gps <- sapply(res.GPS[[2]], rbind)
            bind.gps.bin <- apply(!apply(bind.gps, 2, is.na), 2, as.integer)
            filter.resp <- resp * bind.gps.bin
            resp <- filter.resp / apply(filter.resp, 1, sum)

            if(numTrees(model) > 1) {
              resp[which(apply(bind.gps.bin, 1, sum) == 0),] <- c(1, rep(0, numTrees(model) - 1))
              ## Take the mean for the noise.
              if (sampling.mode == "exponential")  
                res.GPS[[2]][[1]][which(is.na(res.GPS[[2]][[1]]))] <- 1.0/sampling.param
              else 
                res.GPS[[2]][[1]][which(is.na(res.GPS[[2]][[1]]))] <- sampling.param
                      
              for(i in 2:numTrees(model)) 
                res.GPS[[2]][[i]][which(is.na(res.GPS[[2]][[i]]))] <- numeric(1)        
            } else {
              for(i in 1:numTrees(model)) 
                res.GPS[[2]][[i]][which(is.na(res.GPS[[2]][[i]]))] <- numeric(1)
            }
              
            fin.GPS <- apply(t(resp) *  t(sapply(res.GPS[[2]], rbind)), 2, sum)
              
            ## Take the GPS of the patterns present in the given dataset.
            index <- apply(mat[,-1], 1, '%*%', 2^c(0:(ncol(mat) - 2))) + 1              
            res <- fin.GPS[index]
              
            return(new("RtreemixGPS",
                       Data = data,
                       Model = model,
                       SamplingMode = sampling.mode,
                       SamplingParam = sampling.param,
                       GPS = res))
          })

## The dataset is specified as a binary matrix.
setMethod("gps",
          signature(model = "RtreemixModel", data = "matrix"),
          function(model, data, sampling.mode = "exponential", sampling.param = 1,
                    no.sim = 10, seed = (-1)) {

            ## Setting the true type
            sampling.param <- as.numeric(sampling.param)
            no.sim <- as.integer(no.sim)
            seed <- as.integer(seed)           
            ## Check if the necessary parameters are provided and have the correct form.    
            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0)))
              stop("The seed must be a positive integer!")
            
            if (!missing(sampling.mode)) {
              if (!(identical(sampling.mode, "constant") || identical(sampling.mode, "exponential")))
                stop("Please give correct sampling mode (constant or exponential)!")
            }
            
            if (no.sim < 1)
              stop("The number of iterations for the waiting time simulation must be an integer greater than zero!")
            
            mat <- cbind(rep(1, nrow(data)), data)
            ## Check if the patterns have the correct length.            
            if(eventsNum(model) !=  ncol(mat)) {
              cat("The patterns are expected to have length: ", eventsNum(model), ". \n")
              stop("\n")
            }       
                       
            ## Give the node names and edge weights in a list structure. 
            listG <- list()
            for(i in 1:numTrees(model)) {
              listG[[i]] <- list()
              listG[[i]][[1]] <- nodes(Trees(model)[[i]])
              listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
            }
            ## Call to the C function for calculating the GPS.  
            res.GPS <- .Call("R_time", eventsNum(model), listG, Weights(model),
                             as.integer(ifelse(identical(sampling.mode, "constant"), 0, 1)),
                             sampling.param, no.sim, seed)
            ## Generate all possible patterns.
            all.data <- new("RtreemixData", Sample = res.GPS[[1]][, -1])

            resp <- WLikelihoods(likelihoods(model = model, data = all.data))
            bind.gps <- sapply(res.GPS[[2]], rbind)
            bind.gps.bin <- apply(!apply(bind.gps, 2, is.na), 2, as.integer)
            filter.resp <- resp * bind.gps.bin
            resp <- filter.resp / apply(filter.resp, 1, sum)

            if(numTrees(model) > 1) {
              resp[which(apply(bind.gps.bin, 1, sum) == 0),] <- c(1, rep(0, numTrees(model) - 1))
              ## Take the mean for the noise.
              if (sampling.mode == "exponential")  
                res.GPS[[2]][[1]][which(is.na(res.GPS[[2]][[1]]))] <- 1.0/sampling.param
              else 
                res.GPS[[2]][[1]][which(is.na(res.GPS[[2]][[1]]))] <- sampling.param
              
              for(i in 2:numTrees(model)) 
                res.GPS[[2]][[i]][which(is.na(res.GPS[[2]][[i]]))] <- numeric(1)        
            } else {
              for(i in 1:numTrees(model)) 
                res.GPS[[2]][[i]][which(is.na(res.GPS[[2]][[i]]))] <- numeric(1)
            }
            
            fin.GPS <- apply(t(resp) *  t(sapply(res.GPS[[2]], rbind)), 2, sum)
            
            ## Take the GPS of the patterns present in the given dataset.
            index <- apply(mat[,-1], 1, '%*%', 2^c(0:(ncol(mat) - 2))) + 1              
            res <- fin.GPS[index]
            
            return(new("RtreemixGPS",
                       Model = model,
                       Data = new("RtreemixData", Sample = data),
                       SamplingMode = sampling.mode,
                       SamplingParam = sampling.param,
                       GPS = res))                        
          })

## The dataset is not specified. Calculate the GPS for all possible patterns.
setMethod("gps",
          signature(model = "RtreemixModel", data = "missing"),
          function(model, sampling.mode = "exponential", sampling.param = 1,
                    no.sim = 10, seed = (-1)) {
            ## Setting the true type
            sampling.param <- as.numeric(sampling.param)
            no.sim <- as.integer(no.sim)
            seed <- as.integer(seed)           
            ## Check if the necessary parameters are provided and have the correct form.    
            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0 )))
              stop("The seed must be a positive integer!")
            
            if (!missing(sampling.mode)) {
              if (!(identical(sampling.mode, "constant") || identical(sampling.mode, "exponential")))
                stop("Please give correct sampling mode (constant or exponential)!")
            }

            if (no.sim < 1)
              stop("The number of iterations for the waiting time simulation must be an integer greater than zero!")           
                                          
            ## Give the node names and edge weights in a list structure. 
            listG <- list()
            for(i in 1:numTrees(model)) {
              listG[[i]] <- list()
              listG[[i]][[1]] <- nodes(Trees(model)[[i]])
              listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
            }
            ## Call to the C function for calculating the GPS.  
            res.GPS <- .Call("R_time", eventsNum(model), listG, Weights(model),
                             as.integer(ifelse(identical(sampling.mode, "constant"), 0, 1)),
                             sampling.param, no.sim, seed)
            ## Generate all possible patterns.
            all.data <- new("RtreemixData", Sample = res.GPS[[1]][, -1])

            resp <- WLikelihoods(likelihoods(model = model, data = all.data))
            bind.gps <- sapply(res.GPS[[2]], rbind)
            bind.gps.bin <- apply(!apply(bind.gps, 2, is.na), 2, as.integer)
            filter.resp <- resp * bind.gps.bin
            resp <- filter.resp / apply(filter.resp, 1, sum)

            if(numTrees(model) > 1) {
              resp[which(apply(bind.gps.bin, 1, sum) == 0),] <- c(1, rep(0, numTrees(model) - 1))
              ## Take the mean for the noise.
              if (sampling.mode == "exponential")  
                res.GPS[[2]][[1]][which(is.na(res.GPS[[2]][[1]]))] <- 1.0/sampling.param
              else 
                res.GPS[[2]][[1]][which(is.na(res.GPS[[2]][[1]]))] <- sampling.param
                      
              for(i in 2:numTrees(model)) 
                res.GPS[[2]][[i]][which(is.na(res.GPS[[2]][[i]]))] <- numeric(1)        
            } else {
              for(i in 1:numTrees(model)) 
                res.GPS[[2]][[i]][which(is.na(res.GPS[[2]][[i]]))] <- numeric(1)
            }
            ## Compute GPS for all patterns.
            fin.GPS <- apply(t(resp) *  t(sapply(res.GPS[[2]], rbind)), 2, sum)  

            return(new("RtreemixGPS",
                       Model = model,
                       Data = all.data,
                       SamplingMode = sampling.mode,
                       SamplingParam = sampling.param,
                       GPS = fin.GPS))             
          })

###################################################################################
## When the sampling mode and sampling parameter are specified, function for simulating patterns
## and their waiting and sampling times from a given Rtreemix model.
## When only the mixture model is specified, function for drawing data from the mixture model. 
## The model has to be specified in both cases.
if(!isGeneric("sim"))
  setGeneric("sim", function(model, sampling.mode, sampling.param, ...) standardGeneric("sim"))

## Function for drawing patterns from a given Rtreemix model.
setMethod("sim",
          signature(model = "RtreemixModel", sampling.mode = "missing", sampling.param = "missing"),
          function(model, no.draws = 10, seed = (-1)) {
            ## Setting the true type            
            no.draws <- as.integer(no.draws)
            seed <- as.integer(seed)           
            ## Check if the necessary parameters are provided and have the correct form.    
            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0 )))
              stop("The seed must be a positive integer!")
            
            if(no.draws < 1)
              stop("The number of draws must be an integer greater than zero!")
            
            ## Give the node names and edge weights in a list structure.
            listG <- list()
            for(i in 1:numTrees(model)) {
              listG[[i]] <- list()
              listG[[i]][[1]] <- nodes(Trees(model)[[i]])
              listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
            }
            
            ## Call to the C function
            res.sim <- .Call("R_draw", eventsNum(model), listG, Weights(model),
                             no.draws, seed)
      
            return(new("RtreemixData",
                       Sample = res.sim[,-1]))
            
          })

## Function for simulating patterns with their sampling and waiting times.
setMethod("sim",
          signature(model = "RtreemixModel", sampling.mode = "character", sampling.param = "numeric"),
          function(model, sampling.mode, sampling.param, no.sim = 10, seed = (-1)) {
            ## Setting the true type            
            no.sim <- as.integer(no.sim)
            seed <- as.integer(seed)           
            ## Check if the necessary parameters are provided and have the correct form.    
            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0 )))
              stop("The seed must be a positive integer!")

            if (!(identical(sampling.mode, "constant") || identical(sampling.mode, "exponential")))
              stop("Please give correct sampling mode (constant or exponential)!")
            
            if(no.sim < 1)
              stop("The number of iterations for the waiting time simulation must be an integer greater than zero!")            
            
            ## Give the node names and edge weights in a list structure.
            listG <- list()
            for(i in 1:numTrees(model)) {
              listG[[i]] <- list()
              listG[[i]][[1]] <- nodes(Trees(model)[[i]])
              listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
            }
            
            ## Call to the C function
            res.sim <- .Call("R_simulate", eventsNum(model), listG, Weights(model),
                             as.integer(ifelse(identical(sampling.mode, "constant"), 0, 1)),
                             sampling.param, no.sim,
                             seed)
            
            return(new("RtreemixSim",
                       Model = model,
                       SimPatterns = new("RtreemixData", Sample = res.sim$patterns[, -1]),
                       SamplingMode = sampling.mode,
                       SamplingParam = sampling.param,
                       WaitingTimes = res.sim$wtimes,
                       SamplingTimes = res.sim$stimes))    
            
          })

###################################################################################
## Function for generating a random Rtreemix model.
## The number of tree components and the number of genetic events have to be specified.
if(!isGeneric("generate"))
  setGeneric("generate", function(K, no.events, ...) standardGeneric("generate"))

setMethod("generate",
          signature(K = "numeric", no.events = "numeric"),
          function(K, no.events, noise.tree = TRUE,
                   equal.edgeweights = TRUE, prob = c(0.0, 1.0),
                   seed = (-1)) {
            
            ## Setting the true type
            K <- as.integer(K)
            no.events <- as.integer(no.events)
            noise.tree <- as.integer(noise.tree)
            seed <- as.integer(seed)
            equal.edgeweights <- as.integer(equal.edgeweights)
            ## Check if the necessary parameters are provided and have the correct form.                
            if(K < 1)
              stop("The number of trees must be an integer greater than zero!")
            
            if (no.events < 1)
              stop("The number of events must be an integer greater than zero!")
            
            if((seed != -1) &&
               (!(is.numeric(seed)) || !(trunc(seed) == seed) || !(seed >= 0 )))
              stop("The seed must be a positive integer!")
            
            if(!missing(prob)) {
              if(!is.numeric(prob) || (length(prob) != 2))
                stop("Specify the boundaries of the conditional probabilities as a numeric vector
of length two = c(min, max)!")              
              if(prob[2] < prob[1])
                stop("In the probability vector the lower boundary must be smaller than the
upper boundary!")
            }
            ## Call to the C function.  
            res.sim <- .Call("R_random", K, no.events,
                             noise.tree, equal.edgeweights,
                             prob[1], prob[2], seed)

            for(i in 1:K)
              edgeDataDefaults(res.sim$graphs.mixture[[i]], "weight") <- numeric(1)

            res <- new("RtreemixModel",
                       Weights =  as.numeric(res.sim$alpha),
                       Star = as.logical(noise.tree), 
                       Trees = res.sim$graphs.mixture)
  
            return(res)   
          })

###################################################################################
## Function for generating the (scaled) probablility distribution induced with a given RtreemixModel.
## The RtreemixModel has to be specified.
## The sampling mode and the parameters for the sampling times for the observed
## input and output probabilities are optional.
if(!isGeneric("distribution"))
  setGeneric("distribution", function(model, sampling.mode, sampling.param, output.param) standardGeneric("distribution"))

setMethod("distribution",
          signature(model = "RtreemixModel", sampling.mode = "missing",
                    sampling.param = "missing", output.param = "missing"),             
          function(model) {

            if(eventsNum(model) > 30)
              stop("Invalid for more than 30 events!")

            ## There is no specified sampling mode  
            mode = as.integer(-1)

            ## Give the node names and edge weights in a list structure.
            listG <- list()
            for(i in 1:numTrees(model)) {
              listG[[i]] <- list()
              listG[[i]][[1]] <- nodes(Trees(model)[[i]])
              listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
            }
            ## Call to the C function for getting the distribution.
            res.distr <- .Call("R_distr", eventsNum(model), listG, Weights(model),
                               mode, numeric(1), numeric(1))  

            ## All possible patterns
            sample <- as.matrix(res.distr[[1]])
            rownames(sample) <- paste("P", 1:nrow(sample), sep = "")
            ## Build a dataframe of all possible patterns with their corresponding probabilities
            result <- data.frame(sample, res.distr[[2]])
            colnames(result) <- c(paste("E", 0:(ncol(sample) - 1), sep = ""), "probability")
                        
            return(result)  
          })

setMethod("distribution",
          signature(model = "RtreemixModel", sampling.mode = "character",
                    sampling.param = "numeric", output.param = "numeric"),             
          function(model, sampling.mode, sampling.param,
                   output.param) {

            if(eventsNum(model) > 30)
              stop("Invalid for more than 30 events!")
  
            if (!(identical(sampling.mode, "constant") || identical(sampling.mode, "exponential")))
              stop("Please give correct sampling mode (constant or exponential)!")

            mode <- as.integer(ifelse(identical(sampling.mode, "constant"), 0, 1))

            ## Give the node names and edge weights in a list structure.
            listG <- list()
            for(i in 1:numTrees(model)) {
              listG[[i]] <- list()
              listG[[i]][[1]] <- nodes(Trees(model)[[i]])
              listG[[i]][[2]] <- edgeWeights(Trees(model)[[i]])
            }
            ## Call to the C function for getting the distribution.
            res.distr <- .Call("R_distr", eventsNum(model), listG, Weights(model),
                               mode, sampling.param, output.param)
            ## All possible patterns
            sample <- as.matrix(res.distr[[1]])
            rownames(sample) <- paste("P", 1:nrow(sample), sep = "")
            ## Build a dataframe of all possible patterns with their corresponding probabilities
            result <- data.frame(sample, res.distr[[2]])
            colnames(result) <- c(paste("E", 0:(ncol(sample) - 1), sep = ""), "probability")
            ## Print some useful information.
            if (mode == 1) 
              cat("Exponential sampling mode!\n")
            else
              cat("Constant sampling mode!\n")
            
            cat("The parameter for the observed input probabilities:\n")
            print(sampling.param)
            cat("\n The parameter for the observed output probabilities:\n")
            print(output.param)
            print(result)
            
            return(result)  
          })

###################################################################################
## Function for calculating GPS values and their 95% bootstrap 
## confidence intervals for a given dataset.
## the GPS values are calculated with respect to an underlying
## oncogenetic trees mixture model fitted to the given data.
## The number of trees for the model (K) and the dataset have to be specified.

##change to confIntGPS
if(!isGeneric("confIntGPS"))
  setGeneric("confIntGPS", function(data, K, ...) standardGeneric("confIntGPS"))

setMethod("confIntGPS",
          signature(data = "RtreemixData", K = "numeric"),
          function(data, K, sampling.mode = "exponential", sampling.param = 1,
                   no.sim = 10000, B = 1000, equal.star = TRUE) {
            ## Setting the true type
            K <- as.integer(K)
            B <- as.integer(B)
            sampling.param <- as.numeric(sampling.param)            
            no.events <- eventsNum(data)
            no.sim <-  as.integer(no.sim)
            equal.star <- as.logical(equal.star)
            ## Check if the necessary parameters are provided and have the correct form.                
            if(K < 1)
              stop("The number of trees must be an integer greater than zero!")
            
            if(B < 1)
              stop("The number of bootstrap samples integer greater than zero!")

            if (no.sim < 1)
              stop("The number of iterations for the waiting time simulation must be an integer greater than zero!")
            
            if (!missing(sampling.mode)) {
              if (!(identical(sampling.mode, "constant") || identical(sampling.mode, "exponential")))
                stop("Please give correct sampling mode (constant or exponential)!")
            }
  
            ## Fit an Rtreemix fm model to the given data.
            fm <- fit(data = data, K = K, noise = TRUE, equal.edgeweights = equal.star, eps = 0.01)
            ## Compute the GPS values for all possible patterns.
            fit.gps <- GPS(gps(model = fm, no.sim = no.sim))
            ## Find the row indices of the patterns in the given data in the matrix of all possible patterns
            ind <- apply(Sample(data), 1, '%*%', 2^c(0 : (no.events - 2))) + 1
  
            boot.gps <- vector(mode = "list", length = 2 ^ (no.events - 1))  
            for (i in 1:B) {
              ## Draw a bootstrap sample from the data with replacement.
              random.rows <- sample(sampleSize(data), sampleSize(data), replace = TRUE)
              new.data <- new("RtreemixData", Sample = Sample(data)[random.rows, ])
              ## Fit an Rtreemix model to it
              new.fit <- fit(data = new.data, K = K, noise = TRUE, equal.edgeweights = equal.star, eps = 0.01)
              ## Calculate the GPS values of all patterns for the new model
              new.gps <- GPS(gps(model = new.fit, no.sim = no.sim))
              for(j in 1:length(boot.gps))
                boot.gps[[j]] <- c(boot.gps[[j]], new.gps[j])
            }
            sort.boot <- lapply(boot.gps, sort)
            ## Calculate 95% confidence intervals for the GPS values
            boot.intervals <- lapply(sort.boot, function(sort.boot){c(sort.boot[ceiling((B * 2.5) / 100)], sort.boot[as.integer((B * 97.5) / 100)])})
  
            fit.gps.data <- fit.gps[ind]
            result <- matrix(nrow = length(ind), ncol = 2)
  
            for(i in 1:length(ind)) 
              result[i, ] <- boot.intervals[[ind[i]]]
            
            return(new("RtreemixGPS",
                       Model = fm,
                       Data = data,
                       SamplingMode = sampling.mode,
                       SamplingParam = sampling.param,
                       GPS = fit.gps.data,
                       gpsCI = result))
          })
######################################################################          
