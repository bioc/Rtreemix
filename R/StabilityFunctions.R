################################################################################
##  FUNCTIONS FOR STABILITY ANALYSIS OF THE ONCOGENETIC TREES MIXTURE MODEL   ##
################################################################################
## Function for calculating the P value of a given similarity measure.
## What is the probability for obtaining the same or smaller
## value in a random vector of similarity values.
Pval.dist <- function(dist.val, random.vals) {
  return((sum(random.vals <= dist.val) + 1) /(length(random.vals) + 1))
}
######################################################################
## Function for calculating the Kullback-Leibler divergence between
## two discrete probability distributions p and q.
kullback.leibler <- function(p, q) {
  if(length(p) != length(q))
    stop("Error: The distribution vectors have different lengths!")  
  return(sum(p * log(p / q)))
}
######################################################################
## Function for calculating the L1 distance between
## two discrete probability distributions p and q.
L1.dist <- function(p, q) {
  if(length(p) != length(q))
    stop("Error: The distribution vectors have different lengths!")

  return(sum(abs(p - q)))
}
######################################################################
## Function for calculating the cosine distance between
## two discrete probability distributions p and q.
cosin.dist <- function(p, q) {
  if(length(p) != length(q))
    stop("Error: The distribution vectors have different lengths!")

  return(1 - (sum (p * q) / (L2.norm(p) * L2.norm(q)) ))
}
######################################################################
## Function for calculating the L2 norm of a given vector x.
L2.norm <- function(x) {
  return(sqrt(sum(x * x)))
}
######################################################################
## Function for calculating the Euclidian distance between
## two vectors x and y.
euclidian.dist <- function(x, y) {
  if(length(x) != length(y))
    stop("Error: The vectors have different lengths!")

  return(sqrt(drop((x - y) %*% (x - y))))
}
######################################################################
## Function for calculating the rank-correlation distance between
## two vectors x and y.
rank.cor.dist <- function(x, y) {
  if(length(x) != length(y))
    stop("Error: The vectors have different length!")

  return(1 - rcorr(x, y, type = "spearman")[[1]][1, 2])
}
######################################################################
## Function that implements a similarity measure for comparing the topologies
## of the trees of two mixture models mixture1 and mixture2.
## It returns a value that uses the number of different edges for caracterizing
## the difference between the models.
## For the comparisson it is necessary that the two models have the same
## number of tree components build on the same number of genetic events.
## It is assumed that the mixtures have at least two components (the first one is the noise).
comp.models <- function(mixture1, mixture2) {
  if((class(mixture1) != "RtreemixModel") || (class(mixture2) != "RtreemixModel"))
    stop("The specified mixture models should be of type RtreemixModel!")
  
  if(numTrees(mixture1) != numTrees(mixture2))
    stop("The two specified mixture models don't have the same number of trees!")
  
  if(eventsNum(mixture1) != eventsNum(mixture2))
    stop("The two specified mixture models don't have the same number of events!")
  no.trees <- numTrees(mixture1)
  if (no.trees == 1)
    stop("The mixture models should have at least two tree components, where the first one is the noise component!")
  
  no.events <- eventsNum(mixture1)
  
  true.order <- list() ## list specifying the edges in the components of mixture1
  fit.order <- list() ## list specifying the edges in the components of mixture2

  ## Build the vectors that characterize the edges of the components of the mixture models.
  ## The noise components are ignored since they always have the same (star) topology.
  for(k in 2:no.trees) {
    true.vec <- character(0)
    fit.vec <- character(0)
    for(l in 2:no.events) {
      ch <- names(which(sapply(edges(getTree(mixture1, k)), function(x, el) {el %in% x}, nodes(getTree(mixture1, k))[l])))
      true.vec <- c(true.vec, ifelse((identical(ch, character(0)) || is.null(ch)), "e", ch))       
      ch <- names(which(sapply(edges(getTree(mixture2, k)), function(x, el) {el %in% x}, nodes(getTree(mixture2, k))[l])))
      fit.vec <- c(fit.vec, ifelse((identical(ch, character(0)) || is.null(ch)), "e", ch))        
    }
    names(true.vec) <- nodes(getTree(mixture1, k))[-1]
    true.order <- c(true.order, list(true.vec))
    names(fit.vec) <- nodes(getTree(mixture2, k))[-1]
    fit.order <- c(fit.order, list(fit.vec))
  }
  ## Build the comparisson matrix.
  ## The (i, j) element is the number of distinct edges of the
  ## (i + 1)-th and (j + 1)-th component of the models
  ## mixture1 and mixture2 respectively.  
  comp.mat <- matrix(nrow = no.trees - 1, ncol = no.trees - 1)
  for(i in 1:(no.trees - 1))
    for(j in 1:(no.trees - 1)) {
      comp.mat[i, j] <- sum(true.order[[i]] != fit.order[[j]])
    }
  ## Form the similarity pairs and calculate the similarity of the models
  ## as a sum of the number of different edges of the trees in the similarity pairs:
  if(no.trees > 2) {
    rc <- which(comp.mat == min(comp.mat), arr.ind = TRUE)
    diff.sum <- min(comp.mat)
    row.index <- c(rc[1, 1]) ## get the row index of the minimum
    col.index <- c(rc[1, 2]) ## get the column index of the minimum

    for(m in 1:(nrow(comp.mat) - 1)) {
      rc <- which(comp.mat == min(comp.mat[-(row.index), -(col.index)]), arr.ind = TRUE)
      diff.sum <- diff.sum + min(comp.mat[-(row.index), -(col.index)])
      row.index <- c(row.index, rc[1, 1])
      col.index <- c(col.index, rc[1, 2])    
    }
  } else {
    diff.sum <- comp.mat[1, 1]
  }
  ## Return a result between 0 and 1.    
  return(diff.sum/((no.trees - 1) * (no.events - 1)))  

}
######################################################################
## Function that assignes to each node the level at which that node is
## in a specific tree (tree.num) of the trees mixture model.
## The start.val is the maximum depth of the tree with which the tree
## specified with tree.num will be compared.
get.tree.levels <- function(mixture, tree.num, start.val) {
  root <- nodes(getTree(mixture, tree.num))[1]
  levels <- acc(getTree(mixture, tree.num), root)
  ## vec <- rep(max(levels[[1]]), eventsNum(mixture) - 1)
  vec <- rep(start.val, eventsNum(mixture) - 1)
  names(vec) <- nodes(getTree(mixture, tree.num))[-1]
  
  vec[names(levels[[1]])] <- levels[[1]]

  return(vec)
}
######################################################################
## Function that implements a similarity measure for comparing the topologies
## of the trees of two mixture models mixture1 and mixture2.
## It returns a value that uses the number of different edges and the depth at
## which the events are, for caracterizing the difference between the models.
## This measure is more application specific, since the depth at which the
## events are in a tree is very important for disease progression.
## For the comparisson it is necessary that the two models have the same
## number of tree components build on the same number of genetic events.
## It is assumed that the mixtures have at least two components (the first one is the noise).
comp.models.levels <- function(mixture1, mixture2) {
  if((class(mixture1) != "RtreemixModel") || (class(mixture2) != "RtreemixModel"))
    stop("The specified mixture models should be of type RtreemixModel!")
  
  if(numTrees(mixture1) != numTrees(mixture2))
    stop("The two specified mixture models don't have the same number of trees!")

  if(eventsNum(mixture1) != eventsNum(mixture2))
    stop("The two specified mixture models don't have the same number of events!")
  no.trees <- numTrees(mixture1)
  if (no.trees == 1)
    stop("The mixture models should have at least two tree components, where the first one is the noise component!")  
  
  no.events <- eventsNum(mixture1)
  
  true.order <- list() ## list specifying the edges in the components of mixture1
  fit.order <- list() ## list specifying the edges in the components of mixture2

  start.vals1 <- vector(mode = "numeric", length = (no.trees - 1)) ## the depth of the tree components in mixture1 (without the noise component) 
  start.vals2 <- vector(mode = "numeric", length = (no.trees - 1)) ## the depth of the tree components in mixture2 (without the noise component)

  ## Build the vectors that characterize the edges of the components of the mixture models.
  ## The noise components are ignored since they always have the same (star) topology.  
  for(k in 2:no.trees) {
    true.vec <- character(0)
    fit.vec <- character(0)
    for(l in 2:no.events) {
      ch <- names(which(sapply(edges(getTree(mixture1, k)), function(x, el) {el %in% x}, nodes(getTree(mixture1, k))[l])))
      true.vec <- c(true.vec, ifelse((identical(ch, character(0)) || is.null(ch)), "e", ch))       
      ch <- names(which(sapply(edges(getTree(mixture2, k)), function(x, el) {el %in% x}, nodes(getTree(mixture2, k))[l])))
      fit.vec <- c(fit.vec, ifelse((identical(ch, character(0)) || is.null(ch)), "e", ch))        
    }
    names(true.vec) <- nodes(getTree(mixture1, k))[-1]
    true.order <- c(true.order, list(true.vec))
    names(fit.vec) <- nodes(getTree(mixture2, k))[-1]
    fit.order <- c(fit.order, list(fit.vec))

    ## Assign the (maximal depths + 1) for each tree in the mixtures
    start.vals1[k - 1] <- max(acc(getTree(mixture1, k), nodes(getTree(mixture1, k))[1])[[1]] + 1, na.remove = TRUE)
    start.vals2[k - 1] <- max(acc(getTree(mixture2, k), nodes(getTree(mixture2, k))[1])[[1]] + 1, na.remove = TRUE)                          
  }
  ## Build the comparisson matrix.
  ## The (i, j) element is the number of distinct edges of the
  ## (i + 1)-th and (j + 1)-th component of the models
  ## mixture1 and mixture2 respectively.    
  comp.mat <- matrix(nrow = no.trees - 1, ncol = no.trees - 1)
  for(i in 1:(no.trees - 1))
    for(j in 1:(no.trees - 1)) {
      comp.mat[i, j] <- sum(true.order[[i]] != fit.order[[j]])
    }
  ## Form the similarity pairs and calculate the similarity of the models
  ## as a sum of the number of different edges of the trees in the similarity pairs
  ## and their corresponding L1-distance of the levels of the events:
  if(no.trees > 2) {
    rc <- which(comp.mat == min(comp.mat), arr.ind = TRUE)
    diff.sum <- min(comp.mat) +
      L1.dist(get.tree.levels(mixture1, rc[1, 1] + 1, start.vals2[rc[1, 2]]),
              get.tree.levels(mixture2, rc[1, 2] + 1, start.vals1[rc[1, 1]]))
    row.index <- c(rc[1, 1]) ## get the row index of the minimum
    col.index <- c(rc[1, 2]) ## get the column index of the minimum

    for(m in 1:(nrow(comp.mat) - 1)) {
      rc <- which(comp.mat == min(comp.mat[-(row.index), -(col.index)]), arr.ind = TRUE)
      diff.sum <- diff.sum + min(comp.mat[-(row.index), -(col.index)]) +
        L1.dist(get.tree.levels(mixture1, rc[1, 1] + 1, start.vals2[rc[1, 2]]),
                get.tree.levels(mixture2, rc[1, 2] + 1, start.vals1[rc[1, 1]]))
      row.index <- c(row.index, rc[1, 1])
      col.index <- c(col.index, rc[1, 2])   
    }
  } else {
    if(no.trees == 2) {
      diff.sum <- comp.mat[1, 1] +
        L1.dist(get.tree.levels(mixture1, 2, start.vals2[1]),
                get.tree.levels(mixture2, 2, start.vals1[1]))
    } else {
      stop("The specified mixture models must have at least two tree components, where the first one is the noise component!")
    }
  }
    
  return(diff.sum)  

}
######################################################################
## Function that implements a similarity measure for comparing the topologies
## of the nontrivial tree components of a specified mixture model.
## It returns a value that uses the number of different edges for caracterizing
## the difference of the trees in the model.
## For the comparisson it is necessary that the model has at least two
## nontrivial components.
comp.trees <- function(model) {
  no.trees <- numTrees(model)
  if(no.trees <= 2)
    stop("The specified mixture model should have at least two nontrivial tree components.")
  
  no.events <- eventsNum(model)
  true.order <- list() ## list specifying the edges in the nontrivial components of the model  
  ## Build the vectors that characterize the
  ## edges of the nontrivial components of the mixture model.
  for(k in 2:no.trees) {
    true.vec <- character(0)
    for(l in 2:no.events) {
      ch <- names(which(sapply(edges(getTree(model, k)), function(x, el) {el %in% x}, nodes(getTree(model, k))[l])))
      true.vec <- c(true.vec, ifelse((identical(ch, character(0)) || is.null(ch)), "e", ch))             
    }
    names(true.vec) <- nodes(getTree(model, k))[-1]
    true.order <- c(true.order, list(true.vec))    
  }
  ## Calculate the similarity of the nontrivial components of the model 
  ## as a sum of the number of different edges of all combinations of
  ## two different nontrivial trees:
  diff.sum <- 0
  for(k1 in 2:(no.trees - 1)) 
    for(k2 in (k1 + 1):no.trees) 
      diff.sum <- diff.sum + sum(true.order[[k1 - 1]] != true.order[[k2 - 1]])
    
  ## Return a result between 0 and 1.
  return(diff.sum / (choose((no.trees - 1), 2) * (no.events - 1)))
}
######################################################################
## Function that implements a similarity measure for comparing the topologies
## of the nontrivial tree components of a specified mixture model.
## It returns a value that uses the  number of different edges and the depth at
## which the events are, for caracterizing the difference of the trees in the model.
## For the comparisson it is necessary that the model has at least two
## nontrivial components.
comp.trees.levels <- function(model) {
  no.trees <- numTrees(model)
  if(no.trees <= 2)
    stop("The specified mixture model should have at least two nontrivial tree components.")
  
  no.events <- eventsNum(model)
  start.vals <- vector(mode = "numeric", length = (no.trees - 1)) ## the depth of the nontrivial tree components in the model 
  true.order <- list() ## list specifying the edges in the nontrivial components of the model 
  
  ## Build the vectors that characterize the
  ## edges of the nontrivial components of the mixture model.  
  for(k in 2:no.trees) {
    true.vec <- character(0)
    for(l in 2:no.events) {
      ch <- names(which(sapply(edges(getTree(model, k)), function(x, el) {el %in% x}, nodes(getTree(model, k))[l])))
      true.vec <- c(true.vec, ifelse((identical(ch, character(0)) || is.null(ch)), "e", ch))             
    }
    names(true.vec) <- nodes(getTree(model, k))[-1]
    true.order <- c(true.order, list(true.vec))
    
    ## Assign the maximal depths + 1 for each tree in the mixtures
    start.vals[k - 1] <- max(acc(getTree(model, k), nodes(getTree(model, k))[1])[[1]] + 1, na.remove = TRUE)    
  }
  ## Calculate the similarity of the models as a sum of the number of
  ## different edges of all possible combinations of two different
  ## nontrivial trees in the model and their corresponding L1-distance
  ## of the levels of the events:
  diff.sum <- 0
  for(k1 in 2:(no.trees - 1)) {
    for(k2 in (k1 + 1):no.trees) {
      diff.sum <- diff.sum + (sum(true.order[[k1 - 1]] != true.order[[k2 - 1]])
                              +  L1.dist(get.tree.levels(model, k1, start.vals[k2 - 1]),
                                         get.tree.levels(model, k2, start.vals[k1 - 1])))
    }
  }
  return(diff.sum)
}
######################################################################
## Stability analysis of the oncogenetic trees mixture model.
## This includes stability analysis on different levels of the mixture
## model: GPS values, encoded probability distribution, tree topologies.
## The function outputs a list of analysis for the mentioned features.
## Each analysis contains the values of different similarity measures
## with their corresponding p-values.
## The models should have at least two components (the first one is the noise).
stability.sim <- function(no.trees = 3, ## number of tree components
                          no.events = 9, ## number of genetic events
                          prob = c(0.2, 0.8), ## interval for the edgeweights of the random mixture models 
                          no.draws = 300, ## number of samples drawn from the random models used for learning a mixture model
                          no.rands = 100, ## number of rands for calculatin the p-values
                          no.sim = 1 ## number of simulation iterations
                          ) {
  ## Set the true types.
  no.trees <- as.integer(no.trees)
  no.events <- as.integer(no.events)
  no.draws <- as.integer(no.draws)
  no.rands <- as.integer(no.rands)
  no.sim <- as.integer(no.sim)
  ## Check if the necessary parameters are provided and have the correct form.
  if(no.trees < 2)
    stop("The specified mixture model should have at least two
tree components (where the first one is the noise).")
  if (no.events < 1)
    stop("The number of events must be an integer greater than zero!")
  if(no.draws < 1)
    stop("The number of draws (number of samples) must be an integer greater than zero!")  
  if (no.rands < 1)
    stop("The number of random models for the p-values  must be an integer greater than zero!")   
  if (no.sim < 1)
    stop("The number of iterations for the waiting time simulation
must be an integer greater than zero!")
  if(!missing(prob)) {
    if(!is.numeric(prob) || (length(prob) != 2))
      stop("Specify the boundaries of the conditional probabilities as a numeric vector
of length two = c(min, max)!")  
    if(prob[2] < prob[1])
      stop("In the probability vector the lower boundary must be smaller than the
upper boundary!")
  }
  
  simGPS <- numeric(0) ## results from the stability analysis of the GPS values
  topo.dif <- numeric(0) ## results from the stability analysis of the topologies of the tree components of mixture models based on comp.models
  topo.levels.dif <- numeric(0) ## results from the stability analysis of the topologies of the tree components of mixture models based on comp.models.levels
  result.distr <- numeric(0) ## results from the stability analysis of the probability distribution encoded by the model
  
  mat.true.gps <- numeric(0) ## matrix containing the true GPS values from each simulation iteration
  mat.fit.gps <- numeric(0) ## matrix containing the corresponding fitted GPS values from each simulation iteration

  list.true.models <- list() ## list containing the true models from each simulation iteration
  list.fit.models <- list() ## list containing the corresponding fitted models from each simulation iteration
  ## Simulation iterations:
  for(i in 1:as.integer(no.sim)) {
    print(i) ## simulation iteration
    ## Pick a true model from the space of random models and draw data from it.
    true.m <- generate(K = no.trees, no.events = no.events, prob = prob)
    Weights(true.m) <- c(0.05, rep(0.95/(no.trees - 1), (no.trees - 1)))
    sim.data <- sim(model = true.m, no.draws = no.draws)    
    list.true.models <- c(list.true.models, true.m)
    ## Calculate the GPS and the distribution encoded by true.m
    true.gps <- GPS(gps(model = true.m, data = sim.data, no.sim = 10000))
    mat.true.gps <- cbind(mat.true.gps, true.gps)
    true.distr <- distribution(model = true.m)$probability
    ## Fit a mixture model to the data sim.data    
    fm <- fit(data = sim.data, K = no.trees, noise = TRUE, equal.edgeweights = TRUE, eps = 0.01)
    list.fit.models <- c(list.fit.models, fm)
    ## Calculate the GPS and the distribution encoded by fm
    fit.gps <- GPS(gps(model = fm, data = sim.data, no.sim = 10000))
    mat.fit.gps <- cbind(mat.fit.gps, fit.gps)
    fit.distr <- distribution(model = fm)$probability
    ## Compute different distances between the distributions induced
    ## by the true and fitted models:
    true.cosin <- cosin.dist(true.distr, fit.distr)
    true.l1 <- L1.dist(true.distr, fit.distr)
    true.kull <- kullback.leibler(true.distr, fit.distr)
    ## Compute different distances between the GPS values
    ## for the true and fitted models:
    true.euclGPS <- euclidian.dist(true.gps, fit.gps)
    true.rank.dist <-  rank.cor.dist(true.gps, fit.gps)
    ## Compute different similarity measures for comparing the
    ## topologies of the nontrivial components of the true and fitted models:
    true.topo.dif <- comp.models(true.m, fm)
    true.topo.levels.dif <- comp.models.levels(true.m, fm)
    ## Vectors for keeping the values for the different features
    ## of the random mixture models needed for the p-values calculation:
    rand.euclGPS <- numeric(0)
    rand.rankGPS <- numeric(0)

    rand.cos.distr <- numeric(0)
    rand.l1.distr <- numeric(0)
    rand.kull.distr <- numeric(0)

    rand.topo <- numeric(0)
    rand.topo.levels <- numeric(0)
    ## Create the vectors with random values for the p-values calculation:
    for(j in 1:(as.integer(no.rands) - 1)) {
      ## Generate random model and calculate the GPS and distribution:
      model <- generate(K = no.trees, no.events = no.events, prob = prob)
      Weights(model) <- c(0.05, rep(0.95/(no.trees - 1), (no.trees - 1)))      
      random.gps <- GPS(gps(model = model, data = sim.data, no.sim = 10000))
      random.distr <- distribution(model = model)$probability
      ## GPS:
      rand.euclGPS <- c(rand.euclGPS, euclidian.dist(true.gps, random.gps))      
      rand.rankGPS <- c(rand.rankGPS, rank.cor.dist(true.gps, random.gps))
      ## Distribution:
      rand.cos.distr <- c(rand.cos.distr, cosin.dist(true.distr, random.distr))
      rand.l1.distr <- c(rand.l1.distr, L1.dist(true.distr, random.distr))      
      rand.kull.distr <- c(rand.kull.distr, kullback.leibler(true.distr, random.distr))
      ## Tree topologies:
      rand.topo <- c(rand.topo, comp.models(true.m, model))
      rand.topo.levels <- c(rand.topo.levels, comp.models.levels(true.m, model))
    }
    ## Stability analysis of GPS values
    ## (using Euclidian distance and rank correlation distance as similarity measures):
    simGPS <- rbind(simGPS,
                    c(true.euclGPS, Pval.dist(true.euclGPS, rand.euclGPS),   
                      true.rank.dist, Pval.dist(true.rank.dist, rand.rankGPS)))
    ## Stability analysis of induced distributions
    ## (using cosine distance, L1 distance, Kullback-Leibler divergence as similarity measures):
    result.distr <- rbind(result.distr,
                          c(true.cosin, Pval.dist(true.cosin, rand.cos.distr),
                            true.l1, Pval.dist(true.l1, rand.l1.distr),
                            true.kull, Pval.dist(true.kull, rand.kull.distr)))
    ## Stability analysis of tree topologies
    ## (using comp.models and comp.models.levels as similarity measures):
    topo.dif <- rbind(topo.dif,
                      c(true.topo.dif, Pval.dist(true.topo.dif, rand.topo)))
    topo.levels.dif <- rbind(topo.levels.dif,
                             c(true.topo.levels.dif, Pval.dist(true.topo.levels.dif, rand.topo.levels)))
    
  }
  ## Organize the output:
  colnames(simGPS) <- c("eucl.dist", "p.val.eucl", "rank.cor.dist", "p.val.rank")
  colnames(result.distr) <- c("cos", "p.val.cos", "L1", "p.val.L1", "KL", "p.val.KL")
  colnames(topo.dif) <- c("topo.dif", "p.value")
  colnames(topo.levels.dif) <- c("topo.levels.dif", "p.value")
  colnames(mat.true.gps) <- c(1:no.sim)
  colnames(mat.fit.gps) <- c(1:no.sim)
  output <- list(simGPS, result.distr,
                 topo.dif, topo.levels.dif,
                 mat.true.gps, mat.fit.gps,
                 list.true.models, list.fit.models)
  names(output) <- c("GPS", "Distribution", "Tree topologies (distinct edges)",
                     "Tree topologies (distinct edges + levelsL1dist)", "True GPS",
                     "Fitted GPS", "True models", "Fitted models")
    
  return(output)  
}
