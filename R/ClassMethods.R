##################################################################################
###                             CLASS METHODS                                 ###
##################################################################################

## Constructor for creating a new object of the class RtreemixData.
## The data can also be read from a file with a specific format.
setMethod("initialize", "RtreemixData",
          function(.Object, Sample, Patients, Events, Description, File) {
            ## Read from file
            if(!missing(File)) {
              file.in <- file(File, "r", blocking = FALSE)
              while(length(.Object@Description <- scan(file.in, character(0), sep = " ", nlines = 1,
                                           quiet = TRUE, comment.char = "#")) == 0) {}

             .Object@Description <- paste(.Object@Description, collapse = " ")
              
              mat.sample <-  NULL
              
              while(length(line <- scan(file.in, character(0), sep = " ", nlines = 1,
                           quiet = TRUE, comment.char = "#")) == 0) {}
              
              while(length(line) > 0) {
                mat.sample <- rbind(mat.sample, as.integer(line))
                line <- scan(file.in, character(0), sep = " ", nlines = 1,
                             quiet = TRUE, comment.char = "#")
              }
 
              .Object@Sample <- mat.sample
              
              while(length(.Object@Patients <- scan(file.in, character(0), sep = " ", nlines = 1,
                                           quiet = TRUE, comment.char = "#")) == 0) {}              
              
              while(length(.Object@Events <- scan(file.in, character(0), sep = " ", nlines = 1,
                                          quiet = TRUE, comment.char = "#")) == 0) {}
              .Object@Events <- c("0", .Object@Events)
              close(file.in)

              if(.Object@Patients[1] == "NULL")
                .Object@Patients <- paste("P", 1:nrow(.Object@Sample), sep = "")
              else 
                if((length(.Object@Patients) != nrow(.Object@Sample)))
                  stop("From file: Please give the corect number of patient IDs!")
              
              if(.Object@Events[1] == "NULL")
                .Object@Events <- paste("E", 0:ncol(.Object@Sample), sep = "")
              else
                if((length(.Object@Events) != ncol(.Object@Sample) + 1))
                  stop("From file: Please give the corect number of events")               

              return(.Object)

            }
            ## Use given data
            if(missing(Sample))
              stop("The sample matrix is missing!")
             
            if(!is.matrix(Sample))
              stop("The provided Sample must be a binary matrix.")

            .Object@Sample <- Sample

            if(missing(Patients) && (nrow(Sample) != 0))
              Patients <- paste("P", 1:nrow(Sample), sep = "")

            if(missing(Patients) && (nrow(Sample) == 0))
              Patients <- character(0)            

            if((length(Patients) != nrow(Sample)))
              stop("Please give the corect number of patient IDs!")

            .Object@Patients <- Patients

            if(missing(Events))
              Events <- paste("E", 0:ncol(Sample), sep = "")

            if((length(Events) != (ncol(Sample) + 1)))
              stop("Please give the correct number of events!")

            .Object@Events <- Events

            .Object@Description <- ifelse(missing(Description), character(0), Description)
            
            .Object
          })

## Accessor functions for the slots of the class RtreemixData.
if(!isGeneric("Sample"))
  setGeneric("Sample", function(object) standardGeneric("Sample"))

setMethod("Sample", "RtreemixData", function(object) object@Sample)

if(!isGeneric("Patients"))
  setGeneric("Patients", function(object) standardGeneric("Patients"))

setMethod("Patients", "RtreemixData", function(object) object@Patients)

if(!isGeneric("Events"))
  setGeneric("Events", function(object) standardGeneric("Events"))

setMethod("Events", "RtreemixData", function(object) object@Events)

if(!isGeneric("Description"))
  setGeneric("Description", function(object) standardGeneric("Description"))

setMethod("Description", "RtreemixData", function(object) object@Description)

## Replacement functions for the slots of the class RtreemixData.

## We don't allow replacement of the sample.
##if(!isGeneric("Sample<-"))
##  setGeneric("Sample<-", function(object, value) standardGeneric("Sample<-"))

##setReplaceMethod("Sample", "RtreemixData", function(object, value) {object@Sample <- value; object})

if(!isGeneric("Patients<-"))
  setGeneric("Patients<-", function(object, value) standardGeneric("Patients<-"))

setReplaceMethod("Patients", "RtreemixData", function(object, value) {
  if(length(value) != nrow(object@Sample))
    stop("Please give the corect number of patient IDs!")
  else {
    object@Patients <- value;
    object
  }
})

if(!isGeneric("Events<-"))
  setGeneric("Events<-", function(object, value) standardGeneric("Events<-"))

setReplaceMethod("Events", "RtreemixData", function(object, value) {
  if((length(value) != (ncol(object@Sample) + 1)) && (length(object@Sample) != 0))
    stop("Please give the correct number of events!")
  else {
    object@Events <- value;
    object
  }
})

if(!isGeneric("Description<-"))
  setGeneric("Description<-", function(object, value) standardGeneric("Description<-"))

setReplaceMethod("Description", "RtreemixData", function(object, value) {
  if (!is.character(value))
    stop("The specified description must be of type character!")
  else {
    object@Description <- value;
    object
  }
})

## A method for getting the number of patients i.e. the sample size N.
if(!isGeneric("sampleSize"))
  setGeneric("sampleSize", function(object) standardGeneric("sampleSize"))

setMethod("sampleSize", signature = signature(object = "RtreemixData"),
          definition = function(object) 
          return(nrow(object@Sample))
          )

## A method for getting the number of genetic events L.
if(!isGeneric("eventsNum"))
  setGeneric("eventsNum", function(object) standardGeneric("eventsNum"))

setMethod("eventsNum", signature = signature(object = "RtreemixData"),
          definition = function(object) 
          return(length(object@Events))      
          )

## Printing function for the class RtreemixData.
setMethod("print", "RtreemixData", function(x, ...) .printRtreemixData(x))
setMethod("show", "RtreemixData", function(object) .printRtreemixData(x = object))

.printRtreemixData <- function(x) {
  
  if(length(Sample(x)) == 0)               
    return(cat("Empty object!\n"))
  
  cat("RtreemixData object with", eventsNum(x), "genetic events from",
      sampleSize(x), "patients.\n\n")
  if(!identical(Description(x), character(0)))
    cat("Description: ", Description(x), "\n\n")
  
  cat("Patients: ", Patients(x), "\n\n")
  cat("Events: ", Events(x), "\n\n")
  cat("First 10 samples:\n\n") 
  a <- cbind(rep(1, sampleSize(x)), Sample(x))
  rownames(a) <- Patients(x)
  colnames(a) <- Events(x)
  if(nrow(a) > 10) {
    a <- a[1:10, ]
  }
  print(a)
}


################################################################################

## Constructor for creating a new object of the class RtreemixModel.
setMethod("initialize", "RtreemixModel",
          definition = function(.Object, ParentData, Weights, WeightsCI, Resp,
            CompleteMat, Star, Trees) {
            ## Initializing the slots inherited from the parent class
            if(!missing(ParentData)) {
              .Object@Sample <- Sample(ParentData)
              .Object@Patients <- Patients(ParentData)
              .Object@Events <- Events(ParentData)
              .Object@Description <- Description(ParentData)
            }
            ## Initializing the slots from RtreemixModel
            if(missing(Weights))
              .Object@Weights <- numeric(0)
            else {
              .Object@Weights <- Weights              
              if (length(.Object@Weights) != 0)
                names(.Object@Weights) <- paste("Alpha", 1:length(.Object@Weights), sep = "")
            }

            if(missing(WeightsCI))
              .Object@WeightsCI <- list()
            else {
              if(missing(Weights))
                stop("Confidence intervals for the weights cannot be specified when the weights are missing!")
              
              .Object@WeightsCI <- WeightsCI
              if (length(.Object@WeightsCI) != 0) {
                if(length(Weights) != length(WeightsCI))
                  stop("The length of the list specifying the confidence intervals
for the mixture weights have to be equal to their length.")
                
                names(.Object@WeightsCI) <- paste("Alpha", 1:length(.Object@WeightsCI), sep = "")
                
                for(i in 1:length(.Object@WeightsCI)) {
                  if(!(length(.Object@WeightsCI[[i]]) == 2))
                    stop("The confidence intervals for the weights are not correctly specified!")
                else
                  names(.Object@WeightsCI[[i]]) <- c("lower", "upper")
                }
              }
            }
            
            if(missing(Resp)) {
              .Object@Resp <- matrix(numeric(0), 0, 0)
            } else {
              if((nrow(Resp) != length(Trees)) || (ncol(Resp) != nrow(.Object@Sample)))
                stop("Uncorrect assignment of the responsibilities!")
              .Object@Resp <- Resp
              if(nrow(.Object@Resp) != 0)
                rownames(.Object@Resp) <- paste("T", 1:nrow(Resp), sep = "")
            }
            
            if(missing(CompleteMat))
               .Object@CompleteMat <- matrix(integer(0), 0, 0)
            else
              .Object@CompleteMat <- CompleteMat

            if(missing(Star))
              .Object@Star <-  logical(0)
            else
              .Object@Star <- Star

            if(missing(Trees) || (length(Trees) == 0))
              .Object@Trees <- list()
            else {
              if(sum(sapply(Trees, function(x) {class(x) == "graphNEL"})) != length(Trees))
                stop("The tree components are not objects of class graphNEL!")
              if(missing(ParentData)) {
                 .Object@Events <-
                   paste("E", 0:(max(sapply(Trees,
                                            function(x) {return(length(nodes(x)))})) - 1), sep = "")
              }
              .Object@Trees <- Trees
              if(!identical(length(.Object@Trees), as.integer(0)))
                names(.Object@Trees) <-  paste("T", 1:length(.Object@Trees), sep = "")
            }
            
            if(!identical(length(.Object@Trees), length(.Object@Weights)))
              stop("The number of mixture parameters doesn't match the number of tree components!")
            
            .Object
            
          })

## Accessor functions for the class RtreemixModel.
if(!isGeneric("Weights"))
  setGeneric("Weights", function(object) standardGeneric("Weights"))

setMethod("Weights", "RtreemixModel", function(object) object@Weights)

if(!isGeneric("WeightsCI"))
  setGeneric("WeightsCI", function(object) standardGeneric("WeightsCI"))

setMethod("WeightsCI", "RtreemixModel", function(object) object@WeightsCI)

if(!isGeneric("Resp"))
  setGeneric("Resp", function(object) standardGeneric("Resp"))

setMethod("Resp", "RtreemixModel", function(object) object@Resp)

if(!isGeneric("CompleteMat"))
  setGeneric("CompleteMat", function(object) standardGeneric("CompleteMat"))

setMethod("CompleteMat", "RtreemixModel", function(object) object@CompleteMat)

if(!isGeneric("Star"))
  setGeneric("Star", function(object) standardGeneric("Star"))

setMethod("Star", "RtreemixModel", function(object) object@Star)

if(!isGeneric("Trees"))
  setGeneric("Trees", function(object) standardGeneric("Trees"))

setMethod("Trees", "RtreemixModel", function(object) object@Trees)

## Replacement functions for the class RtreemixModel.
if(!isGeneric("Weights<-"))
  setGeneric("Weights<-", function(object, value) standardGeneric("Weights<-"))
## The weights are a numeric vector with components that add up to 1.
setReplaceMethod("Weights", "RtreemixModel",
                 function(object, value) {
                   if((length(value) != numTrees(object)) || !is.numeric(value) || (sum(value) != 1))
                     stop("Uncorrect assignment of the weights of the tree components!")
                   else {
                     object@Weights <- value;
                     object
                   }
                 })

## Replacement functions for the other slots are not allowed.
## if(!isGeneric("WeightsCI<-"))
##   setGeneric("WeightsCI<-", function(object, value) standardGeneric("WeightsCI<-"))

## setReplaceMethod("WeightsCI", "RtreemixModel",
##                  function(object, value) {object@WeightsCI <- value; object})

## if(!isGeneric("Resp<-"))
##   setGeneric("Resp<-", function(object, value) standardGeneric("Resp<-"))

## setReplaceMethod("Resp", "RtreemixModel",
##                  function(object, value) {object@Resp <- value; object})

## if(!isGeneric("CompleteMat<-"))
##   setGeneric("CompleteMat<-", function(object, value) standardGeneric("CompleteMat<-"))

## setReplaceMethod("CompleteMat", "RtreemixModel",
##                  function(object, value) {object@CompleteMat <- value; object})

## if(!isGeneric("Star<-"))
##   setGeneric("Star<-", function(object, value) standardGeneric("Star<-"))

## setReplaceMethod("Star", "RtreemixModel",
##                  function(object, value) {object@Star <- value; object})

## if(!isGeneric("Trees<-"))
##   setGeneric("Trees<-", function(object, value) standardGeneric("Trees<-"))

## setReplaceMethod("Trees", "RtreemixModel",
##                  function(object, value) {object@Trees <- value; object})


## A method for getting the number of trees in the mixture model.
if(!isGeneric("numTrees"))
  setGeneric("numTrees", function(object) standardGeneric("numTrees"))

setMethod("numTrees", signature = signature(object = "RtreemixModel"),
          definition = function(object) 
          return(as.integer(length(object@Trees)))      
          )

## A method for getting a specific tree k from the list of trees in the mixture model.
if(!isGeneric("getTree"))
  setGeneric("getTree", function(object, k) standardGeneric("getTree"))

setMethod("getTree", signature = signature(object = "RtreemixModel", k = "numeric"),
          definition = function(object, k) {
            if(k > numTrees(object))
              stop("The specified k is larger than the number of trees the model contains!")
            if(!(is.numeric(k)) || !(trunc(k) == k) || !(k > 0))
              stop("The tree number must be an integer greater than zero!")
            
            return(object@Trees[[k]])
          }
         )

## A method for getting the data.
if(!isGeneric("getData"))
  setGeneric("getData", function(object) standardGeneric("getData"))
  
setMethod("getData", signature = signature(object = "RtreemixModel"),
          definition = function(object) {
            if (length(Sample(object)) != 0)
              return(new("RtreemixData",
                         Sample = Sample(object),
                         Patients = Patients(object),
                         Events = Events(object)))
            else
              print("Empty object! The random mixture models are not estimated from a given dataset!")
          })

## Printing function for the class RtreemixModel.
setMethod("print", "RtreemixModel", function(x, ...) .printRtreemixModel(x))
setMethod("show", "RtreemixModel", function(object) .printRtreemixModel(x = object))

.printRtreemixModel <- function(x) {
  
  if(numTrees(x) == 0)               
    return(cat("Empty object!\n"))
  
  cat("RtreemixModel object with", numTrees(x), "tree components.\n\n")
  cat("The weights of the trees are:\n")
  print(Weights(x))
  
  cat("\n The trees are:\n")
  print(Trees(x))
}

## Helper function for plotting. Taken from the package Rgraphviz with some changes
## for regulating the size of the text labels in the graphs.
.plott1 <- function (x, y, ..., main = NULL, cex.main = NULL, 
                  col.main = "black", sub = NULL, cex.sub = NULL, col.sub = "black", 
                  drawNode = drawAgNode, xlab, ylab) {

        plot.new()
        old.mai = par(mai = 0.01 + c(0.83 * (!is.null(sub)), 
            0, 0.83 * (!is.null(main)), 0))
        on.exit(par(mai = old.mai), add = TRUE)
        ur <- upRight(boundBox(x))
        plot.window(xlim = c(0, getX(ur)), ylim = c(0, getY(ur)), 
            log = "", asp = NA, ...)
            box()
        xy <- xy.coords(NA, NA)
        def.nodeFontSize <- as.integer(getDefaultAttrs()$node$fontsize)
        def.edgeFontSize <- as.integer(getDefaultAttrs()$edge$fontsize)
        
        plot.xy(xy, type = "n", ...)
        if (!missing(xlab) && !missing(ylab)) 
            stop("Arguments 'xlab' and 'ylab' are not handled.")
        if (!is.null(sub) || !is.null(main)) 
            title(main, sub, cex.main = cex.main, col.main = col.main, 
                cex.sub = cex.sub, col.sub = col.sub)
        agn <- AgNode(x)
        nodeDims <- sapply(agn, function(n) {
            c(getNodeRW(n) + getNodeLW(n), getNodeHeight(n))
        })
        strDims <- sapply(agn, function(n) {
            s <- labelText(txtLabel(n))
            
            if (length(s) == 0) {
                rv <- c(strwidth(" "), strheight(" "))
            }
            else {
                rv <- c(strwidth(s) * (1.1 * def.nodeFontSize / labelFontsize(txtLabel(n))),
                        strheight(s) * (1.4 * def.nodeFontSize / labelFontsize(txtLabel(n))))
            }
            return(rv)
        })
        cex <- min(nodeDims/strDims)
        if (is.finite(cex)) {
            old.cex <- par(cex = cex)
            on.exit(par(cex = old.cex), add = TRUE)
        }
        if (length(drawNode) == 1) {
            lapply(agn, drawNode)
        }
        else {
            if (length(drawNode) == length(AgNode(x))) {
                for (i in seq(along = drawNode)) {
                  drawNode[[i]](agn[[i]])
                }
            }
            else {
                stop(paste("Length of the drawNode parameter is ", 
                  length(drawNode), ", it must be either length 1 or the number of nodes.", 
                  sep = ""))
            }
        }
        arrowLen <- par("pin")[1]/diff(par("usr")[1:2]) * min(nodeDims)/pi
        lapply(AgEdge(x), lines, len = arrowLen, edgemode = edgemode)
        invisible(x)
    }

.plott2 <- function(x, fontSize = 10){       
            
  ## we set the global Graphviz attributes
  graphAttrs <- getDefaultAttrs(layoutType = 'dot')
  graphAttrs <- getDefaultAttrs()
  graphAttrs$cluster <- NULL

  #graphAttrs$graph$splines <- FALSE
  
  ## set the fontsize for the node labels
  graphAttrs$node$fontsize <- fontSize
  graphAttrs$node$fillcolor <- "lightblue"
  #graphAttrs$node$height <- 0.5
  #graphAttrs$node$width <- 0.5

  ## set the local attributes lists
  nodeAttrs <- list()
  edgeAttrs <- list()
  
  
  ## node labels
  node.names <- nodes(x)
  nodeAttrs$label <- nodes(x)
  #nodeAttrs$label <- paste("  ", nodes(x), "  ", sep = "")
  names(nodeAttrs$label) <- node.names

  edgeW <- format(unlist(edgeWeights(x)), digit = 1)
  names(edgeW) <- edgeNames(x)
  toRemove <- removedEdges(x)
  
  if (length(toRemove) > 0)
    edgeW <- edgeW[-toRemove]
    
  edgeAttrs$label <- edgeW

  return(agopen(graph = x, name = "some name", attrs = graphAttrs,
                nodeAttrs = nodeAttrs,  edgeAttrs = edgeAttrs))
}           
   
## Function for plotting a graph x using the two helper functions defined above.
.plotRes <- function(x, fontSize) {
 .plott1(.plott2(x, fontSize))
}

## The plot function for RtreemixModel.
if(!isGeneric("plot"))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod("plot",
          signature(x = "RtreemixModel", y = "missing"),
          function(x, k, fontSize = 8) {  
  
            require('Rgraphviz') || stop('package Rgraphviz is required')     
            if(missing(k)) {
             if(numTrees(x) > 1) {
                par(mfrow = c(ceiling(numTrees(x)/2), 2))
                lapply(Trees(x), .plotRes, fontSize)
             } else
              .plotRes(getTree(x, 1), fontSize)
            } else  
            if ((k > numTrees(x)) || (k < 1)) {
                stop("The tree number must be an integer greater than zero and smaller or equal to the number of tree components in the mixture model!")                
            } else {
                .plotRes(getTree(x, k), fontSize)
            }
  
         })  


## setMethod("plot", "RtreemixModel",
##           definition = function(x, y, ...){      
##          for(i in (1:numTrees(x))) {
##                   ew <- edgeWeights(Trees(x)[[i]])
##              lw <- unlist(unlist(ew))
##              toRemove <- removedEdges(Trees(x)[[i]])
##              if (length(toRemove) > 0) lw <- lw[-toRemove]
##              names(lw) <- edgeNames(Trees(x)[[i]])
##                   eAttrs <- list()
##              eAttrs$label <- lw
##          X11()
##          plot(Trees(x)[[i]], edgeAttrs = eAttrs)
##      } 
##           })


## Constructor for creating a new object of the class RtreemixSim.
setMethod("initialize", "RtreemixSim",
          definition = function(.Object, Model, SimPatterns,
            SamplingMode, SamplingParam, WaitingTimes, SamplingTimes) {
            if(missing(Model) || (class(Model) != "RtreemixModel"))
              stop("Please specify the RtreemixModel from which the simulations are drawn!")
            else {
              ## Initializing the slots inherited from the parent classes (RtreemixData and RtreemixModel) 
              .Object@Sample <- Sample(Model)
              .Object@Patients <- Patients(Model)
              .Object@Events <- Events(Model)
              .Object@Description <- Description(Model)
              .Object@Weights <- Weights(Model)
              .Object@WeightsCI <- WeightsCI(Model)
              .Object@Resp <- Resp(Model)
              .Object@CompleteMat <- CompleteMat(Model)
              .Object@Trees <- Trees(Model)
            } 
            ## Initializing the slots from RtreemixSim
            if(class(SimPatterns) != "RtreemixData")
              stop("The simulated patterns have to be given as an RtreemixData object!")
            
            if(missing(SimPatterns))
              stop("Please specify the simulated patterns!")
            else
              .Object@SimPatterns <- SimPatterns

            if(missing(SamplingMode))
              .Object@SamplingMode <- character(0)
            else {
              if ((SamplingMode != "exponential") && (SamplingMode != "constant"))
                stop("The sampling mode has to be exponential or constant!")
              
              .Object@SamplingMode <- SamplingMode
            }

            if(missing(SamplingParam))
              .Object@SamplingParam <- numeric(0)
            else {        
              .Object@SamplingParam <- SamplingParam
            }            
            
            if(missing(WaitingTimes))
              .Object@WaitingTimes <- numeric(0)
            else {
              if (length(WaitingTimes) != sampleSize(SimPatterns))
                stop("The waiting times are missing for some patterns!")
              
              .Object@WaitingTimes <- WaitingTimes
            }
            
            if(missing(SamplingTimes))
              .Object@SamplingTimes <- numeric(0)
            else {
              if (length(SamplingTimes) != sampleSize(SimPatterns))
                stop("The sampling times are missing for some patterns!")
              
              .Object@SamplingTimes <- SamplingTimes
            }
            
            .Object
          })

## Accessor functions for the slots of the class RtreemixSim.
if(!isGeneric("SimPatterns"))
  setGeneric("SimPatterns", function(object) standardGeneric("SimPatterns"))

setMethod("SimPatterns", "RtreemixSim", function(object) object@SimPatterns)

if(!isGeneric("SamplingMode"))
  setGeneric("SamplingMode", function(object) standardGeneric("SamplingMode"))

setMethod("SamplingMode", "RtreemixSim", function(object) object@SamplingMode)

if(!isGeneric("SamplingParam"))
  setGeneric("SamplingParam", function(object) standardGeneric("SamplingParam"))

setMethod("SamplingParam", "RtreemixSim", function(object) object@SamplingParam)

if(!isGeneric("WaitingTimes"))
  setGeneric("WaitingTimes", function(object) standardGeneric("WaitingTimes"))

setMethod("WaitingTimes", "RtreemixSim", function(object) object@WaitingTimes)

if(!isGeneric("SamplingTimes"))
  setGeneric("SamplingTimes", function(object) standardGeneric("SamplingTimes"))

setMethod("SamplingTimes", "RtreemixSim", function(object) object@SamplingTimes)          

## Replacement functions for the slots of the class RtreemixSim.
## Replacements of the slots is not a good idea. 
## if(!isGeneric("SimPatterns<-"))
##   setGeneric("SimPatterns<-", function(object, value) standardGeneric("SimPatterns<-"))

## setReplaceMethod("SimPatterns", "RtreemixSim",
##                  function(object, value) {object@SimPatterns <- value; object})

## if(!isGeneric("SamplingMode<-"))
##   setGeneric("SamplingMode<-", function(object, value) standardGeneric("SamplingMode<-"))

## setReplaceMethod("SamplingMode", "RtreemixSim",
##                  function(object, value) {object@SamplingMode <- value; object})

## if(!isGeneric("SamplingParam<-"))
##   setGeneric("SamplingParam<-", function(object, value) standardGeneric("SamplingParam<-"))

## setReplaceMethod("SamplingParam", "RtreemixSim",
##                  function(object, value) {object@SamplingParam <- value; object})

## if(!isGeneric("WaitingTimes<-"))
##   setGeneric("WaitingTimes<-", function(object, value) standardGeneric("WaitingTimes<-"))

## setReplaceMethod("WaitingTimes", "RtreemixSim",
##                  function(object, value) {object@WaitingTimes <- value; object})

## if(!isGeneric("SamplingTimes<-"))
##   setGeneric("SamplingTimes<-", function(object, value) standardGeneric("SamplingTimes<-"))

## setReplaceMethod("SamplingTimes", "RtreemixSim",
##                  function(object, value) {object@SamplingTimes <- value; object})

## A method for getting the number of draws in the simulation (The size of the drawn (simulated) sample).
if(!isGeneric("noDraws"))
  setGeneric("noDraws", function(object) standardGeneric("noDraws"))

setMethod("noDraws", signature = signature(object = "RtreemixSim"),
          definition = function(object) 
          return(sampleSize(object@SimPatterns))      
          )

## A method for getting the mixture model. 
if(!isGeneric("getModel"))
  setGeneric("getModel", function(object) standardGeneric("getModel"))
  
setMethod("getModel", signature = signature(object = "RtreemixSim"),
          definition = function(object) {
            if(length(Resp(object)) != 0) {
              return(new("RtreemixModel",
                         ParentData = getData(object),  
                         Weights = Weights(object),
                         WeightsCI = WeightsCI(object),
                         Resp = Resp(object),
                         CompleteMat = CompleteMat(object),
                         Star = Star(object),
                         Trees = Trees(object)))
            } else {
               return(new("RtreemixModel",
                          Weights = Weights(object),
                          WeightsCI = WeightsCI(object),
                          Star = Star(object),
                          Trees = Trees(object)))
             }
          })

## Printing function for the class RtreemixSim.
setMethod("print", "RtreemixSim", function(x, ...) .printRtreemixSim(x))
setMethod("show", "RtreemixSim", function(object) .printRtreemixSim(x = object))

.printRtreemixSim <- function(x) {
  
  if(sampleSize(SimPatterns(x)) == 0)               
    return(cat("Empty object!\n"))
  
  cat("RtreemixSim object with", noDraws(x), "patterns.\n\n")
  
  cat("The sample data is given by: \n\n")
  print(SimPatterns(x))
  
  if(!identical(SamplingMode(x), character(0))) {
    if(SamplingMode(x) == "exponential")
      cat("\nExponential sampling times are used in the simulation of the waiting process.\n")
    else
      cat("\nConstant sampling times (in days) are used in the simulation of the waiting process.\n")
  }  
  
  if(!identical(SamplingParam(x), numeric(0))) {
    cat("The sampling parameter is:\n")
    print(SamplingParam(x))
  }
  
  if(length(WaitingTimes(x)) !=  0) {
    cat("\nThe waiting times of the first 10 patterns are:\n\n")    
    if(length(WaitingTimes(x)) > 10) { 
      a <- WaitingTimes(x)[1:10]
      names(a) <- Patients(SimPatterns(x))[1:10]
    } else {
      a <- WaitingTimes(x)
      names(a) <- Patients(SimPatterns(x))      
    }
    print(a)
  }
  
  if(length(SamplingTimes(x)) != 0) {
    cat("\nThe sampling times of the first 10 patterns are:\n\n")
    if(length(SamplingTimes(x)) > 10) {
      a <- SamplingTimes(x)[1:10]
      names(a) <- Patients(SimPatterns(x))[1:10]
    } else {
      a <- SamplingTimes(x)
      names(a) <- Patients(SimPatterns(x))
    }
    print(a)    
  }            
}

################################################################################

## Constructor for creating a new object of the class RtreemixStats.
setMethod("initialize", "RtreemixStats",
          definition = function(.Object, Data, Model, LogLikelihoods, 
            WLikelihoods) {
            if(missing(Data) || (class(Data) != "RtreemixData"))
              stop("Please specify the data for the likelihood calculation!")
            else {
              ## Initializing the slots inherited from the parent class (RtreemixData)
              .Object@Sample <- Sample(Data)
              .Object@Patients <- Patients(Data)
              .Object@Events <- Events(Data)
              .Object@Description <- Description(Data)
            } 
            ## Initializing the slots from RtreemixStats
            if(missing(Model))
              stop("Please provide the underlying model for calculating the likelihoods.")

            if(class(Model) !=  "RtreemixModel")
              stop("The specified model should be an object of type RtreemixModel.")
            
            .Object@Model <- Model
            
            if(missing(WLikelihoods) || !is.matrix(WLikelihoods))
              stop("Please specify the weighted likelihoods of the provided data as a matrix!")
            else {
              if ((nrow(WLikelihoods) != sampleSize(Data)) || (ncol(WLikelihoods) != numTrees(Model)))
                stop("The dimensions of the weighted likelihoods matrix are not correct!")
              
              .Object@WLikelihoods <- WLikelihoods
              if(ncol(.Object@WLikelihoods) != 0)
                colnames(.Object@WLikelihoods) <-  paste("T", 1:ncol(.Object@WLikelihoods), sep = "")
            }

            if(missing(LogLikelihoods) || !is.numeric(LogLikelihoods))
              stop("Please specify the log-likelihoods of the provided data as a numeric vector!")
            else {
              if (length(LogLikelihoods) != sampleSize(Data))
                stop("The log-likelihoods are missing for some patterns!")
              
              .Object@LogLikelihoods <- LogLikelihoods
            }           
            
            .Object
          })
          
## Accessor functions for the class RtreemixStats.
if(!isGeneric("Model"))
  setGeneric("Model", function(object) standardGeneric("Model"))

setMethod("Model", "RtreemixStats", function(object) object@Model)

if(!isGeneric("LogLikelihoods"))
  setGeneric("LogLikelihoods", function(object) standardGeneric("LogLikelihoods"))

setMethod("LogLikelihoods", "RtreemixStats", function(object) object@LogLikelihoods)

if(!isGeneric("WLikelihoods"))
  setGeneric("WLikelihoods", function(object) standardGeneric("WLikelihoods"))

setMethod("WLikelihoods", "RtreemixStats", function(object) object@WLikelihoods)

## Replacement of the slots of this class is not allowed.
## Replacement functions for the class RtreemixStats.          
## if(!isGeneric("Model<-"))
##   setGeneric("Model<-", function(object, value) standardGeneric("Model<-"))

## setReplaceMethod("Model", "RtreemixStats",
##                  function(object, value) {object@Model <- value; object})

## if(!isGeneric("LogLikelihoods<-"))
##   setGeneric("LogLikelihoods<-", function(object, value) standardGeneric("LogLikelihoods<-"))
  
## setReplaceMethod("LogLikelihoods", "RtreemixStats",
##                  function(object, value) {object@LogLikelihoods <- value; object})

## if(!isGeneric("WLikelihoods<-"))
##   setGeneric("WLikelihoods<-", function(object, value) standardGeneric("WLikelihoods<-"))

## setReplaceMethod("WLikelihoods", "RtreemixStats",
##                  function(object, value) {object@WLikelihoods <- value; object})

## A method for getting the data.
if(!isGeneric("getData"))
  setGeneric("getData", function(object) standardGeneric("getData"))
  
setMethod("getData", signature = signature(object = "RtreemixStats"),
          definition = function(object) 
          return(new("RtreemixData",
                     Sample = Sample(object),
                     Patients = Patients(object),
                     Events = Events(object)))      
          )

## A method for getting the responsibilities as KxN matrix.
## (resp[k][i] == responsibility of tree k for pattern i)
if(!isGeneric("getResp"))
  setGeneric("getResp", function(object) standardGeneric("getResp"))
  
setMethod("getResp", signature = signature(object = "RtreemixStats"),
          definition = function(object) 
          return(as.matrix(t(object@WLikelihoods/apply(object@WLikelihoods, 1, sum))))      
          )

## Printing function for the class RtreemixStats.
setMethod("print", "RtreemixStats", function(x, ...) .printRtreemixStats(x))
setMethod("show", "RtreemixStats", function(object) .printRtreemixStats(x = object))

.printRtreemixStats <- function(x) {
  cat("RtreemixStats object with", sampleSize(x), "samples.\n\n")
  if(sampleSize(x) != 0) {
    cat("The first 10 samples from the dataset: \n\n")
    a <- cbind(rep(1, sampleSize(x)), Sample(x))
    rownames(a) <- Patients(x)
    colnames(a) <- Events(x)
    if(nrow(a) > 10) {
      a <- a[1:10, ]
    }
    print(a)
  } 

  if(length(LogLikelihoods(x)) != 0) {
    cat("\nThe log-likelihoods of the first 10 samples:\n\n")
    a <- LogLikelihoods(x)
    names(a) <- Patients(x)
    if(length(LogLikelihoods(x)) > 10) { 
      a <- LogLikelihoods(x)[1:10]
      names(a) <- Patients(x)[1:10]
    }
    print(a)
  }
  if(length(WLikelihoods(x)) != 0) {   
    cat("\nThe weighted likelihoods of the first 10 samples:\n\n")
    a <- WLikelihoods(x)
    rownames(a) <- Patients(x)  
    if(nrow(WLikelihoods(x)) > 10) { 
      a <- WLikelihoods(x)[1:10, ]
      rownames(a) <- Patients(x)[1:10]
       }
    print(a)
  }
}

################################################################################
## Constructor for creating a new object of the class RtreemixGPS.
setMethod("initialize", "RtreemixGPS",
          definition = function(.Object, Data, Model,
            SamplingMode, SamplingParam, GPS, gpsCI) {
            if(missing(Data) || (class(Data) != "RtreemixData"))
              stop("Please specify the data for the GPS calculation!")
            else {
              ## Initializing the slots inherited from the parent classes (RtreemixData and RtreemixModel) 
              
              ## callNextMethod(.......)
              .Object@Sample <- Sample(Data)
              .Object@Patients <- Patients(Data)
              .Object@Events <- Events(Data)
              .Object@Description <- Description(Data)
            } 
            ## Initializing the slots from RtreemixGPS
            if(class(Model) != "RtreemixModel")
              stop("The model has to be given as an RtreemixModel object!")
            
            if(missing(Model))
              stop("Please specify the underlying model!")
            else
              .Object@Model <- Model
            
            if(missing(SamplingMode))
              .Object@SamplingMode <- character(0)
            else {
              if ((SamplingMode != "exponential") && (SamplingMode != "constant"))
                stop("The sampling mode has to be exponential or constant!")
              
              .Object@SamplingMode <- SamplingMode
            }
            
            if(missing(SamplingParam))
              .Object@SamplingParam <- numeric(0)
            else {        
              .Object@SamplingParam <- SamplingParam
            }            
            
            if(missing(GPS))
              .Object@GPS <- numeric(0)
            else {
              if (length(GPS) != sampleSize(Data))
                stop("The GPS values are missing for some patterns!")
              
              .Object@GPS <- GPS
              names(.Object@GPS) <- Patients(Data)
            }
            
            if(missing(gpsCI))
              .Object@gpsCI <- matrix(numeric(0), 0, 0)
            else {
              if(missing(GPS))
                stop("Confidence intervals for the GPS values cannot be specified when the GPS values are missing!")
              
              if (length(GPS) != nrow(gpsCI))
                stop("The GPS confidence intervals are missing for some GPS values!")
              
              .Object@gpsCI <- gpsCI
              if (nrow(.Object@gpsCI) != 0)
                rownames(.Object@gpsCI) <- paste("GPS", 1:nrow(.Object@gpsCI), sep = "")

              if(ncol(.Object@gpsCI) != 2)
                stop("The confidence intervals for the GPS values are not correctly specified!")
              else 
                colnames(.Object@gpsCI) <- c("lower", "upper")
            }              
            
            .Object
          })  

## Accessor functions for the slots of the class RtreemixGPS.
if(!isGeneric("Model"))
  setGeneric("Model", function(object) standardGeneric("Model"))

setMethod("Model", "RtreemixGPS", function(object) object@Model)

if(!isGeneric("SamplingMode"))
  setGeneric("SamplingMode", function(object) standardGeneric("SamplingMode"))

setMethod("SamplingMode", "RtreemixGPS", function(object) object@SamplingMode)

if(!isGeneric("SamplingParam"))
  setGeneric("SamplingParam", function(object) standardGeneric("SamplingParam"))

setMethod("SamplingParam", "RtreemixGPS", function(object) object@SamplingParam)

if(!isGeneric("GPS"))
  setGeneric("GPS", function(object) standardGeneric("GPS"))

setMethod("GPS", "RtreemixGPS", function(object) object@GPS)

if(!isGeneric("gpsCI"))
  setGeneric("gpsCI", function(object) standardGeneric("gpsCI"))

setMethod("gpsCI", "RtreemixGPS", function(object) object@gpsCI)

## Replacement of the slots of this class is not allowed.
## Replacement functions for the slots of the class RtreemixGPS.          
## if(!isGeneric("Model<-"))
##   setGeneric("Model<-", function(object, value) standardGeneric("Model<-"))

## setReplaceMethod("Model", "RtreemixGPS",
##                  function(object, value) {object@Model <- value; object})

## if(!isGeneric("SamplingMode<-"))
##   setGeneric("SamplingMode<-", function(object, value) standardGeneric("SamplingMode<-"))

## setReplaceMethod("SamplingMode", "RtreemixGPS",
##                  function(object, value) {object@SamplingMode <- value; object})

## if(!isGeneric("SamplingParam<-"))
##   setGeneric("SamplingParam<-", function(object, value) standardGeneric("SamplingParam<-"))

## setReplaceMethod("SamplingParam", "RtreemixGPS",
##                  function(object, value) {object@SamplingParam <- value; object})

## if(!isGeneric("GPS<-"))
##   setGeneric("GPS<-", function(object, value) standardGeneric("GPS<-"))

## setReplaceMethod("GPS", "RtreemixGPS",
##                  function(object, value) {object@GPS <- value; object})

## if(!isGeneric("gpsCI<-"))
##   setGeneric("gpsCI<-", function(object, value) standardGeneric("gpsCI<-"))

## setReplaceMethod("gpsCI", "RtreemixGPS",
##                  function(object, value) {object@gpsCI <- value; object})

## A method for getting the data.
if(!isGeneric("getData"))
  setGeneric("getData", function(object) standardGeneric("getData"))
  
setMethod("getData", signature = signature(object = "RtreemixGPS"),
          definition = function(object) 
          return(new("RtreemixData",
                     Sample = Sample(object),
                     Patients = Patients(object),
                     Events = Events(object)))      
          )

## Printing function for the class RtreemixGPS.
setMethod("print", "RtreemixGPS", function(x, ...) .printRtreemixGPS(x))

setMethod("show", "RtreemixGPS", function(object) .printRtreemixGPS(x = object))

.printRtreemixGPS <- function(x) {

  if((sampleSize(x) == 0) || (length(GPS(x)) == 0))               
    return(cat("Empty object!\n"))
  
  cat("RtreemixGPS object with", sampleSize(x), "patterns.\n")
  
  if(!identical(SamplingMode(x), character(0))) {
    if(SamplingMode(x) == "exponential")
      cat("\nExponential sampling times are used in the simulation of the waiting process.\n")
    else
      cat("\nConstant sampling times (in days) are used in the simulation of the waiting process.\n")
  }  

  if(!identical(SamplingParam(x), numeric(0))) {
    cat("The sampling parameter is:\n")
    print(SamplingParam(x))
  }
  
  if(length(GPS(x)) !=  0) {
    if(length(gpsCI(x)) != 0) {
      cat("\n\nThe GPS values and their confidence intervals of the first 10 patterns are:\n\n")
      
      if(length(GPS(x)) > 10) { 
        g <- GPS(x)[1:10]
        ci <- gpsCI(x)[1:10, ]
      } else {
        g <- GPS(x)
        ci <- gpsCI(x)
      }
      
      a <- cbind(rep(1, sampleSize(x)), Sample(x))
      rownames(a) <- Patients(x)    
      if(nrow(a) > 10) {
        a <- a[1:10, ]
        rownames(a) <- Patients(x)[1:10]
      }
      df <- data.frame(a, g, ci)
      colnames(df) <- c(Events(x), "GPS", "lower", "upper")
      rownames(df) <- rownames(a) 
      print(df)
      
    } else {
      cat("\n\nThe GPS values of the first 10 patterns are:\n\n")
      
      if(length(GPS(x)) > 10) 
        g <- GPS(x)[1:10]
      else
        g <- GPS(x)
      
      a <- cbind(rep(1, sampleSize(x)), Sample(x))
      rownames(a) <- Patients(x)    
      if(nrow(a) > 10) {
        a <- a[1:10, ]
        rownames(a) <- Patients(x)[1:10]
      }
      a <- data.frame(a, g)
      colnames(a) <- c(Events(x), "GPS")
      print(a)
      
    }    
  }
}

################################################################################  
