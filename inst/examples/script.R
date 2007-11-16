library(Rtreemix)


######################################################################
## Create an RtreemixData object from a given file. 
data <- new("RtreemixData", File = "D:/package/Rtreemix/inst/examples/treemix.pat")
show(data)
## Create an RtreemixData object from a randomly generated RtreemixModel object.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
data <- sim(model = rand.mod, no.draws = 300)
show(data)
## See your matrix, the event names, the patient IDs and the description of the object data.
Sample(data)
Events(data)
Patients(data)
Description(data)

## Get the sample size and the number of genetic events.
sampleSize(data)
eventsNum(data)
######################################################################
## Generate a random RtreemixModel object with 3 components.
rand.mod <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
show(rand.mod)

## Create an RtreemixModel object by fitting model to the given data.
mod <- fit(data = data, K = 3, equal.edgeweights = TRUE, noise = TRUE)
show(mod)

## Create an RtreemixData by drawing data from a specified mixture model.
draws <- sim(model = mod, no.draws = 100)
show(draws)

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
mod.boot <- bootstrap(data = data, K = 2, equal.edgeweights = TRUE, B = 100) 

## See the confidence intervals for the mixture parameters (the weights).
WeightsCI(mod.boot)
## See the confidence intervals of the conditional probabilities assigned to the edges.
edgeData(getTree(mod.boot, 2), attr = "ci")

######################################################################
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
######################################################################
## Create an RtreemixGPS object by calculating the GPS for all possible patterns.
modGPS.all <- gps(model = mod, no.sim = 1000)
show(modGPS.all)
## Create an RtreemixGPS object by calculating the GPS for the data based on the model mod.
modGPS <- gps(model = mod, data = data, no.sim = 1000)
show(modGPS)

## See the slots from the RtreemixGPS object. 
#Model(modGPS)
#SamplingMode(modGPS)
#SamplingParam(modGPS)
GPS(modGPS)
## See data.
getData(modGPS)

## Create an RtreemixGPS object by calculating GPS values for a given dataset
## and their 95% confidence intervals using the bootstrap method.
modGPS2 <- confIntGPS(data = data, K = 2, B = 100)
show(modGPS2)

## See the GPS values for the object modGPS2 and their confidence intervals.
GPS(modGPS2)
gpsCI(modGPS2)
######################################################################
## Create an RtreemixSim object by simulating patterns with their sampling and waiting times from a given mixture model.
sim.data <- sim(model = mod, sampling.mode = "exponential", sampling.param = 1, no.sim = 100)
show(sim.data)

## See the slots from the RtreemixSim object.
SimPatterns(sim.data)
SamplingMode(sim.data)
SamplingParam(sim.data)
WaitingTimes(sim.data)
SamplingTimes(sim.data)
## See model.
getModel(sim.data)
######################################################################
## Generate two random mixture models
mix1 <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
mix2 <- generate(K = 3, no.events = 9, noise.tree = TRUE, prob = c(0.2, 0.8))
## Compare the tree topologies of two mixture models using the number
## of distinct edges.
comp.models(mixture1 = mix1, mixture2 = mix2)
## Compare the tree topologies of two mixture models using the number
## of distinct edges and the levels of the events in the trees.
comp.models.levels(mixture1 = mix1, mixture2 = mix2)
## Inspect the diversity of the nontrivial tree components in a given model
## using the number of distinct edges as dissimilarity measure
comp.trees(model = mix1)
comp.trees(model = mix2)
## Inspect the diversity of the nontrivial tree components in a given model
## using the number of distinct edges and the levels of the events in
## the treesas dissimilarity measure
comp.trees.levels(model = mix1)
comp.trees.levels(model = mix2)
## Stability analysis
stability.sim(no.trees = 3, no.rands = 5, no.sim = 4, no.draws = 100)
######################################################################

