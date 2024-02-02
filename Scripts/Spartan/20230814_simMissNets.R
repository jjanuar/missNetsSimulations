## Script to be used for Spartan to generate a lot of simulated missingness matrices
## Grabbing some variables from the SLURM to use for the script
# what do I need?
propMissInd = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ergmCoefSet = as.numeric(Sys.getenv("ergmCoef"))
latentBeta = as.numeric(Sys.getenv("latentBeta"))

# as I'm using the array index for the proportion of missingness, these can't be floats
# so I need to define my proportions of missingness within the script
propMisses = c(0.1, 0.35, 0.6)
chosenPropMiss = propMisses[propMissInd]

# for the ergm coefficient set, this is avector so I'm going to define specfic sets in lists within R
# the same index will be used for both the model specification and the coefficients used.
# this can be revised into two separate indices (environment variables) if necessary
ergmCoefList = list(c(0, 0.4, 0.5, 0.8, 0.2),
                    c(0, 2, 2, 0, 0),
                    c(0, -2, 1, 0, 0),
                    c(0, 2, 2, 1, 1),
                    c(0, -2, 1, -1, -1))
ergmSpecList = list(initMissNet ~ edges + 
                    gwdegree(0.69, fixed = TRUE) +
                    gwesp(0.69, fixed = TRUE) +
                    dyadcov(trueNet) +
                    nodecov('degree'))

chosenErgmCoef = ergmCoefList[[ergmCoefSet]]
chosenErgmSpec = ergmSpecList[[1]]

# the last bit is the beta coefficient value for the latent variables
# an additional variable could be used if I wanted to change the model specification here (e.g., use a different tie function or add attributes)
# but that can be revised if necessary as well.
# doesn't need additional code

## Checks
print(paste("Your chosen proportion of missingness is", chosenPropMiss))
cat("Your chosen ERGM coefficient values are", chosenErgmCoef, "\n")
print(chosenErgmSpec)
print(paste("Your chosen latent model beta parameter value is", latentBeta))

# packages
library(here)
library(stats)
library(testthat)
library(sna)
library(network)
library(ergm)

## Loading functions
source(here("Scripts",'20230811_missNetFunctions.R'))

## Load networks
# this will be a k-long list of matrix objects
# for the 3 methods in this script for generating missingness,
# each will have k outputs (e.g., one per network)
# 6th October 2023 update: made it so both model specifications use altstar instead of gwdeg
load(here("Data", "20231006_missNetsEnMasse.RData"))

# get the vertex sizes
vertexSizes = unlist(lapply(adjMatList, FUN = nrow))

### Independent missingness
indepMissMats = lapply(X = vertexSizes, 
                       FUN = initMissAdj, 
                       propMiss = chosenPropMiss,
                       numSims = 1000)

### ERGM specified missingness
## generate the initial missingness matrix
initErgmAdjs = lapply(X = vertexSizes, 
                      FUN = initMissAdj, 
                      propMiss = chosenPropMiss)

# turn into a network object
initErgmNets = lapply(X = initErgmAdjs,
                      FUN = network::as.network,
                      directed = FALSE)


# put the degrees into each network object
# no idea how to vectorise this since I need an index for the degree
for( i in 1:length(adjMatList) ){
  
  # put in the degree
  initErgmNets[[i]] %v% 'degree' = rowSums(adjMatList[[i]])
  
}

# another loop for the ERGM simulations
# not sure how to vectorise this either
missErgmSimNets = list()

for( nets in 1:length(adjMatList) ){
  
  # choose a specific initial network
  initMissNet = initErgmNets[[nets]]

  # specify which network is the 'true' network for the dyadic covariate
  trueNet = adjMatList[[nets]]
  
  # run the simulation
  missErgmSimNets[[nets]] = simMissNet(model = chosenErgmSpec, coef = chosenErgmCoef)
  
}

### Latent model missingness
# I would need k latent spaces
# then simulate from the distance matrices of these latent spaces

# latent model specification (can change when necessary)
latentMissModel = initMissNet ~ edges + dyadcov(distanceMat)

# generate latent variable coordinates
latentCoords = lapply(X = vertexSizes, FUN = simLatentCoords, noDim = 2)

# calculate distance matrices
distanceMats = lapply(X = latentCoords, FUN = function(x){as.matrix(stats::dist(x))})

## generate the initial missingness matrix
initLatentAdjs = lapply(X = vertexSizes, 
                        FUN = initMissAdj, 
                        propMiss = chosenPropMiss)

# turn into a network object
initLatentNets = lapply(X = initLatentAdjs,
                        FUN = network::as.network,
                        directed = FALSE)

# and then simulate
# not sure how to vectorise this either
missLatentSimNets = list()

for( nets in 1:length(adjMatList) ){
  
  # choose a specific initial network
  initMissNet = initLatentNets[[nets]]

  # specify which distance matrix to use
  distanceMat = distanceMats[[nets]]
  
  # run the simulation
  missLatentSimNets[[nets]] = simMissNet(model = latentMissModel, coef = c(0, latentBeta))
  
}

### Save objects
save(adjMatList, indepMissMats, missErgmSimNets, missLatentSimNets, file = here("Output", "20230814_simMissNets", paste(format(Sys.Date(), "%Y%m%d"), "_simMissNets_","propMiss", chosenPropMiss*10, ".RData", sep = "")))


