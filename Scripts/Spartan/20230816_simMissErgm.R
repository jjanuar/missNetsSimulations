### A version of 20230814_simMissNets that specifically does various ERGM specifications

## Script to be used for Spartan to generate a lot of simulated missingness matrices
## Grabbing some variables from the SLURM to use for the script
# what do I need?
propMissInd = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ergmCoefSet = as.numeric(Sys.getenv("ergmCoef"))
ergmModelSpec = as.numeric(Sys.getenv("ergmSpec"))

# as I'm using the array index for the proportion of missingness, these can't be floats
# so I need to define my proportions of missingness within the script
propMisses = c(0.1, 0.35, 0.6)

# instead of increments of 10%, can divide it up to a high, medium, and low (e.g., 0.6, 0.35, 0.1)

chosenPropMiss = propMisses[propMissInd]

# for the ergm coefficient set, this is avector so I'm going to define specfic sets in lists within R
# the same index will be used for both the model specification and the coefficients used.
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
chosenErgmSpec = ergmSpecList[[ergmModelSpec]]

## Checks
print(paste("Your chosen proportion of missingness is", chosenPropMiss, "\n"))
cat("Your chosen ERGM coefficient values are", chosenErgmCoef, "\n")
print(paste("Your chosen ERGM specification is", chosenErgmSpec, "\n"))

# packages, the necessary ones will be loaded when the functions are ran
library(here)
library(stats)

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

### Save objects
save(adjMatList, missErgmSimNets, chosenErgmCoef, 
     file = here("Output", "20230814_simMissNets", 
       paste(format(Sys.Date(), "%Y%m%d"), "_simMissErgmNets_ergmCoef_",ergmCoefSet,"_propMiss", chosenPropMiss*10, ".RData", sep = "")))


