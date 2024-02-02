## Spartan script to deplete and reestimate the ERGMs
## To be used after generating missing indicator matrices from 20230814_simMissNets.R


# load packages
library(ergm)
library(sna)
library(network)
library(here)
library(testthat)

# loading in some other attributes
# 6th October 2023 update: made it so both model specifications use altstar instead of gwdeg
load(here("Data", "20231006_missNetsEnMasse.RData"))

# loading in the script of miss net functions to degrade the network
source(here("Scripts", "20230811_missNetFunctions.R"))


### Environment variables and command arguments
## Environment variables
# best for setting them to specific values for various runs
missSaveInd = as.numeric(Sys.getenv("missSaveInd"))

# Array index for the trials per chosen model
trialIndex = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# proporton of missingness (0.1/0.35/0.6)
propMissValue = as.numeric(Sys.getenv("propMissVal"))

print(propMissValue)

## Command arguments
# best if I want to use multiple values in multiple runs (e.g., to loop across different networks)
# grabbing command arguments
command_args <- commandArgs(trailingOnly = TRUE)

# print command args
print(command_args)


# specify the necessary arguments
# the network that is used (from 1 to 6)
networkIndex = as.numeric(command_args[1])

# the missingness model used (Independent/ERGM/Latent = 1/2/3)
missModelInd = as.numeric(command_args[2])


# vector of missingness types
missModels = c("Independent", "ERGM", "Latent")
missSave = c(NA, 0)
chosenMissSave = missSave[missSaveInd]

# an additional label for what the missing values are saved as for file management
missSaveLabel = c("Peter", "Todd")
chosenMissLabel = missSaveLabel[missSaveInd]


# print it
print(paste("The trial index is:", trialIndex))
print(paste("The network chosen is:", networkIndex))
print(paste("The missing model chosen is:", missModels[missModelInd]))
print(paste("Missing values are saved as:", chosenMissSave))
print(paste("The proportion of missingness is:", propMissValue))


## Load a specific set of missingness simlulations
load(here("Output", "20230814_simMissNets", paste("20230825_simMissNets_propMiss", propMissValue*10,".RData", sep = "")))

## choose a missingness matrix
# index the network based on the trial value, floored to avoid decimals
trialIndices = floor(seq(from = 501, to = 1000, length.out = 50))
chosenIndex = trialIndices[trialIndex]

# make sure everything is a matrix and put in a list so it can be indexed
chosenIndepMiss = indepMissMats[[networkIndex]][chosenIndex,,]
chosenErgmMiss = as.matrix(missErgmSimNets[[networkIndex]][[chosenIndex]])
chosenLatentMiss = as.matrix(missLatentSimNets[[networkIndex]][[chosenIndex]])

# put all these matrices in a list to be indexed
missAdjList = list(chosenIndepMiss, chosenErgmMiss, chosenLatentMiss)


# choose one, vary this with command arguments/environment variables
chosenMissAdj = missAdjList[[missModelInd]]

## Deplete the network here
chosenNetwork = degradeAdj(trueNet = adjMatList[[networkIndex]],
                           missAdj = chosenMissAdj,
                           missSave = chosenMissSave)


# turn into a network and estimate
chosenNetObj = as.network(chosenNetwork, directed = FALSE)

# a specific branch for the london gangs to add attributes
if(networkIndex == 6){
  
  # add in attributes
  chosenNetObj %v% 'Age' = londonGangsAtt$Age
  chosenNetObj %v% 'Arrests' = londonGangsAtt$Arrests
  chosenNetObj %v% 'Birthplace' = londonGangsAtt$Birthplace
  chosenNetObj %v% 'Convictions' = londonGangsAtt$Convictions
  chosenNetObj %v% 'Prison' = londonGangsAtt$Prison
}


## Choose model
# The thing is, I only use the other model for the london gangs, so...
if(networkIndex == 6){
  modelIndex = 2
} else {
  modelIndex = 1
}

# and choose the model
chosenModel = modelList[[modelIndex]]

# print the chosen network and model
print(paste("The chosen network is", networkLabels[[networkIndex]]))
print(paste("The chosen model is", chosenModel))

# check if the model's got any degenerate parameters
summary(chosenModel)

# initial parameters, obtained from the 'true' model estimations
initList = list(c(-2.9346, -2.1647, 4.696),
                c(-1.7506, -2.0996, 4.0521),
                c(-14.6097, -1.7190, 9.7026),
                c(-0.7729, -1.8808, 3.1026),
                c(-4.0564, -0.1408, 1.4859),
                c(-6.311, -0.509, 1.291, 0.084, -0.177, 0.066, -0.01, -0.08, -0.073, 0.105, 0.835))

chosenInit = initList[[networkIndex]]

# run ergm
modelres = ergm(chosenModel, 
                control = control.ergm(init = chosenInit, MCMC.burnin = 20000, MCMC.interval = 2000))

# save results
save(modelres, file = here("Output", "20230825_simMissReest", chosenMissLabel,
                           paste(format(Sys.Date(), "%Y%m%d"),"_missNetReest_net", networkIndex, "_missModel", missModelInd, "_prop", propMissValue*10,"_trial", trialIndex, ".RData", sep = "")))
