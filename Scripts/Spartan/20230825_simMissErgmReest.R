## Spartan script to deplete and reestimate the ERGMs
## To be used after generating missing indicator matrices from 20230814_simMissNets.R
## This version is adapted to only estimate the ergm-depleted missingness
# because only the ERGM has many varying coefficient sets


### Updated 20 Sep 2023 to include Peter or Todd because at this point I had only
### done peter'd re-estimations (missing as NA, not missing as 0).

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
# because this is the ERGM-only missingness script, this index will always be 2
missModelInd = 2


# an additional index specifically for this script for which ergm coefficient set to use
ergmCoefInd = as.numeric(Sys.getenv("ergmCoef"))
ergmCoefList = list(c(0, 0.4, 0.5, 0.8, 0.2),
                    c(0, 2, 2, 0, 0),
                    c(0, -2, 1, 0, 0),
                    c(0, 2, 2, 1, 1),
                    c(0, -2, 1, -1, -1))

## Reminder: the five parameter here are as follows:
# Parameter 1 = edges, kept to 0 because simulation of the missingness was constrained on its 'density'
# Parameter 2 = GWDegree (not alt star!), high +ve is high centralisation (degree variance), high -ve is high clustering/low centralisation (low degree variance)
# Parameter 3 = GWESP (same as alt triangles), high +ve is high 'clustering' (an edge closes multiple triangles), high -ve is the 'opposite'- long paths, less triangles
# Parameter 4 = Entrainment (true network as dyadic covariate). High +ve means that if an edge exists, it's more likely to be missing, high -ve means null ties more likely to be missing
# Parameter 5 = Degree node covariate. +ve means that high degree nodes (from X) are more likely to have missing tie variables.

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
print(paste("The ergm coefficient set is:", ergmCoefList[[ergmCoefInd]]))

## Load a specific set of missingness simlulations
load(here("Output", "20230814_simMissNets", paste("20231011_simMissErgmNets_ergmCoef_", ergmCoefInd,"_propMiss", propMissValue*10,".RData", sep = "")))

## choose a missingness matrix
# index the network based on the trial value, floored to avoid decimals
trialIndices = floor(seq(from = 501, to = 1000, length.out = 50))
chosenIndex = trialIndices[trialIndex]

# make sure everything is a matrix and put in a list so it can be indexed
chosenErgmMiss = as.matrix(missErgmSimNets[[networkIndex]][[chosenIndex]])


# choose one, vary this with command arguments/environment variables
chosenMissAdj = chosenErgmMiss

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
                           paste(format(Sys.Date(), "%Y%m%d"),"_missErgmNetReest_net", networkIndex, "_coefSet", ergmCoefInd, "_prop", propMissValue*10,"_trial", trialIndex, ".RData", sep = "")))
