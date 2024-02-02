## Spartan script to estimate 'true' models
## To be used in conjunction with the specifications in 20230726_missNetFunctions.Rmd

# load packages
library(ergm)
library(sna)
library(network)
library(here)

## loading in data (networks and model specifications)
# 6th October 2023 update: made it so both model specifications use altstar instead of gwdeg
load(here("Data", "20231006_missNetsEnMasse.RData"))

# specify the necessary arguments
#modelIndex = as.numeric(Sys.getenv("ModelInd"))
trialIndex = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# grabbing command arguments
command_args <- commandArgs(trailingOnly = TRUE)

# print command args
print(command_args)

# trial index because manually running this one at a time sucks.
networkIndex = as.numeric(command_args[1])

# print it
print(paste("The trial index is:", trialIndex))
print(paste("The network index is:", networkIndex))


# choose a network
chosenNetwork = adjMatList[[networkIndex]]

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

# trying out without initial parameters
# network 1, c(-11.6859, 4.0456, 4.7610)
# network 4, c(-6.7, 4.3, 1.5)
# initial parameters
initList = list(0,
                c(-1.3925, -2.0863, 3.8478),
                c(-16.1072, -1.2172, 9.4760),
                c(-0.7524, -1.8705, 3.0747),
                c(-4.0956, -0.1222, 1.4765),
                c(-6.311, -0.509, 1.291, 0.084, -0.177, 0.066, -0.01, -0.08, -0.073, 0.105, 0.835))

chosenInit = initList[[networkIndex]]

# run ergm
modelres = ergm(chosenModel, 
                control = control.ergm(init = chosenInit, MCMC.burnin = 20000, MCMC.interval = 20000)) 

# save results
save(modelres, file = here("Output", "20230726_missNetsTrueModels", 
                           paste("20230726_missNetTrueModel_net", networkIndex, "_model", modelIndex, "_trial", trialIndex, ".RData", sep = "")))
