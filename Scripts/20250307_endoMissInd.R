# packages
library(here)
library(sna)
library(ergm)

# loadinergm
load(here("Data", "20231006_missNetsEnMasse.RData"))

# good ol london gangs
trueAdj = adjMatList[[6]]
trueNet = as.network(trueAdj, directed= FALSE)

# loading some functions
source(here("Scripts", "20230811_missNetFunctions.R"))

# initialising the missingness indicator
initMiss = initMissAdj(n = nrow(trueAdj), propMiss = 0.3)
initMissNet = as.network(initMiss, directed = FALSE)

## setting up a specific missingness indicator

# specifying the endogenous missingness model
missModel = initMissNet ~ edges + dyadcov(trueAdj) + gwdegree(decay = 0.69, fixed = TRUE) + gwesp(decay = 0.69, fixed = TRUE)

# constrained networks
missConNets = simulate(object =  missModel,
         coef = c(0, 0.2, 0.4, 0.5),
         nsim = 100,
          control = control.simulate.formula(MCMC.burnin=20000,MCMC.interval=100),
         constraints =~ edges)

# randomly choose a draw
chosenMissInd = missConNets[[14]]

save(chosenMissInd, file = "20250307_chosenEndoMissInd.RData")