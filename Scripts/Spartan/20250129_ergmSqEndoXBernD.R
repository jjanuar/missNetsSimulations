### A Spartan-friendly version of a script to run the ERGM squared sampler on a single depleted network
## ideally parallelised so

## Grabbing some variables from the SLURM to use for the script
# what do I need?
entrValueInd = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# ergmModelSpec = as.numeric(Sys.getenv("ergmSpec"))


# specify the entrainment values
entrValues = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)

# and specify the chosen entrainment value
chosenEntrValue = entrValues[entrValueInd]

### packages
library(here)
library(ergm)
library(mvtnorm)
library(sna)

# sourcing the ergm squared function
source(here("Scripts", "20241021_ergmSqFunc.R"))

# some functions
source(here("Scripts", "20230811_missNetFunctions.R"))

# Load in some dataset
load(here("Data", "20231006_missNetsEnMasse.RData"))

# London gangs since it's an easy dataset with some attributes
trueAdj = adjMatList[[6]]

# get the true statistics for the given estimation model for reference
trueNet = as.network(trueAdj, directed = FALSE)
trueNet %v% 'Age' = scale(londonGangsAtt$Age)[,1]

# this time, the true model needs to have some markov statistics, so we can use the specification that was previously ran
trueModel = trueNet ~ edges +
  gwdegree(0.69, fixed = TRUE) +
  gwesp(0.69, fixed = TRUE) +
  nodecov("Age") +
  absdiff("Age")

# wanna grab the true params to make the initialisation easier
trueModelEst = ergm(trueModel,
                    control = control.ergm(init = c(-4.409, 0.92919, 1.386, 0.163, -0.489),
                                           MCMC.burnin = 20000, 
                                           MCMC.interval = 20000))

trueModelCoef = coef(trueModelEst)

# the observed statistics for the initial statistics
trueModelStats = summary(trueModel)

# setting up a scheme for the variance of the multivariate normal
tuningConst = 1


# simulate networks given the starting parameters that were either set or previously initialised for the warming phase
# cheating once again with the complete network
tempStats <- simulate( trueModel,                                                     # model formula
                       coef = trueModelCoef,                                          # coefficients
                       output = "stats",                                              # save statistical output (not network objects)
                       basis = trueNet,                                               # starting network
                       nsim = 3000,                                                   # number of saved simulated networks
                       control = control.simulate(MCMC.burnin = 1,                    # various MCMC controls
                                                  MCMC.prop.weights='default',
                                                  MCMC.interval = 2000) )


# calculate the statistics for the model after some burn in
tempStats <- tempStats[2001:3000,]

# calculate the parameter covariance matrix by using the inverse covariance matrix scaled by some tuning constant
propSigma <-  solve(cov(tempStats)) * tuningConst 


# can specify some estimation model
peterModel = peterNet ~ edges +
  gwdegree(0.69, fixed = TRUE) +
  gwesp(0.69, fixed = TRUE) +
  nodecov("Age") +
  absdiff("Age")

# loading just the depleted networks
load(here("Data", "20241031_ergmSqOnlyEntrDepNet.RData"))

## Initialising the network
peterNet = as.network(chosenDepletedNet, directed = FALSE)
peterNet %v% 'Age' = scale(londonGangsAtt$Age)[,1]


# grab the coefficient names
coefNames = names(trueModelCoef)

# now for an endogenous estimation model
ergmSqRes = ergmSquaredSampler(formula = peterModel,
                               naNet = peterNet,
                               initParams = trueModelCoef,
                               initStats = trueModelStats,
                               propSigma = propSigma,
                               iterations = 1500,
                               entrainment = chosenEntrValue,  #gets overwritten by knownMissParams
                               coefNames = coefNames,
                               knownMissParams = c(-0.8472, chosenEntrValue))


# and save the object
save(ergmSqRes, 
     file = here("Output", "20250129_ergmSqEndoXBernD", 
                 paste(format(Sys.Date(), "%Y%m%d"), "_ergmSqOut_entrVal", chosenEntrValue,".RData", sep = "")))