### A Spartan-friendly version of a script to run the ERGM squared sampler on a single depleted network
## Including the log-penalty


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

# get the 'true' values
trueModel = trueNet ~ edges +
  gwdegree(0.69, fixed = TRUE) +
  gwesp(0.69, fixed = TRUE) +
  nodecov("Age") +
  absdiff("Age")

trueStats = as.numeric(summary(trueModel))

print(trueStats)



# for future purposes, the coefs are c(-4.409, 0.92919, 1.386, 0.163, -0.489)
# cheat by using the natural mle parameters of the complete network to get initial parameters
trueRes = ergm(trueModel,
               control = control.ergm(init = c(-4.409, 0.92919, 1.386, 0.163, -0.489),
                                      MCMC.burnin = 20000, 
                                      MCMC.interval = 20000))
initTrueTheta = trueRes$coefficients

initTrueStats = summary(trueModel)

# coefficient names
coefNames = names(coef(trueRes))

# setting up a scheme for the variance of the multivariate normal
tuningConst = 1

# simulate networks given the starting parameters that were either set or previously initialised for the warming phase
# cheating once again with the complete network
tempStats <- simulate( trueModel,                                                     # model formula
                       coef = initTrueTheta,                                           # coefficients
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
peterModel = peterNet ~ edges + gwdegree(0.69, fixed = TRUE) + gwesp(0.69, fixed = TRUE) + nodecov("Age") + absdiff("Age")

# loading just the depleted networks
load(here("Data", "20241212_chosenMissInd.RData"))

## Initialising the network
peterAdjmat =  degradeAdj(trueNet = trueAdj,
                          missAdj = as.matrix(chosenMissInd))

peterNet = as.network(peterAdjmat, directed = FALSE)

# normalise the age so the variance doesn't get too big
# add the age attribute
peterNet %v% 'Age' = scale(londonGangsAtt$Age)[,1]


# add the missingness model
missModel = "edges + dyadcov(noNaPropSubstepNet)"

# try it out
output = chaosErgmSquaredSampler(formula = peterModel, 
                                 naNet = peterNet,
                                 initParams = initTrueTheta,
                                 initStats = initTrueStats,
                                 propSigma = propSigma,
                                 iterations = 5000,
                                 entrainment = chosenEntrValue,
                                 coefNames = coefNames,
                                 missModel = missModel)

# and save the object
save(output, 
     file = here("Output", "20241220_chaosErgmSq", 
                 paste(format(Sys.Date(), "%Y%m%d"), "_chaosErgmSqOut_entrVal", chosenEntrValue,".RData", sep = "")))