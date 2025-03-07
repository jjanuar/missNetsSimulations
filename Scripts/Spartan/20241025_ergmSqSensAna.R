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

# get the 'true' values
trueModel = trueNet ~ edges +
  gwdegree(0.69, fixed = TRUE) +
  gwesp(0.69, fixed = TRUE) +
  nodecov("Age") +
  absdiff("Age")

trueStats = as.numeric(summary(trueModel))

print(trueStats)

# plot (and save coordinates)
coord = gplot(trueAdj, gmode = "graph", xlab = "True network")

# some (randomly specified) proportion of missingness
propMiss = 0.3

## generate the initial missingness matrix
initAdj = initMissAdj(n = nrow(trueAdj), propMiss = propMiss)

# turn into a network object
initMissNet = network::as.network(initAdj, directed = FALSE)
initMissNet %v% 'degree' = rowSums(trueAdj)

# setting coefficients for the ERGM to make the MNAR mechanism
missCoefs = c(0, 0, 0, 0.2)  # everything set to positive is dangerous for mcmc explosions, but fixed density helps.

## specify a model for the missingness
# Note: 0.69 is chosen for the decay value because that's what log(2) is. Rounding lets the diagnostic plots be more interpretable

## this specification was chosen for a couple of reasons:
# a) social circuit assumptions are difficult to justify substantively for missingness
#    but if the point I want to make is that 'dependency maters for missingness', they will suffice
#    because markov model parameters are known to be difficult to work with
# b) dyadCov for the true network (an entrainment effect) to weigh the probability of missingness
#    on the presence of a tie in the true network (either positive or negative weights can be justiied)
# removed the degree covariate since I want a simpler missingness model with only entrainment as the effect conditional on the network

chosenModel = 
  initMissNet ~ edges + 
  gwdegree(0.69, fixed = TRUE) +
  gwesp(0.69, fixed = TRUE) +
  dyadcov(trueAdj) 

# and now we simulate
missErgmSimNets = simMissNet(model = chosenModel, coef = missCoefs)

# and check diagnostics
missSimDiag(simObj = missErgmSimNets)

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
peterModel = peterAdjmat ~ edges + gwdegree(0.69, fixed = TRUE) + gwesp(0.69, fixed = TRUE) + nodecov("Age") + absdiff("Age")

# loading just the depleted networks
load(here("Data", "20241031_ergmSqOnlyEntrDepNet.RData"))

## Initialising the network
peterAdjmat = as.network(chosenDepletedNet, directed = FALSE)

# normalise the age so the variance doesn't get too big
# add the age attribute
peterAdjmat %v% 'Age' = scale(londonGangsAtt$Age)[,1]

# try it out
output = ergmSquaredSampler(formula = peterModel, 
                            adjMat = peterAdjmat,
                            initParams = initTrueTheta,
                            initStats = initTrueStats,
                            propSigma = propSigma,
                            iterations = 5000,
                            entrainment = chosenEntrValue,
                            coefNames = coefNames)

# and save the object
save(output, 
     file = here("Output", "20241025_ergmSqSensAna", 
                 paste(format(Sys.Date(), "%Y%m%d"), "_ergmSqOut_entrVal", chosenEntrValue,".RData", sep = "")))