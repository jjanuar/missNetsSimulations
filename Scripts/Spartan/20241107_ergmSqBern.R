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


# setting up a scheme for the variance of the multivariate normal
tuningConst = 1

# bernoulli example
bernModel = trueNet ~ edges

bernEdge = as.numeric(summary(bernModel))
print(bernEdge)

bernRes = ergm(bernModel,
               control = control.ergm(init = c(-2.27823),
                                      MCMC.burnin = 20000, 
                                      MCMC.interval = 20000))

# grab the initial statistic
bernInitTheta = bernRes$coefficients
bernNames = names(coef(bernRes))

# getting the proposal steps for only the edges
bernTempStats <- simulate( bernModel,                                                     # model formula
                       coef = bernRes$coefficients,                                   # coefficients
                       output = "stats",                                              # save statistical output (not network objects)
                       basis = trueNet,                                               # starting network
                       nsim = 3000,                                                   # number of saved simulated networks
                       control = control.simulate(MCMC.burnin = 1,                    # various MCMC controls
                                                  MCMC.prop.weights='default',
                                                  MCMC.interval = 2000) )

# calculate the statistics for the model after some burn in
burnedBernStats <- bernTempStats[2001:3000,]

bernPropSigma = sd(burnedBernStats)^(-1) * tuningConst

# can specify some estimation model
peterBernModel = peterAdjmat ~ edges

# loading just the depleted networks
load(here("Data", "20241031_ergmSqOnlyEntrDepNet.RData"))

## Initialising the network
peterAdjmat = as.network(chosenDepletedNet, directed = FALSE)

# try it out, not using the function because rmvnorm needs more than one dimension. using rnorm instead

formula = peterBernModel
adjMat = peterAdjmat
initParams = bernInitTheta
initStats = bernEdge
propSigma = bernPropSigma
iterations = 1500
entrainment = chosenEntrValue
coefNames = bernNames


# get some values
numberPara = length(initParams)

## initialising some storage objects
# the sampled parameters
sampledThetas = matrix(data = NA, nrow = iterations, ncol = numberPara) 
colnames(sampledThetas) = coefNames

# auxiliary and imputed network statistics
auxNetStatMat = matrix(data = NA, nrow = iterations, ncol = numberPara)
colnames(auxNetStatMat) = coefNames

impNetStatMat = matrix(data = NA, nrow = iterations, ncol = numberPara)
colnames(impNetStatMat) = c(coefNames)

# parameters with entrainment adjustment to the density
# psi = matrix(data = NA, nrow = iterations, ncol = numberPara)
# colnames(psi) = coefNames

psi = matrix(data = NA, nrow = iterations, ncol = numberPara)
colnames(psi) = c("edges")

# a list for the imputed networks
impNetList = list()

## sequentially printing some iterations to make sure the sampler's progressing
printIter = floor(seq(from = 1, to = iterations, length.out = 5))


# get the list of free dyads
# need the missingness indicator so we can work backwards from the adjmat
missIndMat = (as.matrix(is.na(adjMat)) * 1)
missTiesEdgeList = as.edgelist(as.network(missIndMat, directed=FALSE), n = nrow(missIndMat))


## Starting the sampler
for(iter in 1:iterations){
  
  # if current iteration is in any of the specified printing iterations
  if(iter %in% printIter){
    
    # print something out to show that the sampler's progressing
    message(paste("Current sampler iteration is", iter))
  }
  
  # set something up for the first iteration
  if(iter == 1){
    nextInitTheta = initParams
    nextInitStats = initStats
    
  } else {
    # and for the rest of the iterations
    nextInitTheta = sampledThetas[iter-1, 1]
    nextInitStats = impNetStatMat[iter-1, 1]
  }
  
  ## Generating a proposal theta and auxiliary network
  # auxiliary parameters drawn from a multivariate normal distribution
  auxTheta = rnorm(n = 1, mean = nextInitTheta, sd = propSigma)
  
  # generate a network using the auxiliary parameters
  auxNetStats = simulate( object = formula, 
                          coef = auxTheta,
                          output = "stats",
                          basis = adjMat,
                          nsim = 1,
                          control = control.simulate(MCMC.burnin = 20000,
                                                     MCMC.interval = 2000))
  
  # saving its statistics
  auxNetStatMat[iter, ] = auxNetStats
  auxStats = as.numeric(auxNetStats)
  
  ## Calculate the acceptance ratio
  acceptRatio = (auxTheta - nextInitTheta) %*% (nextInitStats - auxStats)
  
  ## parameter swap
  if(log(runif(1)) < acceptRatio){
    
    # swap parameters if acceptance ratio is larger than random
    sampledThetas[iter, ] = auxTheta
  } else {
    
    # keep the previous params
    sampledThetas[iter, ] = nextInitTheta 
  }
  
 
  
  ## Generation of imputed network and statistics
  # use the most recent thetas with the MNAR model to impute the data
  psi[iter, ] = sampledThetas[iter, ]

  # add the entrainment parameter from the missingness model
  psi[iter, 1] = sampledThetas[iter, 1] + entrainment
  
  # formula for the bernoulli model one with the edge covariate 
  impFormula = peterAdjmat ~ edges
  
  # trying out a different way to do it
  # # adding an edge covariate for the missing tie variables
  # impFormula = peterAdjmat ~ edges + edgecov(is.na(peterAdjmat))
  # 
  # # then adding in the adjustment as an additional parameter instead of adjusting the density
  # psi[iter, ] = cbind(sampledThetas[iter,], entrainment)

  # get the conditional distribution of the missing tie variables given the specified model
  impNets = simulate(object = impFormula,
                     coef = psi[iter, ],
                     output = "network",
                     basis = adjMat,
                     nsim = 1,
                     constraints=~fixallbut(missTiesEdgeList),
                     control = control.simulate(MCMC.burnin = 20000,
                                                MCMC.interval = 2000))
  
  # save the imputed networks
  impNetList[[iter]] = impNets
  
  # get its stats
  impNetStats = attributes(impNets)$stats
  
  # save stuff in the structures specified beforehand
  impNetStatMat[iter,] = impNetStats
  
}



# put all the output in a list
output = list(
  sampledThetas = sampledThetas,
  impNetList = impNetList,
  impNetStatMat = impNetStatMat,
  auxNetStatMat = auxNetStatMat)

# quick plot
plot(ts(impNetStatMat[,1]), ylab = 'edges', main = 'Imputed network statistics')
abline(h = 133, col = 'red')

plot(ts(sampledThetas[,1]), ylab = "edges", main = "Parameter estimate")
abline(h = -2.28, col = 'red')

# and save the object
save(output, 
     file = here("Output", "20241107_ergmSqBern", 
                 paste(format(Sys.Date(), "%Y%m%d"), "_ergmSqOut_entrVal", chosenEntrValue,".RData", sep = "")))