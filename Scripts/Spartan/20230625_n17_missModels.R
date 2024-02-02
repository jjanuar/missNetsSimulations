### Parallelisation/indexing info
# grabbing command arguments
command_args <- commandArgs(trailingOnly = TRUE)

# index for which missingness model is used
missModelIndex = as.numeric(command_args[1])

# array index for number of runs
runIndex = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# printing specified indices  
missModels = c("Indep", "UM", "CM", "ERGM")
paste("The model chosen is", missModels[missModelIndex])
paste("The run is", runIndex)

### script to generate D, deplete the network, and then re-estimate with the missingness
## Loading packages
library(ergm)
library(network)
library(sna)
library(here)


## Loading data
# as taken from UCINET (17 Nov Greek Bombings)
greeceBombs = as.matrix(read.csv(file = here("Data", "17Nov Greece Bombing", "17NOV_GREECE_19751984.csv"), header = T, row.names = "X"))

# symmetrising column and row names because isSymmetric checks this
colnames(greeceBombs) = rownames(greeceBombs)


## Function to deplete independently
# degrade net function from previous script
indepMiss <- function(adjMat, propMiss, directed = FALSE){
  
  ## indepMiss(adjMat, propMiss, directed) takes a graph and replaces the observed graph ties with missing 
  ## values(NA). This function does not handle actor non-response as it only works off the observed
  ## network ties. Works with undirected networks.
  ##
  ## Input:
  ## - adjMat:    An adjacency matrix describing a graph. There shouldn't
  ##             be any missing values in this adjacency matrix yet.
  ##             Should still work if there were though. I think.
  ## - propMiss: A numeric value to indicate the proportion of missingness
  ##             that will be imposed on the observed network.
  ##             Is bounded between 0 and 1.
  ## - directed: A logical value to indicate if it's a directed or
  ##             undirected network. Default is set to undirected.
  ##
  ## Output: 
  ## - A degraded network with missing ties where observed ties once were.
  ##   Note that it only handles item non-response for now. Given the inputs,
  ##   the missing values would be ~'propMiss'% of the observed tie variables in 'graph'.
  
  ## spit out an error if the directed argument is misspecified
  if(directed != TRUE & !isSymmetric(adjMat)){
    stop("The undirected network doesn't have a symmetric matrix")
  }
  
  ## spit out an error if the proportion of missingness exceeds bounds
  if(propMiss > 1 | propMiss < 0){
    stop("The proportion of missingness exceeds bounds")
  }
  
  # grab the number of nodes
  n = nrow(adjMat)
  
  # make an n x n matrix containing randomly generated values between 0 and 1
  missThresh = matrix(data = runif(n^2, min = 0, max = 1),
                      nrow = n,
                      ncol = n)
  
  # set diagonal to 0 for the dircted network
  diag(missThresh) = 0
  
  # branch for undirected network
  if( directed == FALSE ){
    
    # make a matrix
    missThresh = matrix(data = 0, nrow = n, ncol = n)
    
    # take only one triangle
    missThresh[upper.tri(missThresh)] = runif((n*(n-1))/2)
    
    # and symmetrise
    missThresh[lower.tri(missThresh)] = t(missThresh)[lower.tri(missThresh)]
  }
  
  # index which of the ties are going to be missing
  missTies = missThresh <= propMiss
  
  # copy the graph to punch holes in it
  degradedGraph = adjMat
  
  # punch holes
  degradedGraph[missTies] = NA
  
  # diagonal is always fixed to 0, this is an artefact of the way missingness is generated
  diag(degradedGraph) = 0
  
  # return the degraded graph
  return(degradedGraph)
}

## Loading simulated missingness
load(here("Data", "20230622_simulated_missing_nets.RData"))

## Choosing a D
# a loop to draw some (missingness) networks with over 25% missing
# initialise
goodCount = 0

# a while loop becuase i'm spicy
while(goodCount < 4){
  
  # randomly sample an index
  index = sample(1:1000, 1)
  
  # differently for the independent model
  # use the function
  indepMissNet = indepMiss(adjMat = greeceBombs, propMiss = 0.25, directed = FALSE)
  
  # take a network out
  chosenIndep = 1 * is.na(indepMissNet)
  chosenMissNetUM = simMissNetsUM[index]
  chosenMissNetCM = simMissNetsCM[index]
  chosenMissNetERGM = simMissNetsERGM[index]
  
  # grab the edge counts
  edgeCounts = c(sum(chosenIndep)/2, 
                 network.edgecount(chosenMissNetUM[[1]]), 
                 network.edgecount(chosenMissNetCM[[1]]), 
                 network.edgecount(chosenMissNetERGM[[1]]))
  
  # check if in range
  over40 = edgeCounts >= 40
  under49 = edgeCounts < 49
  edgeRange = over40 & under49
  
  # and average it
  goodCount = sum(edgeRange)
}

## Depleting, ignoring Todd for now.
# turn them all into matrix objects
chosenMissMatUM = as.matrix(chosenMissNetUM[[1]])
chosenMissMatCM = as.matrix(chosenMissNetCM[[1]])
chosenMissMatERGM = as.matrix(chosenMissNetERGM[[1]])


# copy over the adjacency matrix
greeceBombsMissIndep = greeceBombs
greeceBombsMissUM = greeceBombs
greeceBombsMissCM = greeceBombs
greeceBombsMissERGM = greeceBombs

## CURRENTLY SET TO Todd
# Now we're actually depleting the network
greeceBombsMissIndep[chosenIndep == 1] = 0
greeceBombsMissUM[chosenMissMatUM == 1] = 0
greeceBombsMissCM[chosenMissMatCM == 1] = 0
greeceBombsMissERGM[chosenMissMatERGM == 1] = 0

# save the depleted network
depletedNetworkList = list(greeceBombsMissIndep, greeceBombsMissUM, greeceBombsMissCM, greeceBombsMissERGM)

# chosen depleted network
chosenDepletedNetwork = depletedNetworkList[[missModelIndex]]

# turn into networks again
bombNetMissIndep = as.network(greeceBombsMissIndep, directed = FALSE)
bombNetMissUM = as.network(greeceBombsMissUM, directed = FALSE)
bombNetMissCM = as.network(greeceBombsMissCM, directed = FALSE)
bombNetMissERGM = as.network(greeceBombsMissERGM, directed = FALSE)


## Re-estimating
# models to be reestimated
bombNetMissIndepModel = bombNetMissIndep ~ edges + gwdegree(decay = log(2), fixed = TRUE) + gwesp(decay = log(2), fixed = TRUE) 
bombNetMissUMModel = bombNetMissUM ~ edges + gwdegree(decay = log(2), fixed = TRUE) + gwesp(decay = log(2), fixed = TRUE) 
bombNetMissCMModel = bombNetMissCM ~ edges + gwdegree(decay = log(2), fixed = TRUE) + gwesp(decay = log(2), fixed = TRUE) 
bombNetMissERGMModel = bombNetMissERGM ~ edges + gwdegree(decay = log(2), fixed = TRUE) + gwesp(decay = log(2), fixed = TRUE) 


# put them all in a list
bombNetMissModelList = list(bombNetMissIndepModel, bombNetMissUMModel, bombNetMissCMModel, bombNetMissERGMModel)
initParaList = list(c(-5.6, 1.039, 2.0331),
                    c(-3.596, 0.2393, 1.123),
                    c(-4.02, 0.419, 1.327),
                    c(-3.104, 0.0114, 0.9232))

# chosen models and initial parameters
chosenBombNetMissModel = bombNetMissModelList[[missModelIndex]]
chosenInitPara = initParaList[[missModelIndex]]

# print to make sure
paste("The index is", index)
paste("The chosen Miss model is", chosenBombNetMissModel)
cat("The initial parameters are", chosenInitPara ,"\n")

# reeeeeestimations
# I'm expecting some of these to fail.
# initial values taken from manual ergms that worked
chosenDepletedMiss = ergm(chosenBombNetMissModel,
                          control = control.ergm(init = chosenInitPara, MCMC.interval = 5024))

## Saving
save(chosenDepletedMiss, chosenDepletedNetwork, chosenIndep, index, file = here("Output", "20230625_n17_missModels", paste("20230625_nov17_missModels_Todd_",missModels[missModelIndex], "_", runIndex, ".RData", sep = "")))