# packages
library(here)
library(sna)
library(stringr)

# load in the BayesALAAM
source(here("Scripts", "bayesNetworkRoutines", "MultivarALAAMalt.R"))


# the functions I'll need to simulate an ALAAM would be
# prepALAAM to format the data in the ALAAMobj format
# getALAAMstats on the ALAAMobj to get the statistics vector (statsvec), simulate.alaam does this if no vector is specified so not necessary
# simulate.alaam because this one takes in the ALAAMobj and the statsvec (can be left NULL since the function does this if statsvec input is NULL)


# let's first use some dataset for an example
# using the R mtcars default and a randomly generated network
# the 'am' variable in mtcars is a binary variable where 0 = automatic transmission, 1 = manual transmission
outcomeVar = mtcars$am

# reverse engineer sample size to get an adjmat
n = length(outcomeVar)
adjMat = matrix(data = runif(n*n),
                nrow = n,
                ncol = n)


# random network that'll get symmetrise because UNDIRECTED
diag(adjMat) = 0
adjMat[lower.tri(adjMat)] = t(adjMat)[lower.tri(adjMat)]
isSymmetric(adjMat)

# binarise for a specified density... say 25%
adjMat[adjMat > 0.25] = 0
adjMat[adjMat != 0] = 1

# plot?
gplot(adjMat, gmode = "graph")

# get some covariates
# just going to use the mtcars variables 
# one float/continuous variable(mpg - miles per gallon) and one binary variable (vs - engine shape)
modelCovs = cbind(mtcars$mpg, mtcars$vs)

# standardising the mpg variable so that the estimation is more stable since its scale seems to be messing with the estimation
centeredMpg = (mtcars$mpg - mean(mtcars$mpg))/sd(mtcars$mpg)
modelCovs = cbind(centeredMpg, mtcars$vs)

# prep ALAAM data
alaamObjExample = prepALAAMdata(y = outcomeVar,
                                ADJ = adjMat,
                                covariates = modelCovs,
                                directed = FALSE,
                                useDegree = FALSE,
                                contagion = "simple")


# and simulate the ALAAM
simAlaamStatsExample = simulate.alaam(ALAAMobj = alaamObjExample,
                                      statsvec = NULL,
                                      theta = c(-0.2, 0.4, 0.2, 0.2), # randomly chosen theta values for this example simulation
                                      thinning = 10,
                                      NumIterations = 3000)

# i'm assuming the output is c(intercept/edges, contagion, mpg, vs) statistics counts

# the other way to do it(?)
simAlaamStatsOutcomeExample = simulateALAAM(y = outcomeVar, 
                                            EdgeList = alaamObjExample$EdgeList, 
                                            RowIn = alaamObjExample$RowIn, 
                                            degree = alaamObjExample$degree, 
                                            covariates = modelCovs, 
                                            NumIterations = 3000,
                                            thinning= 10,
                                            theta = c(-0.2, 0.4, 0.2, 0.2),
                                            statsvec = getALAAMstats(alaamObjExample),
                                            returnNet = TRUE,
                                            canchange = c(1:n))
print(simAlaamStatsOutcomeExample$y)

# would y in the simAlaamStatsOutcomeExample be a predicted outcome variable given the theta values?
# proof of concept-wise this seems successful
# now in terms of getting parameter values
# and specifying different parameters, especially structural parameters, to use as covariates
# is an important question


# what strcutural covariates to look at for undirected networks...?
# let's first try degree/activity
degreeCov = rowSums(adjMat)

# plug that into a combined covariate data frame to use as predictors
combinedCovs = cbind(modelCovs, degreeCov)

# remake the alaam object, now with degree ( I know there's an implicit way to make degree a covariate )
alaamObjDegExample = prepALAAMdata(y = outcomeVar,
                                  ADJ = adjMat,
                                  covariates = combinedCovs,
                                  directed = FALSE,
                                  useDegree = FALSE,
                                  contagion = "simple")

# and simulate
simAlaamOutcome = simulateALAAM(y = outcomeVar, 
                                EdgeList = alaamObjDegExample$EdgeList, 
                                RowIn = alaamObjDegExample$RowIn, 
                                degree = alaamObjDegExample$degree, 
                                covariates = combinedCovs, 
                                NumIterations = 3000,
                                thinning= 10,
                                theta = c(-0.2, 0.4, 0.3, 0.2, -0.2),
                                statsvec = getALAAMstats(alaamObjDegExample),
                                returnNet = TRUE,
                                canchange = c(1:n))
print(simAlaamOutcome)

# this proof of concept works for simulating an ALAAM,
# now let's try doing this with a UCINet empirical covert network
# load in the london gangs dataset (Net 6)
load(here("Data", "20231006_missNetsEnMasse.RData"))

# take the specific adjmat
exampleMat = adjMatList[[6]]


### TODO:: check to make sure the function's actually doing what it's supposed to. Ideally figure out a unit test for it
# a function made in May 2022 to translate an adjmat into a line graph, go 2022 me
lineGraph <- function(adjMat, directed = FALSE){
  
  ## lineGraph(adjMat, directed) takes an adjacency matrix and generates its 
  ## corresponding line graph for its edges. 
  ##
  ## Input:
  ## - adjMat:   An adjacency matrix describing a graph. Not sure if missing data 
  ##             will mess with the conversion. 
  ##
  ## - directed: A logical value to indicate whether the graph is directed or not.
  ##
  ## Output:
  ## - The line graph for the input adjacency matrix, stored as another matrix.
  
  ### count the number of edges
  edgeCount = sum(adjMat)
  
  # if undirected, divide in half.
  if (directed == FALSE){
    edgeCount = edgeCount/2
    adjMat[lower.tri(adjMat)] = 0
  }
  
  # make a matrix of indices for edges (basically an edgelist)
  inds = which(adjMat == 1, arr.ind = TRUE)
  
  # order it for readability
  inds = inds[order(inds),]
  
  # make an empty matrix for the edges (represented as nodes)
  lineMat = matrix(data = 0, 
                   nrow = edgeCount,
                   ncol = edgeCount)
  
  # use the edge list to make the line graph matrix
  for(edgeInd in 1:edgeCount){                              # for every edge in the edge list
    tempMat = matrix(data = c(inds %in% inds[edgeInd,]),    # match the edge with the same nodes
                     ncol = 2)    
    tempMat[edgeInd,] = FALSE                               # remove the selected edge in the matches
    lineMat[edgeInd,] = tempMat[,1] + tempMat[,2]           # and put all the edges with the same nodes in the matrix
  }
  
  
  # make a vector of labels for the nodes
  nodeLabel = apply(inds,
                    MARGIN = 1,
                    FUN = toString)
  
  # remove whitespace
  nodeLabel = gsub(" ", "", nodeLabel, fixed = TRUE)
  
  # name the matrix
  rownames(lineMat) = nodeLabel
  
  # return the line graph matrix
  return(lineMat)
}


# grab the line graph
lineTest = lineGraph(adjMat = exampleMat, directed = FALSE)

# but for the siulation, you don't need an outcome vectorspecified beforehand.. so
emptyLineOutcome = rep(0, times = ncol(lineTest))

# get unique ids for dyadic cov
lineTestEdgeListStr = str_split(rownames(lineTest), pattern = ",")

# reformat into a more friendly data frame
lineTestEdgeList = data.frame(sender = as.numeric(unlist(lapply(lineTestEdgeListStr, FUN = function(x){x[1]}))), 
                              receiver = as.numeric(unlist(lapply(lineTestEdgeListStr, FUN = function(x){x[2]}))))

# some covariates, these need to be specifc to the edge... so they need to be formatted as dyadic covariates
centeredAge = (londonGangsAtt$Age - mean(londonGangsAtt$Age))/sd(londonGangsAtt$Age)

# turn into a 'edge' covariate using the indices in the edge list data frame
centAgeEdgeCov = c()
prisonEdgeCov = c()

# just loop it
for(edgeNo in 1:nrow(lineTestEdgeList)){

  sender = lineTestEdgeList[edgeNo, "sender"]
  receiver = lineTestEdgeList[edgeNo, "receiver"]

  centAgeEdgeCov[edgeNo] = centeredAge[sender] + centeredAge[receiver]
  prisonEdgeCov[edgeNo] = londonGangsAtt$Prison[sender] + londonGangsAtt$Prison[receiver]

}

# I was incorrect about being incorrect- I need a 133-long vector because it's an ALAAM.
# # that was incorrect- I want a dyadic covariate for ALL the senders and receivers
# centAgeDyadCov = matrix(0, nrow = nrow(lineTest), ncol = ncol(lineTest))
# prisonDyadCov = matrix(0, nrow = nrow(lineTest), ncol = ncol(lineTest))
# 
# # just loop it
# for(sender in 1:nrow(lineTest)){
# 
#   # select the sender node id
#   sendNode = lineTestEdgeList[sender,"sender"]
#   
#   # take the receiver vector
#   receivers = lineTestEdgeList[,"receiver"]
# 
#   # these dyadcov are dependent on which edges are in the edge list...
#   centAgeDyadCov[sender, ] = centeredAge[sender] + centeredAge[receivers]
#   prisonDyadCov[sender, ] = londonGangsAtt$Prison[sender] + londonGangsAtt$Prison[receivers]
#   
# }

# since Prison is a binary variable, the dyad cov can be 0, 1 or 2.
# let's keep it to mutual ties for now
prisonEdgeCov[prisonEdgeCov == 2] = 1


# taking prison out for now since it's not that much heterogeneity


# covs
lineAlaamCovs = cbind(centAgeEdgeCov)

# ALAAM stuff
linePrep = prepALAAMdata(y = emptyLineOutcome,
                         ADJ = lineTest,
                         covariates = lineAlaamCovs,
                         directed = FALSE,
                         useDegree = FALSE,
                         contagion = "simple")

# simulate
simAlaamOutcome = simulateALAAM(y = emptyLineOutcome, 
                                EdgeList = linePrep$EdgeList, 
                                RowIn = linePrep$RowIn, 
                                degree = linePrep$degree, 
                                covariates = linePrep$covariates, 
                                NumIterations = 3000,
                                thinning= 10,
                                theta = c(-2, 0.1, 0.1),
                                statsvec = getALAAMstats(linePrep),  # getALAAMstats doesn't work when the y is empty (no obs stats)
                                returnNet = TRUE,
                                canchange = c(1:nrow(lineTest)))
print(simAlaamOutcome)

# check how many edges are sampled
table(simAlaamOutcome$y)


# now do the same for multiple values for contagion
zeroContSim = simulateALAAM(y = emptyLineOutcome, 
                            EdgeList = linePrep$EdgeList, 
                            RowIn = linePrep$RowIn, 
                            degree = linePrep$degree, 
                            covariates = linePrep$covariates, 
                            NumIterations = 1000,
                            thinning= 5,
                            theta = c(-2, 0.1, 0.1),
                            statsvec = getALAAMstats(linePrep),  # getALAAMstats doesn't work when the y is empty (no obs stats)
                            returnNet = TRUE,
                            canchange = c(1:nrow(lineTest)))
print(zeroContSim)

# check how many edges are sampled
table(zeroContSim$y)

# a small contagion value
smallContSim = simulateALAAM(y = emptyLineOutcome, 
                            EdgeList = linePrep$EdgeList, 
                            RowIn = linePrep$RowIn, 
                            degree = linePrep$degree, 
                            covariates = linePrep$covariates, 
                            NumIterations = 1000,
                            thinning= 5,
                            theta = c(-2, 0.3, 0.1),
                            statsvec = getALAAMstats(linePrep),  # getALAAMstats doesn't work when the y is empty (no obs stats)
                            returnNet = TRUE,
                            canchange = c(1:nrow(lineTest)))
print(smallContSim)

# check how many edges are sampled
table(smallContSim$y)

# a larger contagion value
largeContSim = simulateALAAM(y = emptyLineOutcome, 
                             EdgeList = linePrep$EdgeList, 
                             RowIn = linePrep$RowIn, 
                             degree = linePrep$degree, 
                             covariates = linePrep$covariates, 
                             NumIterations = 1000,
                             thinning= 5,
                             theta = c(-2, 0.6, 0.1),
                             statsvec = getALAAMstats(linePrep),  # getALAAMstats doesn't work when the y is empty (no obs stats)
                             returnNet = TRUE,
                             canchange = c(1:nrow(lineTest)))
print(largeContSim)

# check how many edges are sampled
table(largeContSim$y)


# code to get their respective edge lists, the contagion parameter is real finicky
zeroContEdgeList = lineTestEdgeList[zeroContSim$y == 1,]
smallContEdgeList = lineTestEdgeList[smallContSim$y == 1,]
largeContEdgeList = lineTestEdgeList[largeContSim$y == 1,]

# edge list to adjMat
emptyExampleMat = matrix(data = 0, nrow = nrow(exampleMat), ncol = ncol(exampleMat))
zeroContMat = emptyExampleMat
smallContMat = emptyExampleMat
largeContMat = emptyExampleMat

# loops one at a time
for(edgeNo in 1:nrow(zeroContEdgeList)){
  
  # grab the sender
  sender = zeroContEdgeList[edgeNo, 1]
  
  # grab the rceiver
  receiver = zeroContEdgeList[edgeNo, 2]
  
  # fill in the empty matrix
  zeroContMat[sender, receiver] = 1
  
  # do the same for the other end because undirected
  zeroContMat[receiver, sender] = 1
}


# but anyways
for(edgeNo in 1:nrow(smallContEdgeList)){
  
  # grab the sender
  sender = smallContEdgeList[edgeNo, 1]
  
  # grab the rceiver
  receiver = smallContEdgeList[edgeNo, 2]
  
  # fill in the empty matrix
  smallContMat[sender, receiver] = 1
  
  # do the same for the other end because undirected
  smallContMat[receiver, sender] = 1
}

for(edgeNo in 1:nrow(largeContEdgeList)){
  
  # grab the sender
  sender = largeContEdgeList[edgeNo, 1]
  
  # grab the rceiver
  receiver = largeContEdgeList[edgeNo, 2]
  
  # fill in the empty matrix
  largeContMat[sender, receiver] = 1
  
  # do the same for the other end because undirected
  largeContMat[receiver, sender] = 1
}

# and plot
gplot(zeroContMat, gmode = "graph")

gplot(smallContMat, gmode = "graph")

gplot(largeContMat, gmode = "graph")

### two possible avenues with the 'sampled' adjacency matrices
## Turn into network objects
zeroContNet = as.network(zeroContMat, directed = FALSE)
smallContNet = as.network(smallContMat, direced = FALSE)
largeContNet = as.network(largeContMat, directed = FALSE)


## calculating various metrics
getMetrics = function(net, directed = FALSE){
  
  # fork for directed/undirecited
  mode = "digraph"
  
  if(directed == FALSE){
    mode = "graph"
  } 
  
  # density
  density = gden(net, mode = mode)
  
  # clustering coefficient
  clustCoeff = gtrans(net, mode = mode)
  
  # centralisation
  centralisation = centralization(net, FUN = degree, mode = mode)
  
  # average geodesic
  avgGeod = mean(geodist(net, inf.replace= NA)$gdist, na.rm = T)
  
  # diameter
  diameter = max(geodist(net, inf.replace= NA)$gdist, na.rm = T)
  
  # avg degree
  avgDegree = mean(sna::degree(net, gmode = mode))
  
  # avg betweenness
  avgBetw = mean(betweenness(net, gmode = mode))
  
  # slap them all in a list
  metricList = list(
    density = density,
    clustCoeff = clustCoeff,
    centralisation = centralisation,
    avgGeod = avgGeod,
    diameter = diameter,
    avgDegree = avgDegree,
    avgBetw = avgBetw
  )
  
  # reutrn the metric list
  return(metricList)
}

# apply to the networks
zeroContMetrics = getMetrics(zeroContNet)
smallContMetrics = getMetrics(smallContNet)
largeContMetrics = getMetrics(largeContNet)

# NOTE: if I wanted to compare to independently sampled dyads, I would need to match the density... right?




## re-estimation of some ERGMs and compare to the baseline
# can grab code from 




### the 'extra credit' for seeing how the markov assumption in dyads translate to the contagion parameter in the ALAAM
# we first need to get an outcome variable (the sampled edge indicator S_l)
# ... first simulate a sampling matrix so that the outcome variable (the sampling indicator) isn't all 1s
# I can use a missNet generator to generate a sampling matrix
# i.e., 75% sampled = 25% missing
source(here("Scripts", "20230811_missNetFunctions.R"))

## Start with a completely independnet example just for computational ease (just sorting out the code first)
randomSampleInd = initMissAdj(n = nrow(exampleMat), propMiss = 0.75, mode = "graph")

# multiply the sampling indicator with the true matrix t get the sampled adjmat
sampledMat = exampleMat * randomSampleInd

## finding which edges remain (in a very roundabout way)
# turn the sampled mat into a line graph
sampledLines = lineGraph(sampledMat, directed = FALSE)

# grab the edge list
sampledEdgeListStr = str_split(rownames(sampledLines), pattern = ",")

# foramat the edge list
sampledEdgeList = data.frame(sender = as.numeric(unlist(lapply(sampledEdgeListStr, FUN = function(x){x[1]}))), 
                             receiver = as.numeric(unlist(lapply(sampledEdgeListStr, FUN = function(x){x[2]}))))

# find which edges overlap
# can't access internet right now so... i'm combining both edge lists and seeing if there are any duplicates
# since the non-duplicates would be the complement to the sampled edge list, they are 0s in the sampled edge indicator (S_L)

combinedEdgeListStr = c(rownames(lineTest), rownames(sampledLines))

# find which ones only appear once
unsampledEdgeStr = names(which(table(combinedEdgeListStr) == 1))

# then do some testing to ensure the two sets are complementary
combinedEdgeTest = c(rownames(sampledLines), unsampledEdgeStr, rownames(lineTest))

# should be 0
sum(table(combinedEdgeTest) != 2)

# therefore the two sets are complements
# so now let's define an outcome variable to be estimated
sampledOutcomeVar = rep(1, nrow(lineTest))
sampledOutcomeVar[rownames(lineTest) %in% unsampledEdgeStr] = 0

# this is then technically ready to estimate
### Estimation model here.






### Simulation of a snowball sample mechanism
# with the london gangs as the example
# let's say we have a wave size of between 15 for each wave
waveSize = 15

possibleNodes = 1:nrow(exampleMat)

# sample the seed set
seedSet = sample(x = possibleNodes, size = waveSize, replace = FALSE)

# initialise an empty matrix
snowMat = matrix(data = 0, nrow = length(possibleNodes), ncol = length(possibleNodes))

# then sample the ones that have been sampled
snowMat[seedSet,] = exampleMat[seedSet,]

# symmetrising because undirected
snowMat[, seedSet] = exampleMat[, seedSet]

# then prioritise sampling nodes that are incident on the already sampled nodes
seedEdgeList = which(snowMat==1, arr.ind = TRUE)

# then check which columns are incident to the seed set
possibleWave1 = setdiff(as.numeric(names(table(seedEdgeList))), seedSet)

# so sample wave 1
wave1Set = sample(x = possibleWave1, size = waveSize, replace = FALSE)

# sample wave 1
snowMat[wave1Set,] = exampleMat[wave1Set,]
snowMat[,wave1Set] = exampleMat[,wave1Set]

# check which ones can be sampled again
wave1EdgeList = which(snowMat==1, arr.ind = TRUE)

# then check which columns are incident to the seed set
possibleWave2 = setdiff(as.numeric(names(table(wave1EdgeList))), c(seedSet, wave1Set))

# keep going
# so sample wave 2
wave2Set = sample(x = possibleWave2, size = waveSize, replace = FALSE)

# sample wave 1
snowMat[wave2Set,] = exampleMat[wave2Set,]
snowMat[,wave2Set] = exampleMat[,wave2Set]

# check which ones can be sampled again
wave2EdgeList = which(snowMat==1, arr.ind = TRUE)

# then check which columns are incident to the seed set
possibleWave3 = setdiff(as.numeric(names(table(wave2EdgeList))), c(seedSet, wave1Set, wave2Set))

# so sample wave 3, no randomised sampling here because there are not enough nodes
wave3Set = possibleWave3

# sample wave 1
snowMat[wave3Set,] = exampleMat[wave3Set,]
snowMat[,wave3Set] = exampleMat[,wave3Set]

# check which ones can be sampled again
wave3EdgeList = which(snowMat==1, arr.ind = TRUE)
