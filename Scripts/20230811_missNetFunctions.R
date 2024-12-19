## R file containing functions for the simulation of missingness mechanisms with different models
# list of required packages (that have dependencies within them)
requiredPackages = c("testthat", "ergm", "network", "sna")

# install packages if any of the required packages aren't installed
#install.packages(setdiff(requiredPackages, rownames(installed.packages())))

## function to initialise missingness matrix
initMissAdj = function(n, propMiss, numSims = 1, mode = "graph", ...){
  
  
  ## initMissAdj(n, propMiss, numSims, ...) takes a number of nodes, a proportion of missingness, 
  ## a number of simulated networks, and any additional arguments supported by sna::rgnm() to 
  ## generate a uniform network conditioned on the density of the graph.
  ##
  ## Input:
  ## - n:        An integer for the number of nodes for the missingness indicator matrix.
  ##
  ## - propMiss: A numeric value to indicate the proportion of missingness to calculate
  ##             the 'density' of the missingness matrix
  ##             Is bounded between 0 and 1.
  ##
  ## - numSims:  An integer for the number of simulated networks to be generated.
  ##             The default is set to 1 for the function to return a matrix,
  ##             If set to be above one, it will return a numSims x n x n array.
  ##
  ## - mode:     A string to be used for the sna::rgnm() for the type of graph to be generated.
  ##             Default is set to 'graph' for an undirected graph. Can be changed to 'digraph'.
  ##             For directed graphs.
  ##
  ## Output: 
  ## - An n x n matrix (or numSims x n x n array) containing a completely uniform (Bernoulli) graph
  ##   with an exact density (propMiss). Can also be used for independent missingness.
  
  # needs the sna package
  require(sna)
  
  # assume the sna package is read in
  initAdj = sna::rgnm(n = numSims,                       # number of simulated networks
                      nv = n,                            # number of nodes
                      m = round(propMiss*((n*(n-1))/2)), # number of dyads (for density) 
                      mode = mode)                       # type of graph
  
  # return the simulated object
  return(initAdj)
}


### simulation function
simMissNet = function(model, coef, nsim = 300){
  
  ## simMissNet(model, coef) takes a model and a set of coefficients for the model
  ## and plugs it in the simulate.ergm() function to return 1000 (can be modified..) networks
  ## with fixed density (constrained number of edges).
  ##
  ## Input:
  ## - model: A statnet model object. The outcome has to be a network object, the 
  ##          predictors need to be ergm.terms.
  ##
  ## - coef:  A vector of coefficients needed to simulate the networks on. These
  ##          should be chosen corresponding to the local endogenous assumptions made
  ##          about the missingness
  ##
  ## Output: 
  ## - A list of 1000 (modifiable...) networks simulated with the specified model.
  
  # needs ergm or wouldn't work
  require(ergm)
  
  # generate the simulated networks 
  simMissNets <- simulate(object = model, # put in the model object ,
                          coef = coef,  # I *do* need these. edges can be set to 0 when the constraint is in
                          nsim = nsim,
                          control = control.simulate.formula(MCMC.burnin=20000,MCMC.interval=100), # a burnin to prepare the importance sampler
                          constraints = ~edges)   # this constraint fixes the 'density' (i.e., number of missing tie vars)
  
  # return the simulated networks
  return(simMissNets)
}


### Diagnostic function basically
missSimDiag = function(simObj, modelName = "Full ERGM"){
  
  ## missSimDiag(simObj, modelName) takes a simulated matrix of network statistics from 
  ## the simulate.ergm() function and plots the trace plots of the parameters.
  ##
  ## Input:
  ## - simObj:    A matrix object containing the simulated networks from
  ##              the output of a simulate.ergm() function.
  ##
  ## - modelName: A string to indicate the name of the model that is simulated.
  ##              Its default is a 'Full ERGM' simulation.
  ##
  ## Output: 
  ## - For number of parameters p in the simulated object, p trace plots arranged in 
  ##   2 columns.
  
  # needs ergm
  require(ergm)
  
  # take out the statistics from the simulated object
  simObjStats = attr(simObj, "stats")
  
  # number of parameters
  p = ncol(simObjStats)
  
  # setting margins for the trace plots
  # setting it to 2 columns (arbitrarily)
  # and setting the number of rows to the ceiling of the number of parameters, halved
  par(mfrow = c(ceiling(p/2), 2))
  
  # for loop to plot each trace plot
  for(ind in 1:ncol(simObjStats)){
    
    # trace plots to see how the parameter values are doing
    plot(ts(simObjStats[,ind]), ylab = colnames(simObjStats)[ind])
  }
  
  # setting text for the whole plot
  mtext(text = modelName,
        side = 3,
        line = -2,
        outer = TRUE)
}

### Initialising some latent space
simLatentCoords = function(n, noDim = 2){
  
  ## simLatentCoords(n, noDim) takes a simulated matrix of network statistics from 
  ## the simulate.ergm() function and plots the trace plots of the parameters.
  ##
  ## Input:
  ## - n:     An integer for the number of nodes for the missingness indicator matrix.
  ##
  ## - noDim: An integer for the number of dimensions in the latent space. The default is
  ##          set to 2. Upwards of 3 dimensions (4+) might be difficult to visualise. 
  ## Output: 
  ## - An n x noDim matrix with random standardised normal values for each dimension.
  
  # specify the matrix of coordinates
  latentCoords = matrix(0, nrow = n, ncol = noDim)
  
  # using the base constant letters for axis names
  colnames(latentCoords) = letters[1:noDim]
  
  # randomly generating positions with a normal dist
  for(dim in 1:noDim){
    latentCoords[,dim] = rnorm(n, mean = 0, sd = 1)
  }
  
  # return the latent coordinates
  return(latentCoords)
}


# A general purpose function to degrade 'true' matrices 
degradeAdj = function(trueNet, missAdj, missSave = NA){
  
  ## degradeAdj(trueNet, missNet, missSave) takes a 'true network' and then degrades
  ## the network with the missingness adjacency matrix provided. The missingness can then
  ## be saved as either missing values (NA) or any other value (e.g., 1).
  ##
  ## Input:
  ## - trueNet:  An n x n adjacency matrix to indicate the 'true' network.
  ##
  ## - missAdj:  An n x n adjacency matrix to indicate the missingness/observation of tie variables
  ##             for the provided true network.
  ##
  ## - missSave: Any value (usually integer or NA) to be used as the value missingness is saved as.
  ##             The default is set to NA for missing as per R definitions, but can be any value.
  ##
  ## Output: 
  ## - An n x n matrix containing missing values (however saved/stored) for the provided missingness
  ##   matrix.
  
  # make sure the object types are right
  if( inherits(trueNet, "matrix") != inherits(missAdj, "matrix") ){
    stop("Make sure both objects are matrices")
  }
  
  # copy the graph to punch holes in it
  degradedGraph = trueNet
  
  # punch holes
  degradedGraph[missAdj==1] = missSave
  
  # return the degraded graph
  return(degradedGraph)  
  
}
