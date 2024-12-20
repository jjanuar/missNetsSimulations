## Setting up the functional form of the sampler
ergmSquaredSampler = function(formula, naNet, initParams, initStats, propSigma, iterations, entrainment, coefNames){
  
  ## ergmSquaredSampler(formula, naNet, initParams, ...) takes an ergm formula, an adjacency matrix with missingness, 
  ## some initial values, and an entrainment parameter value to sample generate model estimates under 
  ## a specific not-at-random missingness model with a specific entrainment value.
  ##
  ## Input:
  ## - formula:     An R formula object following the ergm specifications of network ~ ergm terms
  ##
  ## - naNet:       A network object with missing edges. Missing values are ideally set to NA values.
  ##                Currently only supports undirected networks.
  ##
  ## - initParams:  A p-long vector containing the initial parameters for the sampler. 
  ##
  ## - initStats:   A p-long vector containing the initial statistics (mean value/observed parameters)
  ##                to calculate the acceptance ratio in the sampler.
  ##
  ## - propSigma:   A p x p matrix reflecting the proposal covariance matrix. When the true network parameters
  ##                are known, use the inverse covariance matrix scaled by some tuning constant.
  ##
  ## - iterations:  The number of iterations to run the sampler.
  ##
  ## - entrainment: A single float corresponding to the entrainment parameter in the missingness model.
  ##                Can be set to 0 to perform 'typical' miss ergm bayes assuming MAR.
  ##
  ## - coefNames:   A p-long vector of coefficient names in the specified estimation model. 
  ##                Purely visual input, there's likely a way to refine this and get names from the formula object.
  ##
  ## Output: 
  ## - output:      A list with 4 items,
  ##                - sampledThetas is an iterations x p matrix containing the sampled parameter values
  ##                - impNetList is a p-long list containing the imputed networks
  ##                - impNetStatMat is an iterations x p matrix containing the imputed network statistics
  ##                - auxNetStatMat is an iterations x p matrix containing the auxiliary network statistics
  
  # requires the ergm and mvtnorm packages
  require(ergm)
  require(mvtnorm)
  
  # have the adjMat for computation purposes
  adjMat = as.matrix(naNet)
  
  # get some values
  numberPara = length(initParams)
  
  # a node size value
  nodeSize = nrow(adjMat)
  
  ## initialising some storage objects
  # the sampled parameters
  sampledThetas = matrix(data = NA, nrow = iterations, ncol = numberPara) 
  colnames(sampledThetas) = coefNames
  
  # auxiliary and imputed network statistics
  auxNetStatMat = matrix(data = NA, nrow = iterations, ncol = numberPara)
  colnames(auxNetStatMat) = coefNames
  
  impNetStatMat = matrix(data = NA, nrow = iterations, ncol = numberPara)
  colnames(impNetStatMat) = coefNames
  
  # parameters with entrainment adjustment to the density
  psi = matrix(data = NA, nrow = iterations, ncol = numberPara)
  colnames(psi) = coefNames
  
  # a list for the imputed networks
  impNetList = list()
  
  ## sequentially printing some iterations to make sure the sampler's progressing
  printIter = floor(seq(from = 1, to = iterations, length.out = 5))
  
  
  # get the list of free dyads
  # need the missingness indicator so we can work backwards from the adjmat
  missIndMat = (as.matrix(is.na(adjMat)) * 1)
  missTiesEdgeList = as.edgelist(as.network(missIndMat, directed=FALSE), n = nrow(missIndMat))
  
  # calculate the alpha (intercept for D) value from the proportion of missing tie variables
  missDens = sum(missIndMat)/(nodeSize * (nodeSize - 1))
  alpha = log(missDens/(1 - missDens))
  
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
      nextInitTheta = sampledThetas[iter-1, ]
      nextInitStats = impNetStatMat[iter-1, ]
    }
    
    ## Generating a proposal theta and auxiliary network
    # auxiliary parameters drawn from a multivariate normal distribution
    auxTheta = rmvnorm(n = 1, mean = nextInitTheta, sigma = propSigma)
    
    # generate a network using the auxiliary parameters
    auxNetStats = simulate( object = formula, 
                            coef = auxTheta,
                            output = "stats",
                            basis = naNet,
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
    # a log ratio penalty
    logPenalty = log((1 + exp(alpha))/(1 + exp(alpha + entrainment)))
    
    ## Generation of imputed network and statistics
    # use the most recent thetas with the MNAR model to impute the data
    psi[iter, ] = sampledThetas[iter, ]
    
    # add the entrainment parameter from the missingness model
    psi[iter, 1] = sampledThetas[iter, 1] + entrainment + logPenalty
    

    # get the conditional distribution of the missing tie variables given the specified model
    impNets = simulate(object = formula,
                       coef = psi[iter, ],
                       output = "network",
                       basis = naNet,
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
  
  # and return it
  return(output)
}


## Setting up the functional form of the sampler
bernErgmSquaredSampler = function(formula, naNet, initParams, initStats, propSigma, iterations, entrainment, coefNames){
  
  ## bernErgmSquaredSampler(formula, naNet, initParams, ...) takes a bernoulli ergm formula, an adjacency matrix with missingness, 
  ## some initial values, and an entrainment parameter value to sample generate model estimates under 
  ## a specific not-at-random missingness model with a specific entrainment value.
  ##
  ## Input:
  ## - formula:     An R formula object following the ergm specifications of network ~ ergm terms
  ##
  ## - naNet:       A network object with missing edges. Missing values are ideally set to NA values.
  ##                Currently only supports undirected networks.
  ##
  ## - initParams:  A p-long vector containing the initial parameters for the sampler. 
  ##
  ## - initStats:   A p-long vector containing the initial statistics (mean value/observed parameters)
  ##                to calculate the acceptance ratio in the sampler.
  ##
  ## - propSigma:   A p x p matrix reflecting the proposal covariance matrix. When the true network parameters
  ##                are known, use the inverse covariance matrix scaled by some tuning constant.
  ##
  ## - iterations:  The number of iterations to run the sampler.
  ##
  ## - entrainment: A single float corresponding to the entrainment parameter in the missingness model.
  ##                Can be set to 0 to perform 'typical' miss ergm bayes assuming MAR.
  ##
  ## - coefNames:   A p-long vector of coefficient names in the specified estimation model. 
  ##                Purely visual input, there's likely a way to refine this and get names from the formula object.
  ##
  ## Output: 
  ## - output:      A list with 4 items,
  ##                - sampledThetas is an iterations x p matrix containing the sampled parameter values
  ##                - impNetList is a p-long list containing the imputed networks
  ##                - impNetStatMat is an iterations x p matrix containing the imputed network statistics
  ##                - auxNetStatMat is an iterations x p matrix containing the auxiliary network statistics
  
  # requires the ergm package
  require(ergm)
  
  # have the adjMat for computation purposes
  adjMat = as.matrix(naNet)
  
  # a node size value
  nodeSize = nrow(adjMat)
  
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
  colnames(impNetStatMat) = coefNames
  
  # parameters with entrainment adjustment to the density
  psi = matrix(data = NA, nrow = iterations, ncol = numberPara)
  colnames(psi) = coefNames
  
  # a list for the imputed networks
  impNetList = list()
  
  ## sequentially printing some iterations to make sure the sampler's progressing
  printIter = floor(seq(from = 1, to = iterations, length.out = 5))
  
  # get the list of free dyads
  # need the missingness indicator so we can work backwards from the adjmat
  missIndMat = (as.matrix(is.na(adjMat)) * 1)
  missTiesEdgeList = as.edgelist(as.network(missIndMat, directed=FALSE), n = nrow(missIndMat))
  
  # calculate the alpha (intercept for D) value from the proportion of missing tie variables
  missDens = sum(missIndMat)/(nodeSize * (nodeSize - 1))
  alpha = log(missDens/(1 - missDens))
  
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
      nextInitTheta = sampledThetas[iter-1, ]
      nextInitStats = impNetStatMat[iter-1, ]
    }
    
    ## Generating a proposal theta and auxiliary network
    # auxiliary parameters drawn from a multivariate normal distribution
    auxTheta = rnorm(n = 1, mean = nextInitTheta, sd = propSigma)
    
    # generate a network using the auxiliary parameters
    auxNetStats = simulate( object = formula, 
                            coef = auxTheta,
                            output = "stats",
                            basis = naNet,
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
    
    
    # a log ratio penalty
    logPenalty = log((1 + exp(alpha))/(1 + exp(alpha + entrainment)))
    
    ## Generation of imputed network and statistics
    # use the most recent thetas with the MNAR model to impute the data
    psi[iter, ] = sampledThetas[iter, ]
    
    # add the entrainment parameter from the missingness model
    psi[iter, 1] = sampledThetas[iter, 1] + entrainment + logPenalty
    
    # get the conditional distribution of the missing tie variables given the specified model
    impNets = simulate(object = formula,
                       coef = psi[iter, ],
                       output = "network",
                       basis = naNet,
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
    #impNetList = impNetList,
    impNetStatMat = impNetStatMat,
    auxNetStatMat = auxNetStatMat)
  
  # and return it
  return(output) 
}

## endogenous missingness mechanism
chaosErgmSquaredSampler = function(formula, naNet, initParams, initStats, propSigma, iterations, entrainment, coefNames, missModel){
  
  ## chaosErgmSquaredSampler(formula, naNet, initParams, ...) takes an ergm formula, an adjacency matrix with missingness, 
  ## some initial values, and an entrainment parameter value to sample generate model estimates under 
  ## a specific not-at-random missingness model with a specific entrainment value.
  ##
  ## Input:
  ## - formula:     An R formula object following the ergm specifications of network ~ ergm terms
  ##
  ## - naNet:       A network object with missing edges. Missing values are ideally set to NA values.
  ##                Currently only supports undirected networks.
  ##
  ## - initParams:  A p-long vector containing the initial parameters for the sampler. 
  ##
  ## - initStats:   A p-long vector containing the initial statistics (mean value/observed parameters)
  ##                to calculate the acceptance ratio in the sampler.
  ##
  ## - propSigma:   A p x p matrix reflecting the proposal covariance matrix. When the true network parameters
  ##                are known, use the inverse covariance matrix scaled by some tuning constant.
  ##
  ## - iterations:  The number of iterations to run the sampler.
  ##
  ## - entrainment: A single float corresponding to the entrainment parameter in the missingness model.
  ##                Can be set to 0 to perform 'typical' miss ergm bayes assuming MAR.
  ##
  ## - coefNames:   A p-long vector of coefficient names in the specified estimation model. 
  ##                Purely visual input, there's likely a way to refine this and get names from the formula object.
  ##
  ## - missModel:   A string object containing the RHS of a formula object for the missingness model.
  ##
  ## Output: 
  ## - output:      A list with 4 items,
  ##                - sampledThetas is an iterations x p matrix containing the sampled parameter values
  ##                - impNetList is a p-long list containing the imputed networks
  ##                - impNetStatMat is an iterations x p matrix containing the imputed network statistics
  ##                - auxNetStatMat is an iterations x p matrix containing the auxiliary network statistics
  
  # requires the ergm and mvtnorm packages
  require(ergm)
  require(mvtnorm)
  
  # have the adjMat for computation purposes
  adjMat = as.matrix(naNet)
  
  # a node size value
  nodeSize = nrow(adjMat)
  
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
  colnames(impNetStatMat) = coefNames
  
  # parameters with entrainment adjustment to the density
  psi = matrix(data = NA, nrow = iterations, ncol = numberPara)
  colnames(psi) = coefNames
  
  # a list for the imputed networks
  impNetList = list()
  
  ## sequentially printing some iterations to make sure the sampler's progressing
  printIter = floor(seq(from = 1, to = iterations, length.out = 5))
  
  
  # get the list of free dyads
  # need the missingness indicator so we can work backwards from the adjmat
  missIndMat = (as.matrix(is.na(adjMat)) * 1)
  missIndEdgeList = as.edgelist(as.network(missIndMat, directed=FALSE), n = nrow(missIndMat))
  
  # miss model coef
  # first param (intercept/proportion of missingness) we can calculate from the proportion of missingness.
  missDens = sum(missIndMat)/(nodeSize * (nodeSize - 1))
  missCoef = c(log(missDens/(1-missDens)), entrainment)
  
  
  # specify the model
  substepModel = as.formula(paste("substepNet ~ ", as.character(formula)[3], sep = "" ))
  
  # another model for the proposed substep model
  propSubstepModel = as.formula(paste("propSubstepNet ~ ", as.character(formula)[3], sep = "" ))
  
  # grab the missingness model
  substepMissInd = as.network(missIndMat, directed = FALSE)
  
  # spceify a model, needs some flexibility eventually
  substepMissModel = as.formula(paste("substepMissInd ~ ", missModel, sep = ""))
  
  
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
      nextInitTheta = sampledThetas[iter-1, ]
      nextInitStats = impNetStatMat[iter-1, ]
    }
    
    ## Generating a proposal theta and auxiliary network
    # auxiliary parameters drawn from a multivariate normal distribution
    auxTheta = rmvnorm(n = 1, mean = nextInitTheta, sigma = propSigma)
    
    # generate a network using the auxiliary parameters
    auxNetStats = simulate( object = formula, 
                            coef = auxTheta,
                            output = "stats",
                            basis = naNet,
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
    
    
    # grab the specific network for the substep to initialise
    # for subsequent steps, use the most recently imputed network
    if(iter != 1){
      
      # use the previous imputed network
      substepNet = impNetList[[iter]]
      
    } else {
      # if it's the first step, use the adjmat with missingness
      substepNet = naNet
      
    }
    
    # prep the impnet object
    tempImpNet = substepNet
    
    ## Generation of imputed networks and statistics
    # starting the nested sampler
    # for all the missing tie variables, (missInd is 'r' in the notation)
    # since we're starting with the missing tie variables, we will ALWAYS attempt to impute (d_r = 1; d^*_r = 0)
    for(missInd in 1:nrow(missIndEdgeList)){
      
      # grab the dyad's nodes
      head = missIndEdgeList[missInd, 1]
      tail = missIndEdgeList[missInd, 2]
      
      
      # use a proposal substep network
      propSubstepNet = substepNet
      
      # compute the change statistic?
      startSubstepNetStats = summary(substepModel)
      
      # first evaluate the tie variable
      # if the value is missing, then we need to handle it differently.
      # this should only apply for the first iteration since we'll have imputed networks in the subsequent ones
      
      # the following if branch evaluates three possible cases with two possible outcomes
      # due to logical operations not working with missing values, we're going with a nested branch
      # we first consider the case where the selected tie variable is missing
      if(is.na(propSubstepNet[head, tail]) == TRUE){
        
        # realistically only taking place in the first iteration where x_r = NA
        propSubstepNet = add.edge(x = propSubstepNet,
                                  tail = tail,
                                  head = head)
      } else {
        # I despite nested branches, but I can't have the code rerun after adding the edge.
        if(propSubstepNet[head, tail] == 1){
          
          # grab the edge id
          deleteNetEdgeId = get.edgeIDs(x = propSubstepNet, v = head, alter = tail)
          
          # this is the only case where we would want to delete the edge to get a change statistic
          propSubstepNet = delete.edges(x = propSubstepNet,
                                        eid = deleteNetEdgeId)
        } else {
          # this branch evaluates the last case when the tie variable is a null tie (x^*_r = 0)
          
          # add an edge so we HAVE a change statistic
          propSubstepNet = add.edge(x = propSubstepNet,
                                    tail = tail,
                                    head = head)
        }
        
      }
      
      # re-compute the statistics
      endSubstepNetStats = summary(propSubstepModel)
      
      # grab the change statistics
      substepChangeNetStats = startSubstepNetStats - endSubstepNetStats
      
      # depending on whether an edge is added or deleted, compute the conditional probability
      if(sign(substepChangeNetStats[1]) == -1){
        
        # if we've deleted an edge then we can just compute the probability as usual
        condTieProb = 1/(1 + exp(sampledThetas[iter,] %*% substepChangeNetStats))
      } else {
        
        # if an edge is added, we need to subtract it 
        condTieProb = 1 - (1 / (1 + exp(sampledThetas[iter,] %*% substepChangeNetStats)))
      }
      
      # drawing x^*_r | x_{-r], thetas
      substepNetDraw = sample(c(0, 1), size = 1, prob = c((1 - condTieProb), condTieProb))
      
      # do the same for D, except a couple of differences
      # we always want to 'impute', so d_r = 1, d^*_r can be whatever
      # we use the imputed, or in this case, the substep network currently undergoing imputation as the 'true' network
      
      
      # computers can't sum across missing values so we need to temporarily zero-out the missingness to compute the starting
      # missingness model statisics
      # needs to be a matrix to be use as a dyadic covariate
      noNaSubstepNet = as.matrix(propSubstepNet)
      noNaSubstepNet[is.na(noNaSubstepNet)] = 0
      
      # overwriting the change stat'd tie variable with the drawn tie variable for the next step
      noNaSubstepNet[head, tail] = substepNetDraw 
      noNaSubstepNet[tail, head] = substepNetDraw 
      
      # grab the starting statistics
      startSubstepMissStats = summary(substepMissModel)
      
      # the following might not need to be a branch since we're assuming d_r = 1 at all times
      # no additional branch for missingness required since the missingness indicator can't have missing values
      if(substepMissInd[head, tail] == 1){
        
        # grab the edge id, this one is likely equivalent to the missInd index, but I don't know if it's consistent.
        deleteMissEdgeId = get.edgeIDs(x = substepMissInd, v = head, alter = tail)
        
        # delete the missingness indicator since that's the change
        substepMissInd = delete.edges(x = substepMissInd,
                                      eid = deleteMissEdgeId)
      } else {
        
        # this code SHOULD be useless since it's only if the chosen tie variable isn't missing (d_r = 0)
        # add an edge to have a change statistic
        substepMissInd = add.edge(x = substepMissInd,
                                  tail = tail,
                                  head = head)
      }
      # re-compute the statistics
      endSubstepMissStats = summary(substepMissModel)
      
      # grab the change statistics
      substepChangeMissStats = startSubstepMissStats - endSubstepMissStats
      
      # since we're assuming d_r = 1 all the time, this should ALWAYS be a positive value (reducing a missingness)
      if(sign(substepChangeMissStats[1]) == -1){
        
        # this step happens when we add a missingness indicator. It should never happen.
        condMissProb = 1/(1 + exp(missCoef %*% substepChangeMissStats))
      } else {
        
        # this step is computed when a missingness indicator is removed
        condMissProb = 1 - (1 / (1 + exp(missCoef %*% substepChangeMissStats)))
      }
      
      # drawing d^*_r | d_{-r], gamma/psi/whatever, x
      substepMissDraw = sample(c(0, 1), size = 1, prob = c((1 - condMissProb), condMissProb))
      
      # plug the draw in
      # both sides cause undirected.
      propSubstepNet[head, tail] = substepNetDraw
      propSubstepNet[tail, head] = substepNetDraw 
      
      # swap with some probability
      acceptRatio = (missCoef[2] * (substepMissDraw - substepMissInd[head, tail]) * (substepNet[head, tail] - substepTieDraw))
      + sampledThetas[iter, ] %*% (summary(substepModel) - summary(propSubstepModel))
      
      ## parameter swap
      if(log(runif(1)) < acceptRatio){
        
        # plug in the drawn tie variables (x_r <= x^*_r)
        tempImpNet[head, tail] = substepNetDraw
        tempImpNet[tail, head] = substepNetDraw
      } else {
        # nothing happens (x_r = x_r)
      }
    }
    
    # save the imputed networks
    impNetList[[iter]] = tempImpNet
    
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
  
  # and return it
  return(output)
}
