## Plotting functions for degraded network simulations and estimations
# predominantly copied over from 20230821_missModResults

### Function to turn the data from the extracted lists into plot-friendly format

# start from the list of coefficients and standard errors
# these can have varying rows per model because not all re-estimations converge


prepMissPlot = function(modelResList, modelSeList, missModelList, propMissList, chosenProp, chosenPara){
  
  ## prepMissPlot(modelResList, modelSeList, missModelList, propMissList, chosenProp, chosenPara) takes 
  ## four different containing various estimated values and extracted indices to be used 
  ## in making a plot-friendly dataset.
  ##
  ## Input:
  ## - modelResList: A k-long list containing the estimated model parameters for k re-estimated models.
  ##
  ## - modelSeList: A k-long list containing the estimated standard errors for k re-estimated models.
  ##                
  ## - missModelList: A k-long list containing the chosen missingness model index for k re-estimated models.
  ##                  
  ## - propMissList: A k-long list containing the missingness proportion index for k re-estimated models.
  ##
  ## - chosenProp: The specified proportion of missingness to plot. Three possible values here:
  ##               1 = 10%, 3 = 35%, 6 = 60%.
  ##
  ## - chosenPara: The specified parameter to plot. Three possible values here:
  ##               "edges" for edges, "altStar" for alternating stars, "gwesp" for gwesp.
  ##
  ## Output: 
  ## - A k x 3 dataframe containing the parameter values, standard error values, and missing model index.
  ##   The dataframe is ordered by parameter value (smallest to largest) for each missingness model.
  
  
  ## Check to make sure all the lists are of the same length
  listLengths = lapply(list(modelResList, modelSeList, missModelList, propMissList), FUN = length)
  
  ## A unit test that spits out an error if there's more than one unique length in the list of lengths
  if( length(unique(listLengths)) != 1 ){
    stop("Lists have differing lengths")
  }
  
  # subset the lists to only be a single chosen proportion
  chosenPropIndex = propMissList == chosenProp
  
  modelResChosenPropList = modelResList[chosenPropIndex]
  modelSeChosenPropList = modelSeList[chosenPropIndex]
  missModelChosenPropList = missModelList[chosenPropIndex]
  
  # Specify indices for the chosen parameter
  # doing it this way because parameter names get very lengthy
  paraIndex = c(1:3)
  modelParaNames = c("edges", "altStar", "gwesp")
  chosenParaIndex = paraIndex[modelParaNames == chosenPara]
  
  
  # making the dataset
  plotData = data.frame(parameter = c(unlist(lapply(modelResChosenPropList, `[[`, chosenParaIndex) )), 
                        SE = c(unlist(lapply(modelSeChosenPropList, `[[`, chosenParaIndex))),
                        missModel = factor(c(unlist(missModelChosenPropList)), levels = c(1,2,3), labels = c("Indep", "ERGM", "Latent")))
  
  # order the plot data
  plotOrd = plotData[order(plotData[,"parameter"]),]
  
  # one more ordering with the missing model type
  plotOrd = plotOrd[order(plotOrd[,"missModel"]),]
  
  # return the ordered plot data
  return(plotOrd)
}

## actual plotting function

missCaterpillarPlot = function(plotData, truePara, trueSe, chosenPara, chosenNetInd, chosenProp, chosenMiss){
  
  ## missCaterpillarPlot(plotData, truePara, trueSe, chosenPara, chosenNetInd, chosenProp) is a very specialised
  ## plotting function that produces a caterpillar plot to compare parameter estimates for very specific data.
  ##
  ## Input:
  ## - plotData: Needs to be a k x 3 data frame for k re-estimated models containing parameter values,
  ##             standard error values, and a missing model index. Ordering is optional, but visually
  ##             useful. Ideally from prepMissPlot().
  ##
  ## - truePara: A single float (can technically be integer, but unlikely). Represents the true model estimate.
  ##                
  ## - trueSe: A single float (can technically be integer, but unlikely). Represents the true standard error.
  ##
  ## - chosenPara: A string with three (but really two) options. The chosen string is the parameter that is 
  ##               plotted. The options, in order, are "edges", "altStar", and "gwesp"
  ##
  ## - chosenNetInd: An integer between 1 to 6 for current purposes. Each integer represents one of the six
  ##                 datasets from UCINet chosen for the current round of re-estimations.
  ##                
  ## - chosenProp: A single integer or float with three options representing the proportion of missingness.
  ##               Options are 1 = 10%, 3.5 = 35%, 6 = 60%.
  ##
  ## - chosenMiss: A value to indicate what missingness is saved as (e.g., NA is Peter, 0 is Todd).
  ##               Input should be strings that are directly used for the label.
  ##                  
  ## Output: 
  ## - The output should be a caterpillar plot containing with k datapoints for k re-estimated models.
  ##   Colours are labelled and represent the missingness model used. True values are represented as a 
  ##   grey rectangle. All error bars or other depictions of spread are 95% confidence intervals. Plot
  ##   title should also be adaptive to the parameters in the function.
  
  
  # caterpillar
  caterPlot = ggplot( data = plotData,
                      aes( x = 1:nrow(plotData), y = parameter, col = missModel)) + 
    geom_errorbar(aes (ymin = (parameter - 1.96*SE), ymax = (parameter + 1.96*SE), width = 0.4)) + 
    xlab("") + 
    ylab(paste(chosenPara,"estimate (95% Confidence Int)", sep = " ")) + 
    labs(col = "missType") + 
    geom_hline(yintercept = truePara, col = "darkblue") + 
    geom_hline(yintercept = 0, col = "black", lty = 2) +
    geom_point() + 
    ggtitle(paste("Diagnostics - Net", chosenNetInd, " - Prop", chosenProp, " - Miss", chosenMiss ,sep = "")) + 
    theme_classic() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    geom_rect(aes(xmin = -Inf,
                  xmax = Inf, 
                  ymin = (truePara - 1.96*trueSe),
                  ymax = (truePara + 1.96*trueSe)),
              alpha = 0.01,
              colour = NA,
              fill = "grey") +
    ylim(min(min(plotData$parameter) - 2*max(abs(plotData$SE)), (truePara - 1.96*trueSe)), 
         max(max(plotData$parameter) + 2*max(abs(plotData$SE)), (truePara + 1.96*trueSe)))
  
  # print out the plot
  caterPlot
  
}

# function to calculate metrics
netMetrics = function(net, directed = FALSE){
  
  ## netMetrics(net, directed) calculates some metrics for a given network. Needs the 'sna' package
  ## and the 'network' package to work.
  ##
  ## Input:
  ## - net: Ideally a network object for the network whose metrics will be calculated. Technically
  ##        can be an adjacency matrix since it will be converted to a network object.
  ##
  ## - directed: A logical value to indicate whether the network is directed or not.
  ##                  
  ## Output: 
  ## - The output is a list of different network metrics. These include the clustering coefficient,
  ##   the global centralisation (degree variance) value, the median geodesic distance, and median
  ##   betweenness centrality. More metrics can be added on as necessary.
  
  # needs two packages
  require(network)
  require(sna)
  
  # specify if undirected or not
  mode = 'graph'
  
  # if instead directed network are the input, turn it directed
  if(directed == TRUE){
    mode = 'digraph'
  } 
  
  # if input is an adjacency matrix, turn it into a network object
  if(inherits(net, "matrix")){
    net = as.network(net, directed = directed) 
    print("Matrix coerced to network object")
    
  }
  
  # make sure there's no missingness in it
  if(network.naedgecount(net) != 0){
    
    tempMat = as.matrix(net)
    tempMat[is.na(tempMat)] = 0
    net = as.network(tempMat, directed = directed)
    print("Missing values in the network converted to 0s")
    
  }
  
  # getting a distribution of clustering coefficients for the 80 simulated graphs
  tempClust <- gtrans(net, mode = mode)
  tempCent <- centralization(net, degree, mode = mode)
  tempGeodistMedian <- median(geodist(net, inf.replace = NA)$gdist, na.rm = T)
  tempBetwMedian <- median(betweenness(net, gmode = mode))
  
  ## specify the chosen metric
  metricList = list("clustering" = tempClust,
                    "centralisation" = tempCent,
                    "MedianGeodesic" = tempGeodistMedian,
                    "MedianBetween" = tempBetwMedian)
  
  # return the list of metrics
  return(metricList)
  
}


## get mean value parameterisations
# function to check mean value parameterisations
meanValPar = function(net, model){
  
  ## net is the outcome network, usually the degraded adjacency matrix
  ## model is a STRING for the specified model
  
  # paste together the input network object and the specified model
  eval(parse(text = paste("chosenModel = net", "~", model )))
  
  # check mean value par
  summary(chosenModel)
}

## don't think this function below is used very much
## but *can* be used simulate various missingness mechanisms with a specfiied model
geneSimNets = function(net, modelSpec, initParas, chosenPara, paraRange, trueNet){
  
  # takes a network, like a depleted network without missingness (i.e., saved as 0)
  # use different models to simulate networks with the same amount of edges
  # i.e., simulating the effects of various parameter values (and model specification)
  # on a depleted network with constrained edges/density
  # ideally accommodates a range of different model specifications and parameter values
  # the parameters are varied one at a time to view its effect holding constant on other parameters....
  # so the function will look like
  # simNet(toddNet, modelSpec, initParas, chosenPara, paraRange, trueNet)
  # output = list of simulated networks, 
  # mean value parameterisations can be applied in a separate chunk of code for the list
  # of simulated networks since it'd be conditional to the model specified for the simulated networks.
  
  ### NO:
  # take a simulated missingness matrix with a fixed proportion
  # specify some model, with varying specification or parameter values,
  # then simulate some other possible missingness matrices
  # while constraining the edges (fixed prportion of missingness)
  # use the simulated missingness matrices to deplete true adjmat X
  # then calculate metrics and mean value parameterisations on the depleted network
  
  ## checks
  # make sure it's a network object
  if(inherits(net, "matrix")){
    net = as.network(net)
  }
  
  # make sure there's no missingness in it
  if(network.naedgecount(testNet) != 0){
    stop("Network has missing values in it.")
  }
  
  # degree covariate calculation, should I retroactively update the reestimation code to attach the degree covariate...
  # might be better to specify the true network here...
  net %v% 'degree' = rowSums(trueNet)
  
  # initialising the simulation
  netStart = net
  simNetList = list()
  
  # loop to simulate different gwdeg values
  for (paraIndex in c(1:length(paraRange)) ){
    
    # a conditional after the first loop for making the starting network iterate between gwdeg parameters
    if ( paraIndex > 1 ){
      netStart <- tempSimNets[[1]]
      netStart %v% 'degree' = rowSums(trueNet)
    }
    
    # choose coefs
    chosenCoefs = initParas
    chosenCoefs[chosenPara] = paraRange[paraIndex]
    
    # simulating 80 networks given the specified model and their parameter values
    tempSimNets <- simulate(object = modelSpec,
                            constraints =~edges,
                            coef= chosenCoefs,
                            control=control.simulate(
                              MCMC.burnin=20000,
                              MCMC.interval=100), nsim=80 )
    
    # save the simulated networks
    simNetList[[paraIndex]] = tempSimNets
    
  }
  
  # return the simulated 'networks' (these are simulated missingness matrices saved as networks)
  return(simNetList)
  
}

## Function to prep metric plot data
prepMetricPlot = function(degradedNetMetricList, chosenPara){
  
  ## prepMetricPlot(degradedNetMetricList, chosenPara) takes a list of metric calculated from degraded networks
  ## and a specified parameter name that controlled the degradation to subset them into a data frame to be used
  ## for plotting purposes.
  ##
  ## Input:
  ## - degradedNetMetricList: A list of (nested) network metrics calculated using the netMetrics function on
  ##                          a list of degraded networks. Missingness needs to be saved as 0 for values to calculable.
  ##
  ## - chosenPara: A string to describe which parameter is being varied in the degradation of the networks.
  ##                  
  ## Output: 
  ## - A dataframe containing the metric names, which parameter is being varied and what its value is for the degradation,
  ##   and three quantiles (5, 50, and 95) to describe the range and median of the degraded metric.
  
  
  # extracting specific metrics
  uniqueMetrics = unique(names(unlist(degradedNetMetricList[[1]])))
  
  # initialise plot data
  metricPlotData = data.frame(metricName = NA, varyingParameter = NA, parameterValue = NA, quant0.05 = NA, median = NA, quant0.95 = NA)
  
  # additional list for metric ranges
  rawMetricList = list()
  
  # loop to choose different metrics in calculating the different quantiles
  for(chosenMetric in 1:length(uniqueMetrics)){
    
    # a value for the row index based on the iteration of the loop
    tempRowInd = seq(from = 1 + (length(degradedNetMetricList) * (chosenMetric - 1)), to = chosenMetric * length(degradedNetMetricList))
    
    # fix this value for now, but this would depend on the model spec and which chosen para it is
    chosenParaLabel = chosenPara
    metricPlotData[tempRowInd,"varyingParameter"] = chosenParaLabel
    
    # put in the metric name
    metricPlotData[tempRowInd, "metricName"] = uniqueMetrics[chosenMetric]
    
    # extract the chosen metrics from the nested list
    rawMetricList[[chosenMetric]] = lapply(degradedNetMetricList, function(chosenParaVal){unlist(chosenParaVal)[names(unlist(chosenParaVal)) == uniqueMetrics[[chosenMetric]]]})
    
    # fill in the quantiles
    metricPlotData[tempRowInd, c("quant0.05", "median", "quant0.95")] = do.call(rbind, lapply(rawMetricList[[chosenMetric]], function(chosenParaVal){quantile(chosenParaVal, c(0.05, 0.5, 0.95))}))
    
    # lastly, give the parameter values
    metricPlotData[tempRowInd, "parameterValue"] = paraRange
  }
  
  # return the plot data
  return(metricPlotData)  
}

## Turn this base metric plot code into a function-
missMetricPlot = function(data, chosenMetric, chosenPara, chosenProp){
  
  ## missMetricPlot(data, chosenMetric, chosenPara, chosenProp) is a very specialised plotting function that
  ## produces a caterpillar plot to compare parameter estimates for very specific data.
  ##
  ## Input:
  ## - data: A data frame containing the metric name, varying parameter, the varying parameter's value, and 
  ##         5th, 50th, and 95th quantiles of the metric. Obtained from prepMetricPlot().
  ##
  ## - chosenPara: A string to describe which parameter is being varied in the degradation of the networks.
  ##                
  ## - chosenProp: A single integer or float with three options representing the proportion of missingness.
  ##               Options are 1 = 10%, 3.5 = 35%, 6 = 60%.
  ##
  ## Output: 
  ## - The output should be a line plot with a shaded area representing the 90% confidence interval of the metric.
  ##   Its labels should be adaptive to the chosen metric, (varying) parameter, and proportion of missingness.
  
  
  # requires ggplot2 and dplyr
  require(ggplot2)
  require(dplyr)
  
  # and plot
  data %>%
    filter(metricName == chosenMetric) %>% 
    ggplot(aes(x = parameterValue, y = median)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = quant0.05, ymax = quant0.95), alpha = 0.2) + 
    labs(title = paste("Depleted network diagnostics - varying", chosenPara," - Prop", chosenProp, sep = "")) + 
    xlab(paste("Varying", chosenPara, "value")) + 
    ylab(paste(toupper(substr(chosenMetric, 1, 1)), substr(chosenMetric, 2, nchar(chosenMetric)), sep = ""))
  
}