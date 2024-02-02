# Spartan script to subset the chosen networks
# read in the array index for the proportion of missingness chosen
propIndex = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# set the proportions of missingness available
propChoices = c(1, 3.5, 6)

# and choose
chosenProp = propChoices[propIndex]

# here
library(here)

# the index used to select networks for its reestimation
trialIndices = floor(seq(from = 501, to = 1000, length.out = 50))

# prepare lists to save the missingness matrices in
chosenIndep = list()
chosenMarErgm = list()
chosenLatent = list()
chosenMnarErgm = list()

# read those chunky files for all the various missingness models
# these are indep/mar/latent
load(here("Output", "20230814_simMissNets", paste("20230825_simMissNets_propMiss", chosenProp,".RData", sep = "")))

# loop for all 6 networks
for(networkIndex in 1:6){

  chosenIndep[[networkIndex]] = indepMissMats[[networkIndex]][trialIndices,,]
  chosenMarErgm[[networkIndex]] = missErgmSimNets[[networkIndex]][trialIndices]
  chosenLatent[[networkIndex]] = missLatentSimNets[[networkIndex]][trialIndices]

}

# and this one is MNAR, done separately because the object has the same name as the marErgm
load(here("Output", "20230814_simMissNets", paste("20231011_simMissErgmNets_ergmCoef_1_propMiss", chosenProp,".RData", sep = "")))

# loop for all 6 networks
for(networkIndex in 1:6){

  chosenMnarErgm[[networkIndex]] = missErgmSimNets[[networkIndex]][trialIndices]
}

# save all as matrices
indepMissNets = lapply(chosenIndep, function(x){apply(x, MARGIN = c(1), FUN = network::as.network, directed = FALSE)})
indepMissAdjs = lapply(indepMissNets, function(x){lapply(x, FUN = as.matrix)})
marErgmMissAdjs = lapply(chosenMarErgm, function(x){lapply(x, FUN = as.matrix)})
latentMissAdjs = lapply(chosenLatent, function(x){lapply(x, FUN = as.matrix)})
mnarErgmMissAdjs = lapply(chosenMnarErgm, function(x){lapply(x, FUN = as.matrix)})


# save the output as the chosen subset of simMissNets
compiledMissNets = list(Indep = indepMissAdjs,
                        marErgm = marErgmMissAdjs,
                        latent = latentMissAdjs,
                        mnarErgm = mnarErgmMissAdjs)

save(compiledMissNets, file = here("Output", "20230814_simMissNets", paste("20231024_chosenMissNets_prop",chosenProp,".RData", sep = "")))