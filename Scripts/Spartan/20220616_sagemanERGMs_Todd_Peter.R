##### Script to run ERGMs on peter and todd networks
### Slurm array index
netInd <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# loading packages
library('ergm')
library('network')
library('here')

# loading network
load(here("Data", "20220616_sagemanERGMs_todd_peter_nets.RData"))       # this one is after my depletion
chosenTodd = toddNets[[netInd]]
chosenPeter = peterNets[[netInd]]

# loading dyadic covariates
load(here("Data", "sagemanClumpDyadcovs.RData")) 
load(here("Data", "sagemanSchoolDyadcovs.RData")) 

# label the dyad covs
clump.1.1.dyadcov = clumpMatrices[[1]]
clump.1.2.dyadcov = clumpMatrices[[2]]
clump.1.3.dyadcov = clumpMatrices[[3]]
clump.2.2.dyadcov = clumpMatrices[[5]]
clump.2.3.dyadcov = clumpMatrices[[6]]
clump.3.3.dyadcov = clumpMatrices[[8]]

school.1.1.dyadcov = schoolMatrices[[1]]
school.1.2.dyadcov = schoolMatrices[[2]]
school.2.2.dyadcov = schoolMatrices[[4]]

# model for Todd
# same configuration from before, now with the degraded network
model.2 <- chosenTodd ~ edges + isolates +gwesp(decay=log(2), fixed=TRUE)+gwdegree(decay=log(2), fixed=TRUE)+# density
  nodecov( 'place.joined') + # 
  nodematch( 'place.joined')+ # 
  nodematch( 'country.joined')+ # 
  nodecov('educ') + #
  absdiff('educ', pow = 2) +# 
  nodecov('age') + #
  absdiff('age', pow = 2) +#
  nodecov('year.joined') + #
  absdiff('year.joined', pow = 2) +# 
  nodecov('ses') + #
  absdiff('ses', pow = 2) +# 
  #  nodemix(  'clump',levels=c(1,2,3))+
  #  nodemix( 'school',levels=c(1,2))+
  dyadcov(clump.1.1.dyadcov) + 
  dyadcov(clump.1.2.dyadcov) + 
  dyadcov(clump.1.3.dyadcov) + 
  dyadcov(clump.2.2.dyadcov) + 
  dyadcov(clump.2.3.dyadcov) + 
  dyadcov(clump.3.3.dyadcov) + 
  dyadcov(school.1.1.dyadcov) + 
  dyadcov(school.1.2.dyadcov) + 
  dyadcov(school.2.2.dyadcov) + 
  nodematch( 'plbirth')

# ergm
mod2ergm = ergm(model.2, control = control.ergm(init.method = "CD", MCMC.burnin = 50000, MCMC.interval = 50000))
mod2summary = summary(mod2ergm)

# model for Peter
# same configuration from before, now with the degraded network
model.3 <- chosenPeter ~ edges + isolates +gwesp(decay=log(2), fixed=TRUE)+gwdegree(decay=log(2), fixed=TRUE)+# density
  nodecov( 'place.joined') + # 
  nodematch( 'place.joined')+ # 
  nodematch( 'country.joined')+ # 
  nodecov('educ') + #
  absdiff('educ', pow = 2) +# 
  nodecov('age') + #
  absdiff('age', pow = 2) +#
  nodecov('year.joined') + #
  absdiff('year.joined', pow = 2) +# 
  nodecov('ses') + #
  absdiff('ses', pow = 2) +# 
  #  nodemix(  'clump',levels=c(1,2,3))+
  #  nodemix( 'school',levels=c(1,2))+
  dyadcov(clump.1.1.dyadcov) + 
  dyadcov(clump.1.2.dyadcov) + 
  dyadcov(clump.1.3.dyadcov) + 
  dyadcov(clump.2.2.dyadcov) + 
  dyadcov(clump.2.3.dyadcov) + 
  dyadcov(clump.3.3.dyadcov) + 
  dyadcov(school.1.1.dyadcov) + 
  dyadcov(school.1.2.dyadcov) + 
  dyadcov(school.2.2.dyadcov) + 
  nodematch( 'plbirth')

## If necessary, initialise with the full model's coefficients, but first try without it.
### I DO run into the initialisation error with Peter
## I'll just read in an old copy of the full ergm
load(here("Output", "20220531_sagemanFullSummary_CD.RData"))

# ergm
mod3ergm = ergm(model.3, control = control.ergm(init = mod1ergm$coef, init.method = "CD", MCMC.burnin = 50000, MCMC.interval = 50000))
mod3summary = summary(mod3ergm)

save(mod2ergm, mod2summary, mod3ergm, mod3summary,  file = here("Output", "20220616_sagemanERGMs_Todd_Peter", paste("20220616_sagemanERGMs_Todd_Peter_", netInd, ".RData",sep = "") ))
