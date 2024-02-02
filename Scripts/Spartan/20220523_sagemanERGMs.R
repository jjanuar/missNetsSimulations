#### Sageman ERGMs

# loading packages
library('ergm')
library('network')
library('here')

# loading network
load(here("Data", "sagemanNetworks.RData"))       # this one is after my depletion
load(here("Data", "SagemanImputedData.RData"))    # this one comes from Johan's script.

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


# using the non-depleted network
# models
model.1 <- alqnet ~ edges + isolates +gwesp(decay=log(2), fixed=TRUE)+gwdegree(decay=log(2), fixed=TRUE)+# density
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
mod1ergm = ergm(model.1, control = control.ergm(init.method = "CD", MCMC.burnin = 50000, MCMC.interval = 50000))
mod1summary = summary(mod1ergm)

save(mod1ergm, mod1summary, file = here("Output", "20220531_sagemanFullSummary_CD.RData"))

# models
# same configuration from before, now with the degraded network
model.2 <- sagemanNetDegraded0 ~ edges + isolates +gwesp(decay=log(2), fixed=TRUE)+gwdegree(decay=log(2), fixed=TRUE)+# density
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

save(mod2ergm, mod2summary, file = here("Output", "20220531_sagemanToddSummary_CD.RData"))

# same configuration from before, now with the degraded network
model.3 <- sagemanNetDegraded ~ edges + isolates +gwesp(decay=log(2), fixed=TRUE)+gwdegree(decay=log(2), fixed=TRUE)+# density
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
mod3ergm = ergm(model.3, control = control.ergm(init = mod1ergm$coef,init.method = "CD", MCMC.burnin = 50000, MCMC.interval = 50000))
mod3summary = summary(mod3ergm)

save(mod3ergm, mod3summary, file = here("Output", "20220531_sagemanPeterSummary_CD.RData"))
