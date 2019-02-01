setwd("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/gam")
source("./makeMod.R");source("./getBreaks.R")
fLevs <- read.csv('C:/Users/Maia Kapur/Dropbox/UW/sab-growth/iPopSim/inputs/scenarios.csv',na.strings = 'NA') ## manual file

require(mgcv);require(dplyr)

for(l in c(1,35)){
  
  ## get scenario name
  scen0 <- paste0(fLevs[l,'DESC'])
  scen1 <-  ifelse(is.na(fLevs[l, 5]),
                   paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
                   paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))
  scen <- paste(scen0,scen1,sep = "_")
  dat <-read.csv(paste0("C:/users/maia kapur/dropbox/uw/sab-growth/ipopsim/gendata/",
                           scen,'.csv')) %>% filter(Age == 4)
  mod <- makeMod(scenario = scen,dat)
  getBreaks(mod,dat, scen)
  
}
