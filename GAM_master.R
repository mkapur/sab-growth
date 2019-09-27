## Code to test GAM method on datasets generated via iPopSim (IBM)
## This will fit a gam for each replicate, identify break points, fit new life history params
## and plot comparisons
## Maia Sosa Kapur kapurm@uw.edu Forever... 2019

require(mgcv);require(dplyr);require(ggplot2); require(TMB); library(reshape2)
library(gridExtra); library(grid); library(lattice)

compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))
compname <- c("Maia Kapur","mkapur")[1]

## Build GAMS, get breakpoints ----
scenarios <- read.csv("./input_data/scenarios.csv",na.strings = 'NA')## manual file
age <- 6; nboot <- 100; #testrows <- unique(scenarios$DESC)
# source("./functions/bootBreaks.R") ## about 10 mins -- don't need to do this >1x, but yes if re-gen IBM
source("./functions/getGR.R");source("./functions/fitMod.R");source("./functions/missby.R")


# Fit VB and calc coverage ----
ldfprop <-  read.csv( paste0("./GAM_output/ldf_raw_a6.csv")) ## dataframe of detected breaks made above○
rm(cdf); cdf <- data.frame(); idx <-1 ## storage coverage prob totals, rbind each scen
for(l in 1:length(unique(ldfprop$scen))){
  # for(b in sample(1:nboot,5)){
  for(b in 1:nboot){
    ## TMB FITTING ----
    scen <- unique(ldfprop$scen)[l]
    tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_",b,".csv")) #%>% filter(Age == 6)
    if(scen == 'NoBreaks') tempdf$REG <- as.factor('R1')
    breaksdf <- read.csv( paste0("./GAM_output/ldf_raw_a",age,".csv")) %>% filter(scen ==  unique(ldfprop$scen)[l] & boot == b)
    dat<-getGR(tempdf,breaksdf);rm(tempdf)
    
    ## generate DES matrix of vectors and a KEY for later comparison
    DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat))
    dat$cREG <-  paste0(dat$gamREG,"_",dat$Period)
    # DES <- ifelse(!is.na(dat$gamREG), as.numeric(factor(dat$gamREG)),1)-1 ## overwrite NA for nobreaks
    DES <- ifelse(!is.na(dat$cREG), as.numeric(factor(dat$cREG)),1)-1 ## overwrite NA for nobreaks
    KEY <- paste(scen,DES,sep = "_")
    # keybase <- unique(ifelse(!is.na(dat$gamREG), as.character(dat$gamREG),"R1")) ## text regions
    
    temp <- data.frame(DES = as.numeric(as.character(DES)), cREG = dat$cREG); temp <- unique(temp)
    keybase <- paste(as.factor(temp[order(temp$DES),'cREG']));rm(temp)
    
    ## run TMB parest----
    nStrata <- length(unique(DES))
    data <-
      list(
        Length_cm = dat[,"Length_cm"],
        Age = dat[,"Age"],
        DES = as.vector(DES),
        nStrata = nStrata,
        a1 = 3,
        a2 = 30
      )

    parameters <-
      list(
        log_Linf = rep(log(75), nStrata),
        log_k = rep(log(0.3), nStrata),
        t0 = rep(3, nStrata),
        log_Sigma = 0.1
      )

    # Now estimate everything
    map <- NULL
    mod <- fitMod(data, parameters, modversion = 'sptlvb',map)
    dat0 <- mod[[1]]
    rep0 <- mod[[2]]
    
   ## reformat outputs ----
    write.csv(rep0, file = paste0("./output_data/",scen,"_parEst_gam_",b,"_",Sys.Date(),'.csv'),row.names = F)
    ypreds0 <- data.frame(dat0); names(ypreds0)[1] <- c('Predicted')
    write.csv(ypreds0,  paste0("./output_data/",scen,"_predicts_",b,Sys.Date(),".csv"),row.names = F)
    cat("Fit TMB model ",paste(scen)," GAM boot ",b," & saved outputs \n")

    ## check and save coverage ----
    parms <- c('L1','L2','k','Linf')[1:2]
    vars <- paste0('R',1:4,"_");   vis <- c("pooled",'early','late'); Rlevs <- c("R1", "R2","R3","R4",apply(expand.grid(vars, vis), 1, paste, collapse=""))
    source("./functions/compareBreaks2.R")
    cat("Tabulated Breaks ",paste(scen)," GAM boot ",b," CDF row ",idx," \n")
    
  } ## end boots
} ## end ldfrows

cdf %>%
  mutate(method = "GAM") %>%
  write.csv(.,file = paste0('./gam_output/cdf_',Sys.Date(),'.csv'),row.names = F)

## Get coverage probs -- for all params. 
## Note N varies because of varying # regions per boot. 
cdfprop <- cdf %>%
  filter(!is.na(scen) & !is.na(L2))  %>% 
  select(scen, L1, L2) %>% 
  # mutate(both = (L1 == T & L2 == T)) %>%
  melt(id = c('scen')) %>%
  group_by(scen,variable) %>%
  dplyr::summarise(denom = n(), n = sum(value)) %>% 
  mutate(prop = round(n/denom,2))%>%
  mutate(method = "GAM")
write.csv(cdfprop,file = paste0('./gam_output/cdf_prop_',Sys.Date(),'.csv'),row.names = F)

## When did regional designation go right? (original analysis)
cdfaccu <- cdf %>% 
  select(scen, LAT, LON, YEAR) %>%
  # mutate(both = (LAT == T & LON == T), all = (both == T & YEAR == T)) %>%
  melt(id = c('scen')) %>%
  group_by(scen,variable) %>%
  dplyr::summarise(denom = n(), n = sum(value)) %>% 
  mutate(prop = round(n/denom,2)) %>%
  mutate(method = "GAM")
write.csv(cdfaccu,file = paste0('./gam_output/cdf_accu_',Sys.Date(),'.csv'),row.names = F)
