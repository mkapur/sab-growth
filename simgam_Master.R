## Code to test GAM method on datasets generated via iPopSim (IBM)
## This will fit a gam for each replicate, identify break points, fit new life history params
## and plot comparisons
## Maia Sosa Kapur kapurm@uw.edu WiSp 2019

require(mgcv);require(dplyr);require(ggplot2); require(TMB); library(reshape2)
library(gridExtra); library(grid); library(lattice)

compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))
# compile("sptlschnute.cpp") 
# dyn.load(dynlib("sptlschnute"))

compname <- c("Maia Kapur","mkapur")[1]

## Build mods, get breakpoints
scenarios <- read.csv("./input_data/scenarios.csv",na.strings = 'NA') ## manual file
age <- 6; nboot <- 100; testrows <- unique(scenarios$DESC) 
# source("./functions/bootBreaks.R") ## about 10 mins -- don't need to do this >1x
source("./functions/getGR.R");source("./functions/fitMod.R")
# now loop boots and fit VB to get coverage probs
ldfprop <-  read.csv( paste0("./GAM_output/ldf_raw_a6.csv"))
cdf <- data.frame();idx <-1 ## storage coverage prob totals, rbind each scen

compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}


for(l in 1:length(unique(ldfprop$scen))){ 
  for(b in 1:nboot){
    
    ## get specific scenario and re-aggregate
    scen <- unique(ldfprop$scen)[l]
    tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_",b,".csv")) #%>% filter(Age == 6)
    if(scen == 'NoBreaks') tempdf$REG <- as.factor('R1')
    breaksdf <- read.csv( paste0("./GAM_output/ldf_raw_a",age,".csv")) %>% filter(scen ==  unique(ldfprop$scen)[l] & boot == b)
    dat<-getGR(tempdf,breaksdf);rm(tempdf)
    
    ## now re-aggregate the data ----
    ## generate DES matrix of vectors and a KEY for later comparison
    DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat)) 
    DES <- ifelse(!is.na(dat$gamREG), as.numeric(factor(dat$gamREG)),1)-1 ## overwrite NA for nobreaks
    KEY <- paste(scen,DES,sep = "_")
    keybase <- unique(ifelse(!is.na(dat$gamREG), as.character(dat$gamREG),"R1")) ## text regions
    
    ## run TMB parest----
    nStrata <- length(unique(DES))
    
    data <-
      list(
        Length_cm = dat[,"Length_cm"],
        Age = dat[,"Age"],
        DES = as.vector(DES),
        nStrata = nStrata
      )
    
    # parameters <-
    #   list(
    #     log_Lone = rep(log(60), nStrata),
    #     log_Ltwo = rep(log(275), nStrata),
    #     log_k = rep(log(0.3), nStrata),
    #     t0 = rep(0.1, nStrata),
    #     log_Sigma = 0
    #   )
    
    parameters <-
      list(
        log_Linf = rep(log(150), nStrata),
        log_k = rep(log(0.3), nStrata),
        t0 = rep(0.1, nStrata),
        log_Sigma = 0
      )
    
    # Now estimate everything -- DO PHASES
    map <- NULL
    mod <- fitMod(data, parameters, modversion = 'sptlvb',map)
    dat0 <- mod[[1]]
    rep0 <- mod[[2]]
    # map <- list(log_Lone = rep(factor(NA),nStrata)) ## get K, other stuff shouldn't change much
    # model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)
    # phase1 <- fitMod(data, parameters, modversion = 'sptlvb',map)[[2]]
    # parameters <-
    #   list(
    #     log_Lone = rep(log(60), nStrata),
    #     log_Ltwo = log(phase1[phase1$variable == 'Ltwo','value']),
    #     log_k = rep(log(0.3), nStrata),
    #     t0 = rep(0.1, nStrata),
    #     log_Sigma = 0
    #   )
    # map <- list(log_Ltwo = rep(factor(NA),nStrata)) ## get K, other stuff shouldn't change much
    # phase2 <- fitMod(data, parameters, modversion = 'sptlschnute',map)
    # 
    # dat0 <- phase2[[1]]
    # rep0 <-phase2[[2]]
    # rep0[rep0$variable == 'Ltwo','sd'] <-     phase1[phase1$variable == 'Ltwo','sd']

    
   ## reformat outputs ----
    write.csv(rep0, file = paste0("./output_data/",scen,"_parEst_gam_",b,"_",Sys.Date(),'.csv'),row.names = F)

    # ypreds0 <- cbind(dat0,dat) %>% data.frame()
    ypreds0 <- data.frame(dat0)
    names(ypreds0)[1] <- c('Predicted')

    write.csv(ypreds0,  paste0("./output_data/",scen,"_predicts_",b,Sys.Date(),".csv"),row.names = F)
    ypreds <-  read.csv(paste0("./output_data/",scen,"_predicts_",b,Sys.Date(),".csv"))
    
    cat("Fit TMB model ",paste(scen)," boot ",b," & saved outputs \n")

    ## check and save coverage --
    ## plot estimates
    gamREGS <- factor(unique(dat$gamREG), levels=c("R1", "R2","R3","R4"))
    REGS <-  factor(unique(dat$REG), levels=c("R1", "R2","R3","R4")) ## what was used to GENERATE
    parms <- c('L1','L2','k','Linf')[1:2]
    
    parest <-  read.csv(paste0("./output_data/",scen,"_parEst_gam_",b,"_",Sys.Date(),'.csv')) %>%
      filter(variable %in% parms) %>% mutate(source = 'Estimated')
    
    parest <- rbind(parest, (read.csv("./input_data/true_ibm_vals.csv") %>% 
                               filter(REG %in% REGS & variable %in% parms)))    ## only extract regions from 'actual' that were present in original

    parest$REG <- factor(parest$REG, levels=c("R1", "R2","R3","R4"))
    parest$lwr <- parest$value - 1.96*parest$sd
    parest$upr <- parest$value + 1.96*parest$sd

   
    for(iq in 1:length(unique(parest$REG))){ ## loop regions
      cdf[idx,'scen'] <- scen
      cdf[idx,'boot'] <- b
      ## qualitative compare for regional
      cdf[idx,'gamREG'] <- if(!is.na(gamREGS[iq])){gamREGS[iq]}else{last(gamREGS)}
      cdf[idx,'REG'] <- if(!is.na(REGS[iq])){REGS[iq]}else{last(REGS)}
      
      ## direct comp for break locations
      cdf[idx,'gamLAT'] <- breaksdf$lat_breaks
      cdf[idx,'gamLON'] <- breaksdf$lon_breaks
      
      ptf <- NULL
      for(v in  1:length(parms)){ ## loop variables
        tmp <- subset(parest, variable == parms[v] & source == 'Estimated'  | variable == parms[v] & source == 'Actual' & REG == cdf[idx,'REG'])
        bounds <- tmp %>% filter(source == 'Estimated' & REG ==  cdf[idx,'gamREG'] ) %>% select(lwr,upr)
        ptf[v] <- tmp$value[tmp$source == 'Actual' &  tmp$REG ==    cdf[idx,'REG']] >= bounds[1] &
          tmp$value[tmp$source == 'Actual'&  tmp$REG ==    cdf[idx,'REG']] <= bounds[2]
      } ## end vars
      cdf[idx,parms[1]] <- ptf[1];  cdf[idx,parms[2]] <- ptf[2]
      # cdf[idx,"LAT"] <- cdf[idx,"gamLAT"] == scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1]
      # cdf[idx,"LON"] <- cdf[idx,"gamLON"] == scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1]
   
      # if(is.na(cdf[idx,"gamLAT"])  | is.na(cdf[idx,"gamLON"])){ ## won't match NA
        # cdf[idx,"LAT"] <- ifelse(is.na(cdf[idx,"gamLAT"]) & is.na(scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1]),TRUE,FALSE )
        # cdf[idx,"LON"] <- ifelse(is.na(cdf[idx,"gamLON"]) & is.na(scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1]),TRUE,FALSE )
        cdf[idx,"LAT"] <- compareNA(cdf[idx,"gamLAT"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
        cdf[idx,"LON"] <- compareNA(cdf[idx,"gamLON"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
      # }
        if(cdf[idx,'scen']=="F0LMW")  { # won't match "overlap"
          cdf[idx,"LAT"] <- cdf[idx,"gamLAT"] >= 20 & cdf[idx,"gamLAT"] <= 25
          cdf[idx,"LON"] <- cdf[idx,"gamLON"] >= 20 & cdf[idx,"gamLON"] <= 25
        }
      if(is.na(cdf[idx,"LAT"])) stop(paste0("Throwing NAs in spatial comparison ", scen, " ", idx, "\n"))

      idx <- idx+1
    } ## end gamREGS
    idx <- idx+1
   
  } ## end boots
} ## end ldfrows

cdf0 <- cdf %>% filter(!is.na(scen)) # & !is.na(L1) & !is.na(L2)) 
cdf0 %>% write.csv(.,file = paste0('./gam_output/cdf_',Sys.Date(),'.csv'),row.names = F)
# cdf <- cdf[,c(1:5,11)]
# names(cdf)[5:6] <- c('k','Linf')

## When did regional designation go right? (original analysis)
cdfaccu <- cdf0 %>% 
  select(scen, LAT, LON) %>% 
  mutate(both = (LAT == T & LON == T)) %>%
  melt(id = c('scen')) %>%
  group_by(scen,variable) %>%
  dplyr::summarise(denom = n(), n = sum(value)) %>% 
  # dplyr::summarise(denom = n(), n = sum(gamREG == REG)) %>% 
  mutate(prop = n/denom)
# summarise(propa = sum(gamREG == REG)/n())
write.csv(cdfaccu,file = paste0('./gam_output/cdf_accu_',Sys.Date(),'.csv'),row.names = F)

## Get coverage probs -- for all params. 
## Note N varies because of varying # regions per boot. 
cdfprop <- cdf0 %>% filter(!is.na(scen))  %>% 
  # melt(id = c('scen','boot','gamREG','REG')) %>%
  select(scen, L1, L2) %>% 
  melt(id = c('scen')) %>%
  # group_by(scen) %>%
  group_by(scen,variable) %>%
  dplyr::summarise(denom = n(), n = sum(value)) %>% 
  mutate(prop = n/denom)
write.csv(cdfprop,file = paste0('./gam_output/cdf_prop_',Sys.Date(),'.csv'),row.names = F)


