## Code to test GAM method on datasets generated via iPopSim (IBM)
## This will fit a gam for each replicate, identify break points, fit new life history params
## and plot comparisons
## Maia Sosa Kapur kapurm@uw.edu WiSp 2019

require(mgcv);require(dplyr);require(ggplot2); require(TMB)
library(gridExtra); library(grid); library(lattice)

compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))

# comparedf<- expand.grid(paste0("R",c(1:4)),paste0("R",c(1:4))) %>% mutate('compare' = ifelse(Var1 == Var2 | 
#                                                                                    Var1 == 'R2' & Var2 == 'R3'  |
#                                                                                    Var1 == 'R3' & Var2 == 'R2', TRUE,FALSE)) %>% 
#   filter(!(Var1 == 'R4'))
# names(comparedf)[1:2] <- c('Actual','Estimated')
compname <- c("Maia Kapur","mkapur")[1]

## Build mods, get breakpoints
scenarios <- read.csv("./input_data/scenarios.csv",na.strings = 'NA') ## manual file
age <- 6; nboot <- 100; testrows <- unique(scenarios$DESC) 
# source("./functions/bootBreaks.R") ## about 10 mins -- don't need to do this >1x
source("./functions/getGR.R")
# now loop boots and fit VB to get coverage probs
ldfprop <-  read.csv( paste0("./GAM_output/ldf_raw_a6.csv"))
cdf <- data.frame();idx <-1 ## storage coverage prob totals, rbind each scen

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
    dat0 <- rep0 <- NULL ## later storage
    nStrata <- length(unique(DES))
    
    data <-
      list(
        Length_cm = dat[,"Length_cm"],
        Age = dat[,"Age"],
        DES = as.vector(DES),
        nStrata = nStrata
      )
    
    parameters <-
      list(
        log_Linf = rep(log(200), nStrata),
        log_k = rep(log(0.3), nStrata),
        t0 = rep(0.1, nStrata),
        log_Sigma = 0
      )
    
    # Now estimate everything -- DO PHASES
    map <- NULL
    # map <- list(log_Linf = factor(NA)) ## get K, other stuff shouldn't change much
    model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)
    fit <- nlminb(
      model$par,
      model$fn,
      model$gr,
      control = list(
        rel.tol = 1e-12,
        eval.max = 100000,
        iter.max = 10000
      )
    )
    # best <- model$env$last.par.best
    # rep1 <- cbind(names(sdreport(model)$value)[2],
    #               data.frame(sdreport(model)$value['log_k'], data.frame(sdreport(model)$sd[2])))
    # names(rep1) <- c('variable', 'value', 'sd')
    
    # 
    # map <- list(log_k = factor(NA))
    # model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)
    # fit <- nlminb(
    #   model$par,
    #   model$fn,
    #   model$gr,
    #   control = list(
    #     rel.tol = 1e-12,
    #     eval.max = 100000,
    #     iter.max = 10000
    #   )
    # )
    # 
    best <- model$env$last.par.best
    rep <- sdreport(model)
    
    dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
    
    rep0 <- bind_rows(rep0,
                      bind_cols(
                        data.frame(names(rep$value)),
                        data.frame(rep$value),
                        data.frame(rep$sd),
                        
                        data.frame(c(rep(keybase
                                         , 3), rep("ALL", 1)))
                        
                      ))
    names(rep0) <- c('variable', 'value','sd', 'REG')
    
    rep0$value[rep0$variable %in% c('log_k','log_Linf')] <- exp(rep0$value[rep0$variable %in% c('log_k','log_Linf')])
    rep0$sd[rep0$variable %in% c('log_k','log_Linf')] <- exp(rep0$sd[rep0$variable %in% c('log_k','log_Linf')]) 
    rep0$variable <- factor(rep0$variable, levels = c("k","Linf","Sigma","t0","log_k","log_Linf")) ## enable new levels
    rep0$variable[rep0$variable == 'log_k'] <- 'k'
    rep0$variable[rep0$variable == 'log_Linf'] <- 'Linf'
    
    # rep3<-   bind_cols(
    #   data.frame(names(rep2$value[c(1,3,4)])),
    #   data.frame(rep2$value[c(1,3,4)]),
    #   data.frame(rep2$sd[c(1,3,4)]))
    # names(rep3) <- c('variable', 'value','sd')
    # 
    # rep0 <- rbind(rep0, cbind(rbind(rep1,rep3),data.frame(c(rep(keybase
    #                                                             , 3), rep("ALL", 1)))))
                  
                  # ## reformat outputs ----
    write.csv(rep0, file = paste0("./output_data/",scen,"_parEst_gam_",b,"_",Sys.Date(),'.csv'),row.names = F)

    ypreds0 <- cbind(dat0,dat) %>% data.frame()
    names(ypreds0)[1] <- c('Predicted')

    write.csv(ypreds0,  paste0("./output_data/",scen,"_predicts_",b,Sys.Date(),".csv"),row.names = F)
    ypreds <-  read.csv(paste0("./output_data/",scen,"_predicts_",b,Sys.Date(),".csv"))
    
    cat("Fit TMB model ",paste(scen)," boot ",b," & saved outputs \n")

    ## check and save coverage --
    ## plot estimates
    gamREGS <- factor(unique(dat$gamREG), levels=c("R1", "R2","R3","R4"))
    REGS <-  factor(unique(dat$REG), levels=c("R1", "R2","R3","R4")) ## what was used to GENERATE
    
    parest <-  read.csv(paste0("./output_data/",scen,"_parEst_gam_",b,"_",Sys.Date(),'.csv')) %>%
      filter(variable != "Sigma" & variable != 't0') %>% mutate(source = 'Estimated')
    parest <- rbind(parest, (read.csv("./input_data/true_ibm_vals.csv") %>% filter(REG %in% REGS)))    ## only extract regions from 'actual' that were present in original

    parest$REG <- factor(parest$REG, levels=c("R1", "R2","R3","R4"))
    parest$lwr <- parest$value - 1.96*parest$sd
    parest$upr <- parest$value + 1.96*parest$sd

    parms <- c('k','Linf')
   
    for(iq in 1:length(gamREGS)){ ## loop regions
      cdf[idx,'scen'] <- scen
      cdf[idx,'boot'] <- b
      cdf[idx,'gamREG'] <- gamREGS[iq]
      cdf[idx,'REG'] <- if(!is.na(REGS[iq])){REGS[iq]}else{last(REGS)}
      ptf <- NULL
      for(v in  1:length(parms)){ ## loop variables
        tmp <- subset(parest, variable == parms[v])
        bounds <- tmp %>% filter(source == 'Estimated' & REG == gamREGS[iq]) %>% select(lwr,upr)
        # cdf[idx,'k'] <- ptf[1];  cdf[idx,'Linf'] <- ptf[2]
        ptf[v] <- tmp$value[tmp$source == 'Actual' &  tmp$REG ==    cdf[idx,'REG']] >= bounds[1] &
          tmp$value[tmp$source == 'Actual'&  tmp$REG ==    cdf[idx,'REG']] <= bounds[2]
      } ## end vars
      cdf[idx,'k'] <- ptf[1];  cdf[idx,'Linf'] <- ptf[2]
      idx <- idx+1
    } ## end gamREGS
    idx <- idx+1
    # cdf <- rbind(cdf, cdf)
    # linfs <- subset(parest, variable == 'Linf')%>% filter(compare == TRUE)
    # ks <- subset(parest, variable == 'k')%>% filter(compare == TRUE)
    # ## loop regions
    #   
    #   cdf[idx,'scen'] <- scen
    #   
    #   # cdf[idx,'quad'] <- keybase[q]
    #   
    #   cdf[,'gamREG'] <- factor( cdf[,'gamREG'], levels=c("R1", "R2","R3","R4"))
    #   
    #   cdf[,'REG'] <- factor( cdf[,'REG'], levels=c("R1", "R2","R3","R4"))
    #   
    #  linf.actual <- merge(linfs,comparedf)  %>% filter( Actual ==  cdf[idx,'REG'] & Estimated ==  cdf[idx,'gamREG'] & source == 'Actual') %>% select(value)
    #   linfTF <- linf.actual %
    #   
    #   if(subset(comparedf, Estimated == cdf[idx,'gamREG'] & Actual == cdf[idx,'REG'])$compare){
    #     linfTF <- subset(parest, variable == 'Linf' & REG == cdf[idx,'REG'] & source == 'Actual')$value %in%
    #       subset(parest, variable == 'Linf' & REG == cdf[idx,'gamREG'] & source == 'Estimated')$value
    #       
    #       }
    # 
    #   idx <- idx+1
    # }
    
    # cdf <- rbind(cdf, cdf)
  } ## end boots
} ## end ldfrows

cdf
# cdf <- cdf[,c(1:5,11)]
# names(cdf)[5:6] <- c('k','Linf')

cdf %>% melt(id = c('scen','boot','gamREG','REG')) %>% group_by(scen) %>% dplyr::summarise(n = sum(value)) %>% mutate(prop = n/nrow(.))

