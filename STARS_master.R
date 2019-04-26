## Comparison with STARS method a la Rodionov 2003. 

require(dplyr);require(ggplot2); require(TMB); library(reshape2)
compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))
compname <- c("Maia Kapur","mkapur")[1]

scenarios <- read.csv("./input_data/scenarios.csv",na.strings = 'NA') ## manual file
nboot <- 100

## Create STARS breakpoints (only need to run once) ----
# source("./functions/runSTARS.R") ## slightly mod stars functionality
# rm(breaksdf);  idx <- 1 ## storage coverage prob totals, rbind each scen
# breaksdf <- data.frame()
# Terms <- c("Year","Latitude_dd","Longitude_dd")

# for(l in 1:length(unique(scenarios$DESC))){ 
#   scen <- unique(scenarios$DESC)[l]
#   for(b in 1:nboot){
#     cat(paste0(scen," boot ",b,"\n"))
#     tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_",b,".csv")) %>% filter(Age == 6)
#     if(scen == 'NoBreaks') tempdf$REG <- as.factor('R1')
#     breaksdf[idx,'scen'] = scen; breaksdf[idx,'boot'] = b
#     
#     for(t in 1:length(Terms)){
#       tempstar <- stars(y= tempdf$Length_cm, time  = tempdf[,Terms[t]], L=5, p=0.05, h=1,  AR1red="none", prewhitening = FALSE)
#       
#       if(is.na(tempstar)  ){ breaksdf[idx,Terms[t]] <- NA ; next()
#       } else {
#         tempstar$starsResult[which(tempstar$starsResult[,3]==0),3] = NA
#         if(all(is.na(tempstar$starsResult[,3]))) { breaksdf[idx,Terms[t]] <- NA ; next()
#         } else{ breaksdf[idx,Terms[t]] <- round(as.numeric(rownames(tempstar$starsResult)[which.max(abs(tempstar$starsResult[,3]))])) }
#         # breaksdf[idx,Terms[t]] <- ifelse(is.na(tempstar),NA,which.max(abs(tempstar$starsResult[,3]))[1])
#       }## end else
#     } ## end loop terms
#     idx <- idx+1 ## bump down for next boot
#   } ## end boots
# } ## end scen
# write.csv(breaksdf, paste0("./GAM_output/STARS_breaksdf.csv"),row.names=F)

## Read in and re-fit based on STARS breaks ----
source("./functions/getGR.R");source("./functions/fitMod.R");source("./functions/makematchdf.R")
cdf <- data.frame(); idx <-1 ## storage coverage prob totals, rbind each scen
for(l in 1:length(unique(scenarios$DESC))){ 
  for(b in 1:nboot){
    ## TMB FITTING ----
    scen <- unique(scenarios$DESC)[l]
    tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_",b,".csv")) #%>% filter(Age == 6)
    if(scen == 'NoBreaks') tempdf$REG <- as.factor('R1')
    breaksdf <- read.csv("./GAM_output/STARS_breaksdf.csv") %>% filter(scen ==  unique(scenarios$DESC)[l] & boot == b)
    names(breaksdf)[3:5] <- c("yr_breaks","lat_breaks2","lon_breaks2")
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
        nStrata = nStrata
      )
    
    parameters <-
      list(
        log_Linf = rep(log(150), nStrata),
        log_k = rep(log(0.3), nStrata),
        t0 = rep(0.1, nStrata),
        log_Sigma = 0
      )
    
    # Now estimate everything
    map <- NULL
    mod <- fitMod(data, parameters, modversion = 'sptlvb',map)
    dat0 <- mod[[1]]
    rep0 <- mod[[2]]
    
    
    ## reformat outputs ----
    write.csv(rep0, file = paste0("./output_data/STARS_",scen,"_parEst_gam_",b,"_",Sys.Date(),'.csv'),row.names = F)
    ypreds0 <- data.frame(dat0); names(ypreds0)[1] <- c('Predicted')
    write.csv(ypreds0,  paste0("./output_data/STARS_",scen,"_predicts_",b,Sys.Date(),".csv"),row.names = F)
    
    cat("Fit TMB model ",paste(scen)," boot ",b," & saved outputs \n")
    
    ## check and save coverage ----
    # vars <- paste0('R',1:4,"_");   vis <- c("pooled",'early','late');
    # Rlevs <- c("R1", "R2","R3","R4") #,apply(expand.grid(vars, vis), 1, paste, collapse=""))
    vars <- paste0('R',1:4,"_");   vis <- c("pooled",'early','late'); Rlevs <- c("R1", "R2","R3","R4",apply(expand.grid(vars, vis), 1, paste, collapse=""))
    source("./functions/compareBreaks2.R")
    
  } ## end boots
} ## end ldfrows
 	
 	

cdf %>% write.csv(.,file = paste0('./stars_output/STARS_cdf_',Sys.Date(),'.csv'),row.names = F)

## Get coverage probs -- for all params. 
## Note N varies because of varying # regions per boot. 
cdfprop <- cdf %>%
  filter(!is.na(scen) & !is.na(L2))  %>% 
  select(scen, L1, L2) %>% 
  mutate(both = (L1 == T & L2 == T)) %>%
  melt(id = c('scen')) %>%
  group_by(scen,variable) %>%
  dplyr::summarise(denom = n(), n = sum(value)) %>% 
  mutate(prop = round(n/denom,2))
write.csv(cdfprop,file = paste0('./stars_output/STARS_cdf_prop_',Sys.Date(),'.csv'),row.names = F)

## When did regional designation go right? (original analysis)
cdfaccu <- cdf %>% 
  select(scen, LAT, LON, YEAR) %>%
  mutate(both = (LAT == T & LON == T), all = (both == T & YEAR == T)) %>%
  melt(id = c('scen')) %>%
  group_by(scen,variable) %>%
  dplyr::summarise(denom = n(), n = sum(value)) %>% 
  mutate(prop = round(n/denom,2))
write.csv(cdfaccu,file = paste0('./stars_output/STARS_cdf_accu_',Sys.Date(),'.csv'),row.names = F)

	
		
	