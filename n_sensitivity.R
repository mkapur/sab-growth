## sensitivity to sample size
## just...run the whole simulation again with 3/4 and half samples

require(mgcv)
## wrapper function for getBreaks and makeMod, given Age
source("./functions/getBreaks.R")
require(mgcv);require(dplyr);require(ggplot2); require(TMB); library(reshape2)
library(gridExtra); library(grid); library(lattice)

compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))
compname <- c("Maia Kapur","mkapur")[1]

## Build mods, get breakpoints
scenarios <- read.csv("./input_data/scenarios.csv",na.strings = 'NA')## manual file
age <- 6; nboot <- 100; testrows <- unique(scenarios$DESC)

## bootBreaks.R modified ----
## loop through boots, build GAMS and ID Breakpoints ----
for(size in c('half','threequ')){
  ldf <- data.frame()
  for(l in 1:length(testrows)){
    
    scen0 <- testrows[l]
    scen <- scen0
    outdir0 <-  paste0("C:/users/",compname,"/dropbox/uw/sab-growth/gam/plots/", 
                       scen)
    if(!exists(outdir0)) dir.create(outdir0)
    
    for(b in 1:nboot){ ## loop boots
      dat0 <- read.csv(paste0("./IBM_output/datasets/",scen,"_",b,'.csv'))  %>% filter(Age == age)
      if(size == 'half'){ ## make vector of unique rows
        temprows <- sample(1:nrow(dat0),nrow(dat0)/2)
      } else if(size == 'threequ'){
        temprows <- sample(1:nrow(dat0),nrow(dat0)*0.75)
      }
      dat <- dat0[temprows,]
      
      cat( testrows[l],' boot ',b,' Age ',age,' sample size ',size,'\n')
      
      mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd), data = dat) 
      save(mod, file = paste0("./GAM_output/mod_",scen,"_",b,size,".rds")) ## save on the fly
      bdf <-  getBreaks(mod,dat, scen)
      
      ## fill NAs in bdf for binding
      for(i in 1:length(bdf)){
        if (length(bdf[[i]]) == 0){ ## fill NA for empty
          bdf[[i]] <- NA
        }
      }
      
      tldf <- cbind(as.numeric(bdf[[1]]), as.numeric(bdf[[2]]), as.numeric(bdf[[3]]), as.character(rep(scen, length(bdf[[1]]))), rep(b, length(bdf[[1]])))
      
      ldf <- rbind(ldf, tldf)
      
    } ## end boots
    
  } ## end scenarios
  ## save
  names(ldf) <- c('yr_breaks','lat_breaks','lon_breaks','scen','boot')
  ldf$lat_breaks2 <- factor(ldf$lat_breaks, levels=c(paste(1:50),NA))
  ldf$lon_breaks2 <- factor(ldf$lon_breaks, levels=c(paste(1:50),NA))
  write.csv(ldf, file = paste0("./GAM_output/ldf_raw_a",age,"_",size,".csv"),row.names = F)
  cat('plotted & saved break tabulations in GAM_output \n')
}
source("./functions/getGR.R");source("./functions/fitMod.R");source("./functions/missby.R")
for(size in c('half','threequ')){
  
  # now loop boots and fit VB to get coverage probs
  ldfprop <-  read.csv( paste0("./GAM_output/ldf_raw_a6","_",size,".csv")) 
  rm(cdf); cdf <- data.frame(); idx <-1 ## storage coverage prob totals, rbind each scen
  for(l in 1:length(unique(ldfprop$scen))){
    # for(b in sample(1:nboot,5)){
    for(b in 1:nboot){
      ## TMB FITTING ----
      scen <- unique(ldfprop$scen)[l]
      tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_",b,".csv")) ## ok to do full dataset
      if(scen == 'NoBreaks') tempdf$REG <- as.factor('R1')
      breaksdf <- read.csv( paste0("./GAM_output/ldf_raw_a",age,"_",size,".csv")) %>% filter(scen ==  unique(ldfprop$scen)[l] & boot == b)
     breaksdf <- breaksdf[1,]
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
          a2 = 15
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
      write.csv(rep0, file = paste0("./output_data/",scen,"_parEst_gam_",b,"_",Sys.Date(),size,'.csv'),row.names = F)
      ypreds0 <- data.frame(dat0); names(ypreds0)[1] <- c('Predicted')
      write.csv(ypreds0,  paste0("./output_data/",scen,"_predicts_",b,Sys.Date(),size,".csv"),row.names = F)
      cat("Fit TMB model ",paste(scen)," GAM boot ",b," & saved outputs \n")
      
      ## check and save coverage ----
      parms <- c('L1','L2','k','Linf')[1:2]
      vars <- paste0('R',1:4,"_");   vis <- c("pooled",'early','late'); Rlevs <- c("R1", "R2","R3","R4",apply(expand.grid(vars, vis), 1, paste, collapse=""))
      source("./functions/compareBreaks2.R")
      cat("Tabulated Breaks ",paste(scen)," GAM boot ",b," CDF row ",idx," \n")
      
    } ## end boots
  } ## end ldfrows
  
  cdf %>%
    mutate(method = paste0('GAM_',size)) %>%
    write.csv(.,file = paste0('./gam_output/cdf_',Sys.Date(),size,'.csv'),row.names = F)
  
  
  ## Get coverage probs -- for all params. 
  ## Note N varies because of varying # regions per boot. 
  cdfprop <- cdf %>%
    filter(!is.na(scen) & !is.na(L2))  %>% 
    select(scen, L1, L2) %>% 
    mutate(both = (L1 == T & L2 == T)) %>%
    melt(id = c('scen')) %>%
    group_by(scen,variable) %>%
    dplyr::summarise(denom = n(), n = sum(value)) %>% 
    mutate(prop = round(n/denom,2))%>%
    mutate(method = "GAM")
  write.csv(cdfprop,file = paste0('./gam_output/cdf_prop_',Sys.Date(),size,'.csv'),row.names = F)
  
  ## When did regional designation go right? (original analysis)
  cdfaccu <- cdf %>% 
    select(scen, LAT, LON, YEAR) %>%
    mutate(both = (LAT == T & LON == T), all = (both == T & YEAR == T)) %>%
    melt(id = c('scen')) %>%
    group_by(scen,variable) %>%
    dplyr::summarise(denom = n(), n = sum(value)) %>% 
    mutate(prop = round(n/denom,2))%>%
    mutate(method = "GAM")
  write.csv(cdfaccu,file = paste0('./gam_output/cdf_accu_',Sys.Date(),size,'.csv'),row.names = F)
  
  
} ## end size loop 2

for(size in c('half','threequ')){
  
  
  
  ## plotmaster
  ## GAM propagg----
  cdfprop <- read.csv(paste0('./gam_output/cdf_prop_',Sys.Date(),size,'.csv')) 
  levels(cdfprop$scen) <- c("Break at 25 deg.", "Break at 49 deg.",
                            "20% Higher k, Break at 25 deg.", "Overlap 20-25 deg.","No Breaks","Temporal Break Year 50")
  
  # cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
  levels(cdfprop$variable) <- c('Both L1 and L2','L1','L2' )
  cdfprop2 <- cdfprop %>% filter(!(variable %in% c('Both L1 and L2')) & scen != "20% Higher k, Break at 25 deg.")
  
  # cdfprop$variable <- factor(cdfprop$variable, levels=c('L1','L2','Both L1 and L2'))
  plist1 <- list()
  plist1[[1]] <- ggplot(cdfprop2, aes(x = scen, y = prop, fill = scen)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(),
          legend.position = c(0.9,0.75)) +
    scale_fill_brewer(palette = ifelse(size == 'half','Blues','Reds'))+
    scale_y_continuous(limits = c(0,1)) +
    labs(x = '',y = 'Coverage Probability', fill = 'Scenario', 
         title = 'a) Coverage Probability for Endpoints of Growth Curve') +
    geom_bar(stat = 'identity',width=0.6, position = position_dodge(width=0.7)) +
    facet_wrap(~variable)
  # ggsave(plot = last_plot(),  file = paste0("./figures/cdfprop.png"), width = 9, height = 6, units = 'in', dpi = 480)
  
  cdfaccu <- read.csv(paste0('./gam_output/cdf_accu_',Sys.Date(),size,'.csv'))
  levels(cdfaccu$scen) <-c("Break at 25 deg.", "Break at 49 deg.",
                           "20% Higher k, Break at 25 deg.", "Overlap 20-25 deg.","No Breaks","Temporal Break Year 50")
  # cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
  levels(cdfaccu$variable) <- c('Lat, Long and Year','Both Latitude and Longitude','Latitude', 'Longitude' ,'Year')
  # cdfaccu$variable <- factor(cdfaccu$variable, levels=c('L1','L2','Both L1 and L2'))
  
  cdfaccu2 <- cdfaccu %>% filter(!(variable %in% c('Lat, Long and Year','Both Latitude and Longitude'))& scen != "20% Higher k, Break at 25 deg.")
  # cdfaccu$scen  <- factor(cdfaccu$scen , levels = cdfaccu$scen [order(cdfprop$prop  )])
  plist1[[2]] <- ggplot(cdfaccu2, aes(x = scen, y = prop, fill = scen)) +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_blank(),
          legend.position = 'none') +
    scale_fill_brewer(palette = ifelse(size == 'half','Blues','Reds'))+
    scale_y_continuous(limits = c(0,1)) +
    labs(x = '',y = 'Proportion Detected Accurate Spatial Breaks', fill = 'Scenario', title = 'b) Proportion Detected Accurate Breaks') +
    geom_bar(stat = 'identity',width=0.5, position = position_dodge(width=0.5)) +
    facet_wrap(~variable, nrow = 1)
  
  ggarrange(plotlist = plist1, ncol=1, nrow=2, common.legend = TRUE, legend="bottom") %>%
    ggsave(plot = .,  file = paste0("./figures/GAM_cdfprob_",Sys.Date(),size,".png"), width = 11, height = 8, units = 'in', dpi = 480)
  
} #end size loop 3
