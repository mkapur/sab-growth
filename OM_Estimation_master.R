## Update Apr 2020
## Re-fit @ breakpoints generated in Kapur et al but estimate sigma on a per-region basis
## using sptlvbSel2.cpp

require(mgcv);require(dplyr);require(ggplot2); require(TMB); library(reshape2)
library(gridExtra); library(grid); library(lattice)
compname <- c("Maia Kapur","mkapur")[2]
source("./functions/getGR.R")
# source(paste0(getwd(),"/functions/Deriv.R")); source("./functions/getGR.R")

compile("sptlVB_Sel_Sigma.cpp") ## pretty sure THIS is the most recent one due to link; added sigma as param_vector index via DES
dyn.load(dynlib("sptlVB_Sel_Sigma"))


## V1: In prep for meeting, refit at pub strata w sig @ strata ----

load("./sabdat_Apr2020_phase2_formatted.Rda") ## phase 2 pooling
## loads as "dat" saved up to this point since getGR is slow.
# dat <- sample_n(dat, nrow(dat)*0.25)  %>% filter(selType == 2) ## testing denom
# DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat))
DES <- ifelse(!is.na(dat$cREG), as.numeric(as.factor(dat$cREG)),1)-1 ## this is now numeric index, R3=slot 3 (idx 2)
# KEY <- paste("sab",DES,sep = "_")
temp <- data.frame(DES = as.numeric(as.character(DES)), cREG = dat$cREG); temp <- unique(temp)
keybase <- paste(as.factor(temp[order(temp$DES),'cREG']))

dat0 <- rep0 <- NULL ## later storage
nStrata <- length(unique(DES))



## this will assign a unique DES depending on period X sex X region -- whatever is in DES
data <-
  list(
    Length_cm = dat[,"Length_cm"],
    Age = dat[,"Age"],
    DES = as.vector(DES), ## keep this for master iterations
    selType = dat[,'selType'],
    Sel = dat[,'Sel'],
    nStrata = nStrata,
    a2 = 30 
  )

parameters <-
  list(
    log_Linf = rep(log(70), nStrata),
    log_k = rep(log(0.5), nStrata),
    t0 = rep(0.1, nStrata),
    log_Sigma = rep(log(2), nStrata)
  )

# Now estimate everything
map <- NULL
model <- MakeADFun(data, parameters,  DLL="it10_2",silent=T,map=map)
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
# for (k in 1:3)  fit <- nlminb(model$env$last.par.best, model$fn, model$gr) ## start at last-best call, for stability
model$report()$denominator ## if we only ran seltype 2 points, this should NOT be 1.0
best <- model$env$last.par.best
rep <- sdreport(model)

dat0 <- rep0 <- NULL ##  storage
dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
rep0 <- bind_rows(rep0,
                  bind_cols(
                    data.frame(names(rep$value)),
                    data.frame(rep$value),
                    data.frame(rep$sd),
                    data.frame(c(rep(keybase, 6)))))

## reformat outputs 
## also saving in MSE folder
names() <- c('variable', 'value','sd', 'REG')
write.csv(rep0, file = paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),'_',2,'.csv'),row.names = F) 
write.csv(rep0, file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/SAB_parEst_gam_",Sys.Date(),'.csv'),row.names = F) ## also saving in MSE folder

ypreds0 <- cbind(dat0,dat) %>% data.frame()  
names(ypreds0)[1] <- c('Predicted')
write.csv(ypreds0,  paste0("./GAM_output/SAB_predicts_",Sys.Date(),'_',2,".csv"),row.names = F)
write.csv(ypreds0, file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/SAB_predicts_",Sys.Date(),'.csv'),row.names = F) ## also saving in MSE folder

## reshape and save again
LR <- function(Linf, k, t0, AREF){
  lreg <- Linf*(1-exp(-k*(AREF - t0)))
  return(lreg)
}
rep0 %>%
  mutate(REG = as.character(REG)) %>%
  select(-sd) %>%
  # filter(variable != 'Sigma') %>%
  mutate(
    Region = substr(REG,1,2),
    Period =  sapply(strsplit(REG, "_"), function(x) x[2]),
    Sex = sub('(^[^_]+_[^_]+)_(.*
    )$', '\\2', REG) ) %>%
  tidyr::spread(key = variable, value = value) %>%
  mutate(   L1_0.5 = round(LR(Linf, k, t0, AREF = 0.5),2)) %>%
  mutate(L1 = ifelse(L1_0.5 <0, L1, L1_0.5)) %>%
  mutate(Linf =round(Linf,2), k=  round(k,2), Sigma = round(Sigma,2)) %>%
  select(Region, Sex, Period, Linf,k,Sigma) %>%
  write.csv(., file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/Table3_",Sys.Date(),".csv"),row.names = F)

require(DescTools)
parest <- read_csv("GAM_output/SAB_parEst_gam_2020-04-17_2.csv")%>%
  filter(variable ==  "k") %>%
  mutate(source = 'Estimated') %>%
  mutate(REG2 = gsub("_.*", "\\1", REG),
         REG3 = as.numeric(substr(REG,2,2)),
         REG4 = sub('_([^_]*)$', '',REG),
         Sex = gsub(".*_", "\\1", REG),
         lwr = value - 1.96 * sd,
         upr = value + 1.96 * sd,
         matchcol = 'black')
# parest$value <- exp(parest$value - (parest$sd^2)/2) ## bias correction
parest <- parest[order(parest$REG2,parest$Sex),]
for(i in 1:nrow(parest)){ ## we want to compare w adjacent regions
  if(parest$REG3 < 5){
    tmp <- subset(parest, Sex == parest$Sex[i] & REG3 == parest$REG3[i]+1)
    parest$match[i] <-ifelse(any(c(tmp$lwr, tmp$upr) %overlaps% c(parest$lwr[i], parest$upr[i])),paste0('OVERLAP ', tmp$REG2[1]),'NO OVERLAP')
  } ## end R < 5
}  

parest$pd <- substr(parest$REG4,4,nchar(parest$REG4))
parest$matchcol <- parest$match != 'NO OVERLAP'
parest$Sex <- factor(parest$Sex)
levels(parest$Sex) <- c('Females','Males')
ggplot(parest, aes(x = pd, y = value, col = REG2))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text = element_text(size = 10,angle = 45),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size=14))+
  # scale_y_continuous(limits = c(0,100)) +
  # scale_color_manual(values = c('black','red')) +
  geom_point() +
  geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
  labs(x = 'Spatiotemporal x Sex Stratum', y = "",
       col = ifelse(phase == 'phase1', 
                    'CI overlap within Region + Sex',
                    'CI Overlap Adjacent Region x Sex'),
       title = paste0(phase," K Estimates")) +
  facet_wrap(~Sex )

## V2: Luke's feedback: re-fit at orthogonal OR NPB break ----

## DON'T NEED TO RECREATE DATA UNLESS BREAK CHANGES
# load(paste0("./input_data/gam_data_sab_0415.rda")) ## full_data -- made using gam_dataprep NOT 15k subsample
# full_data$Longitude_dd[full_data$Longitude_dd > 0] <- full_data$Longitude_dd[full_data$Longitude_dd > 0]*-1
# # full_data%>% filter(REG == 'BC') %>% summarise(min(Latitude_dd), min(Longitude_dd), max(Latitude_dd), max(Longitude_dd))
# BC_ydiff =  54.34175 -   50 #48.13133 ## or 50
# BC_xdiff =  -125.7076 - -130 # -133.9758  ## or - 130
# BC_intercept = 50 - BC_ydiff/BC_xdiff *    -130
# ## quick plot of domain change
# with(subset(full_data, REG == 'BC'), 
#      plot(Latitude_dd ~ Longitude_dd,
#           ylim = c(48, 55)))
# lines(x = -134:-126, y = (BC_ydiff/BC_xdiff)*(-134:-126)+BC_intercept, lwd = 2, col = 'blue')
# 
# splitLat <- function(sample_x){
#   split_y <- (BC_ydiff/BC_xdiff)*sample_x + BC_intercept
#   return(split_y)
# }
# 
# splitBC <- function(sample_dat){
#   newReg <- ifelse(sample_dat$Latitude_dd > splitLat(sample_dat$Longitude_dd), "R4","R2")
#   return(newReg)
# }
# 
# splitBC(sample_dat = full_data[23310,]) ## should be below line (R2)
# splitBC(sample_dat = full_data[24377,]) ## should be above line (R3)
# splitBC(sample_dat = full_data[23898,]) ## should be above line (R3)
# 
# 
# breaksdf <- data.frame(yr_breaks = 2010, lat_breaks2 = c(36,50), lon_breaks2 =c(-130,-145)) ## auto
# dat <- getGR(tempdf = full_data, breaksdf) ## assign breakpoints original way
# 
# ## OVERRIDE R3 based on new split - slow
# for(i in 1:nrow(dat)){
#   if(dat$REG[i] != 'BC') next()
#   else{
#     if(i %% 100 == 0) cat(paste(i, "\n"))
#     dat$gamREG[i] <- splitBC(sample_dat = dat[i,])
#   }
# }
# 
# 
# unique(dat$gamREG)## there should be no more R3
# 
# ## check that the grouping makes sense
ggplot(subset(dat, REG == 'BC'  | REG == 'AK'), aes(x = Longitude_dd, y = Latitude_dd, col = gamREG)) +
  geom_point() #+
  # geom_abline(slope = BC_ydiff/BC_xdiff, intercept = BC_intercept, col = 'blue', lwd =1.1)
# 
# save(dat, file =  paste0("./input_data/gam_data_sab_",Sys.Date(),"_Luke.rda"))
# save(dat, file =  paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/gam_data_sab_",
#                          Sys.Date(),"_Luke.rda"))


load(paste0("./input_data/gam_data_sab_2020-05-01_Luke.rda")) ## dat made above
for(phase in c("phase1","phase2")){
  
  
  dat$selType <- ifelse(dat$REG == 'BC', 2, 1) ## 2 does truncation, 1 leaves as-is
  
  # ## add in selex for BC data, using DFO values
  selectivity <- function(Length){
    selex <- 1/(1+exp( 52.976 - Length)) ## email from Sam Johnson/
    return(selex)
  }
  dat$Sel <- ifelse(dat$REG == 'BC', selectivity(dat$Length_cm), 1)
  dat$Sel <- ifelse(dat$Sel == 0, 1E-5, dat$Sel) ## to avoid div by zero
  
  DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat))
  cat(phase," Built Dat \n")
  
  ## sanity check
  # dat %>% group_by(gamREG) %>% summarise(mnlat = mean(Latitude_dd), mnlon = mean(Longitude_dd))
  # dat %>% group_by(gamREG,REG) %>% summarise(n = n())
  
  dat$cREG <- paste0(dat$gamREG,"_",dat$Period,"_",dat$Sex)
  
  ## manual override based on Linf ambiguity
  ## phase 2

  if(phase == 'phase2'){
    ## pooling R5 females, and all except R4 males
    dat$cREG[dat$cREG == "R1_early_M"  | dat$cREG == "R1_late_M" ] <- "R1_pool_M"
    dat$cREG[dat$cREG == "R2_early_M"  | dat$cREG == "R2_late_M" ] <- "R2_pool_M"
    ## both sexes r3/R4//r5
    # dat$cREG[dat$cREG == "R3_early_M"  | dat$cREG == "R3_late_M" ] <- "R3_pool_M"
    # dat$cREG[dat$cREG == "R3_early_F"  | dat$cREG == "R3_late_F" ] <- "R3_pool_F"
    
    dat$cREG[dat$cREG == "R4_early_M"  | dat$cREG == "R4_late_M" ] <- "R4_pool_M"
    # dat$cREG[dat$cREG == "R4_early_F"  | dat$cREG == "R4_late_F" ] <- "R4_pool_F"
    
    dat$cREG[dat$cREG == "R5_early_M"  | dat$cREG == "R5_late_M" ] <- "R5_pool_M"
    dat$cREG[dat$cREG == "R5_early_F"  | dat$cREG == "R5_late_F" ] <- "R5_pool_F"
    ## females r4
    length(unique(dat$cREG)) ## 12 for phase2
  }
  
  DES <- ifelse(!is.na(dat$cREG), as.numeric(as.factor(dat$cREG)),1)-1 ## this is now numeric index, R3=slot 3 (idx 2)
  # KEY <- paste("sab",DES,sep = "_")
  temp <- data.frame(DES = as.numeric(as.character(DES)), cREG = dat$cREG); temp <- unique(temp)
  keybase <- paste(as.factor(temp[order(temp$DES),'cREG']))
  
  dat0 <- rep0 <- NULL ## later storage
  nStrata <- length(unique(DES))
  
  ## this will assign a unique DES depending on period X sex X region -- whatever is in DES
  data <-
    list(
      Length_cm = dat[,"Length_cm"],
      Age = dat[,"Age"],
      DES = as.vector(DES), ## keep this for master iterations
      selType = dat[,'selType'],
      Sel = dat[,'Sel'],
      nStrata = nStrata,
      a2 = 30 
    )
  
  parameters <-
    list(
      log_Linf = rep(log(70), nStrata),
      log_k = rep(log(0.5), nStrata),
      t0 = rep(0.1, nStrata),
      log_Sigma = rep(log(2), nStrata)
    )
  
  # Now estimate everything
  map <- NULL
  model <- MakeADFun(data, parameters,  
                     DLL="sptlVB_Sel_Sigma",
                     silent=T,map=map)
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
  # for (k in 1:3)  fit <- nlminb(model$env$last.par.best, model$fn, model$gr) ## start at last-best call, for stability
  model$report()$denominator ## if we only ran seltype 2 points, this should NOT be 1.0
  best <- model$env$last.par.best
  rep <- sdreport(model)
  
  dat0 <- rep0 <- NULL ##  storage
  dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
  rep0 <- bind_rows(rep0,
                    bind_cols(
                      data.frame(names(rep$value)),
                      data.frame(rep$value),
                      data.frame(rep$sd),
                      data.frame(c(rep(keybase, 6)))))
  
  
  ## reformat outputs 
  ## also saving in MSE folder
  names(rep0) <- c('variable', 'value','sd', 'REG')
  write.csv(rep0, file = paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),'_',phase,'.csv'),row.names = F) 
  write.csv(rep0, file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/SAB_parEst_gam_",Sys.Date(),phase,'.csv'),row.names = F) ## also saving in MSE folder
  
  ypreds0 <- cbind(dat0,dat) %>% data.frame()  
  names(ypreds0)[1] <- c('Predicted')
  write.csv(ypreds0,  paste0("./GAM_output/SAB_predicts_",Sys.Date(),'_',phase,".csv"),row.names = F)
  write.csv(ypreds0, file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/SAB_predicts_",Sys.Date(),phase,'.csv'),row.names = F) ## also saving in MSE folder
  
  ## reshape and save again
  LR <- function(Linf, k, t0, AREF){
    lreg <- Linf*(1-exp(-k*(AREF - t0)))
    return(lreg)
  }
  rep0 %>%
    mutate(REG = as.character(REG)) %>%
    select(-sd) %>%
    # filter(variable != 'Sigma') %>%
    mutate(
      Region = substr(REG,1,2),
      Period =  sapply(strsplit(REG, "_"), function(x) x[2]),
      Sex = sub('(^[^_]+_[^_]+)_(.*
      )$', '\\2', REG) ) %>%
  tidyr::spread(key = variable, value = value) %>%
    mutate(   L1_0.5 = round(LR(Linf, k, t0, AREF = 0.5),2)) %>%
    mutate(L1 = ifelse(L1_0.5 <0, L1, L1_0.5)) %>%
    mutate(Linf =round(Linf,2), k=  round(k,2), Sigma = round(Sigma,2)) %>%
    select(Region, Sex, Period, Linf,k,Sigma) %>%
    write.csv(., file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/Table3_",Sys.Date(),phase,".csv"),row.names = F)
  
 
  cat(phase," Fit TMB model  & saved outputs \n")
}


require(DescTools)
for(par in c('Linf','k')){
  parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),'_',phase,'.csv'))%>%
    filter(variable ==  par) %>%
    mutate(source = 'Estimated') %>%
    mutate(REG2 = gsub("_.*", "\\1", REG),
           REG3 = as.numeric(substr(REG,2,2)),
           REG4 = sub('_([^_]*)$', '',REG),
           Sex = gsub(".*_", "\\1", REG),
           lwr = value - 1.96 * sd,
           upr = value + 1.96 * sd,
           matchcol = 'black')
  # parest$value <- exp(parest$value - (parest$sd^2)/2) ## bias correction
  parest <- parest[order(parest$REG2,parest$Sex),]
  
  ## check if CIs overlap
  if(phase == 'phase1'){
    for(i in seq(1,nrow(parest),2)) { ## will go by sex
      parest$match[i:(i+1)] <- ifelse(c(parest$lwr[i], parest$upr[i]) %overlaps%  c(parest$lwr[i+1], parest$upr[i+1]),
                                      paste0('OVERLAP ', parest$REG2[i+1]), 'NO OVERLAP')
    }
    
  } else if(phase == 'phase2') {
    for(i in 1:nrow(parest)){ ## we want to compare w adjacent regions
      if(parest$REG3 < 5){
        tmp <- subset(parest, Sex == parest$Sex[i] & 
                        REG3 == ifelse(parest$REG3[i] == 2,4, parest$REG3[i]+1))
        parest$match[i] <-ifelse(any(c(tmp$lwr, tmp$upr) %overlaps% c(parest$lwr[i], parest$upr[i])),paste0('OVERLAP ', tmp$REG2[1]),'NO OVERLAP')
      } ## end R < 5
    } ## end parest rows
  } ## end else phase2
  
  
  parest$pd <- substr(parest$REG4,4,nchar(parest$REG4))
  parest$matchcol <- parest$match != 'NO OVERLAP'
  parest$Sex <- factor(parest$Sex)
  levels(parest$Sex) <- c('Females','Males')
  ggplot(parest, aes(x = REG4, y = value, col = matchcol))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'right',
          legend.background = element_blank(),
          axis.text = element_text(size = 10,angle = 45),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size=14))+
    # scale_y_continuous(limits = c(0,100)) +
    scale_color_manual(values = c('black','red')) +
    geom_point() +
    geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    labs(x = 'Spatiotemporal x Sex Stratum', y = "",
         col = ifelse(phase == 'phase1', 
                      'CI overlap within Region + Sex',
                      'CI Overlap Adjacent Region x Sex'),
         title = paste0(phase,par," Estimates")) +
    facet_wrap(~Sex )
  ggsave(plot = last_plot(),  
         file = paste0("./_figures/sab_parest_",par,"_",Sys.Date(),"_",phase,".png"),
         width = 10, height = 8, units = 'in', dpi = 480)
}


## V3 50 degress only ----

## DON'T NEED TO RECREATE DATA UNLESS BREAK CHANGES
load(paste0("./input_data/gam_data_sab_0415.rda")) ## full_data -- made using gam_dataprep NOT 15k subsample
full_data$Longitude_dd[full_data$Longitude_dd > 0] <- full_data$Longitude_dd[full_data$Longitude_dd > 0]*-1

splitBC <- function(sample_dat){
  newReg <- ifelse(sample_dat$Latitude_dd > 49.999, "R4","R2")
  return(newReg)
}

splitBC(sample_dat = full_data[23433,]) ## should be below line (R2)
splitBC(sample_dat = full_data[53906,]) ## should be above line (R4)



breaksdf <- data.frame(yr_breaks = 2010, lat_breaks2 = c(36,50), lon_breaks2 =c(-130,-145)) ## auto
dat <- getGR(tempdf = full_data, breaksdf) ## assign breakpoints original way
# 
# ## OVERRIDE R3 based on new split - slow
for(i in 1:nrow(dat)){
  if(dat$REG[i] != 'BC') next()
  else{
    if(i %% 100 == 0) cat(paste(i, "\n"))
    dat$gamREG[i] <- splitBC(sample_dat = dat[i,])
  }
}

# unique(dat$gamREG)## there should be no more R3

# ## check that the grouping makes sense
ggplot(subset(dat, REG == 'BC'  | REG == 'AK'), 
       aes(x = Longitude_dd, y = Latitude_dd, col = gamREG)) +
  geom_point() + geom_hline(yintercept = 50)

## save in 2 places
save(dat, file =  paste0("./input_data/gam_data_sab_",Sys.Date(),"_50N.rda"))
save(dat, file =  paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/gam_data_sab_",
                         Sys.Date(),"_50N.rda"))


load(paste0("./input_data/gam_data_sab_2020-05-06_50N.rda")) ## dat made above
for(phase in c("phase1","phase2")[2]){
  
  
  dat$selType <- ifelse(dat$REG == 'BC', 2, 1) ## 2 does truncation, 1 leaves as-is
  
  # ## add in selex for BC data, using DFO values
  selectivity <- function(Length){
    selex <- 1/(1+exp( 52.976 - Length)) ## email from Sam Johnson/
    return(selex)
  }
  dat$Sel <- ifelse(dat$REG == 'BC', selectivity(dat$Length_cm), 1)
  dat$Sel <- ifelse(dat$Sel == 0, 1E-5, dat$Sel) ## to avoid div by zero
  
  DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat))
  cat(phase," Built Dat \n")
  
  ## sanity check
  # dat %>% group_by(gamREG) %>% summarise(mnlat = mean(Latitude_dd), mnlon = mean(Longitude_dd))
  # dat %>% group_by(gamREG,REG) %>% summarise(n = n())
  
  dat$cREG <- paste0(dat$gamREG,"_",dat$Period,"_",dat$Sex)
  
  ## manual override based on Linf ambiguity
  ## phase 2
  
  if(phase == 'phase2'){
    ## pooling R5 females, and all except R4 males
    dat$cREG[dat$cREG == "R1_early_M"  | dat$cREG == "R1_late_M" ] <- "R1_pool_M"
    dat$cREG[dat$cREG == "R2_early_M"  | dat$cREG == "R2_late_M" ] <- "R2_pool_M"
    ## both sexes r3/R4//r5
    # dat$cREG[dat$cREG == "R3_early_M"  | dat$cREG == "R3_late_M" ] <- "R3_pool_M"
    # dat$cREG[dat$cREG == "R3_early_F"  | dat$cREG == "R3_late_F" ] <- "R3_pool_F"
    
    dat$cREG[dat$cREG == "R4_early_M"  | dat$cREG == "R4_late_M" ] <- "R4_pool_M"
    # dat$cREG[dat$cREG == "R4_early_F"  | dat$cREG == "R4_late_F" ] <- "R4_pool_F"
    
    dat$cREG[dat$cREG == "R5_early_M"  | dat$cREG == "R5_late_M" ] <- "R5_pool_M"
    dat$cREG[dat$cREG == "R5_early_F"  | dat$cREG == "R5_late_F" ] <- "R5_pool_F"
    ## females r4
    length(unique(dat$cREG)) ## 12 for phase2
  }
  
  DES <- ifelse(!is.na(dat$cREG), as.numeric(as.factor(dat$cREG)),1)-1 ## this is now numeric index, R3=slot 3 (idx 2)
  # KEY <- paste("sab",DES,sep = "_")
  temp <- data.frame(DES = as.numeric(as.character(DES)), cREG = dat$cREG); temp <- unique(temp)
  keybase <- paste(as.factor(temp[order(temp$DES),'cREG']))
  
  dat0 <- rep0 <- NULL ## later storage
  nStrata <- length(unique(DES))
  
  ## this will assign a unique DES depending on period X sex X region -- whatever is in DES
  data <-
    list(
      Length_cm = dat[,"Length_cm"],
      Age = dat[,"Age"],
      DES = as.vector(DES), ## keep this for master iterations
      selType = dat[,'selType'],
      Sel = dat[,'Sel'],
      nStrata = nStrata,
      a2 = 30 
    )
  
  parameters <-
    list(
      log_Linf = rep(log(70), nStrata),
      log_k = rep(log(0.5), nStrata),
      t0 = rep(0.1, nStrata),
      log_Sigma = rep(log(2), nStrata)
    )
  
  # Now estimate everything
  map <- NULL
  model <- MakeADFun(data, parameters,  
                     DLL="sptlVB_Sel_Sigma",
                     silent=T,map=map)
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
  # for (k in 1:3)  fit <- nlminb(model$env$last.par.best, model$fn, model$gr) ## start at last-best call, for stability
  model$report()$denominator ## if we only ran seltype 2 points, this should NOT be 1.0
  best <- model$env$last.par.best
  rep <- sdreport(model)
  
  dat0 <- rep0 <- NULL ##  storage
  dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
  rep0 <- bind_rows(rep0,
                    bind_cols(
                      data.frame(names(rep$value)),
                      data.frame(rep$value),
                      data.frame(rep$sd),
                      data.frame(c(rep(keybase, 6)))))
  
  
  ## reformat outputs 
  ## also saving in MSE folder
  names(rep0) <- c('variable', 'value','sd', 'REG')
  write.csv(rep0, file = paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),'_',phase,'.csv'),row.names = F) 
  write.csv(rep0, file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/SAB_parEst_gam_",Sys.Date(),phase,'.csv'),row.names = F) ## also saving in MSE folder
  
  ypreds0 <- cbind(dat0,dat) %>% data.frame()  
  names(ypreds0)[1] <- c('Predicted')
  write.csv(ypreds0,  paste0("./GAM_output/SAB_predicts_",Sys.Date(),'_',phase,".csv"),row.names = F)
  write.csv(ypreds0, file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/SAB_predicts_",Sys.Date(),phase,'.csv'),row.names = F) ## also saving in MSE folder
  
  ## reshape and save again
  LR <- function(Linf, k, t0, AREF){
    lreg <- Linf*(1-exp(-k*(AREF - t0)))
    return(lreg)
  }
  rep0 %>%
    mutate(REG = as.character(REG)) %>%
    select(-sd) %>%
    # filter(variable != 'Sigma') %>%
    mutate(
      Region = substr(REG,1,2),
      Period =  sapply(strsplit(REG, "_"), function(x) x[2]),
      Sex = sub('(^[^_]+_[^_]+)_(.*
      )$', '\\2', REG) ) %>%
  tidyr::spread(key = variable, value = value) %>%
    mutate(   L1_0.5 = round(LR(Linf, k, t0, AREF = 0.5),2)) %>%
    mutate(L1 = ifelse(L1_0.5 <0, L1, L1_0.5)) %>%
    mutate(Linf =round(Linf,2), k=  round(k,2), Sigma = round(Sigma,2)) %>%
    select(Region, Sex, Period, Linf,k,Sigma) %>%
    write.csv(., file = paste0("C:/Users/mkapur/Dropbox/UW/sab-mse/input/raw/Table3_",Sys.Date(),phase,".csv"),row.names = F)
  
  
  cat(phase," Fit TMB model  & saved outputs \n")
} ## end phase


require(DescTools)
for(par in c('Linf','k')){
  parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),'_',phase,'.csv'))%>%
    filter(variable ==  par) %>%
    mutate(source = 'Estimated') %>%
    mutate(REG2 = gsub("_.*", "\\1", REG),
           REG3 = as.numeric(substr(REG,2,2)),
           REG4 = sub('_([^_]*)$', '',REG),
           Sex = gsub(".*_", "\\1", REG),
           lwr = value - 1.96 * sd,
           upr = value + 1.96 * sd,
           matchcol = 'black')
  # parest$value <- exp(parest$value - (parest$sd^2)/2) ## bias correction
  parest <- parest[order(parest$REG2,parest$Sex),]
  
  ## check if CIs overlap
  if(phase == 'phase1'){
    for(i in seq(1,nrow(parest),2)) { ## will go by sex
      parest$match[i:(i+1)] <- ifelse(c(parest$lwr[i], parest$upr[i]) %overlaps%  c(parest$lwr[i+1], parest$upr[i+1]),
                                      paste0('OVERLAP ', parest$REG2[i+1]), 'NO OVERLAP')
    }
    
  } else if(phase == 'phase2') {
    for(i in 1:nrow(parest)){ ## we want to compare w adjacent regions
      if(parest$REG3 < 5){
        tmp <- subset(parest, Sex == parest$Sex[i] & 
                        REG3 == ifelse(parest$REG3[i] == 2,4, parest$REG3[i]+1))
        parest$match[i] <-ifelse(any(c(tmp$lwr, tmp$upr) %overlaps% c(parest$lwr[i], parest$upr[i])),paste0('OVERLAP ', tmp$REG2[1]),'NO OVERLAP')
      } ## end R < 5
    } ## end parest rows
  } ## end else phase2
  
  
  parest$pd <- substr(parest$REG4,4,nchar(parest$REG4))
  parest$matchcol <- parest$match != 'NO OVERLAP'
  parest$Sex <- factor(parest$Sex)
  levels(parest$Sex) <- c('Females','Males')
  ggplot(parest, aes(x = REG4, y = value, col = matchcol))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'right',
          legend.background = element_blank(),
          axis.text = element_text(size = 10,angle = 45),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size=14))+
    # scale_y_continuous(limits = c(0,100)) +
    scale_color_manual(values = c('black','red')) +
    geom_point() +
    geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    labs(x = 'Spatiotemporal x Sex Stratum', y = "",
         col = ifelse(phase == 'phase1', 
                      'CI overlap within Region + Sex',
                      'CI Overlap Adjacent Region x Sex'),
         title = paste0(phase,par," Estimates")) +
    facet_wrap(~Sex )
  ggsave(plot = last_plot(),  
         file = paste0("./_figures/sab_parest_",par,"_",Sys.Date(),"_",phase,".png"),
         width = 10, height = 8, units = 'in', dpi = 480)
}
