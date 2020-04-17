## Update Apr 2020
## Re-fit @ breakpoints generated in Kapur et al but estimate sigma on a per-region basis
## using sptlvbSel2.cpp

require(mgcv);require(dplyr);require(ggplot2); require(TMB); library(reshape2)
library(gridExtra); library(grid); library(lattice)
compname <- c("Maia Kapur","mkapur")[2]

source(paste0(getwd(),"/functions/Deriv.R")); source("./functions/getGR.R")

compile("it10_2.cpp") ## pretty sure THIS is the most recent one due to link; added sigma as param_vector index via DES
dyn.load(dynlib("it10_2"))


load(paste0("./input_data/gam_data_sab_0415.rda")) ## full_data -- made using gam_dataprep NOT 15k subsample
full_data$Longitude_dd[full_data$Longitude_dd > 0] <- full_data$Longitude_dd[full_data$Longitude_dd > 0]*-1


# breaksdf <- data.frame(yr_breaks = 2010, lat_breaks2 = c(36,50), lon_breaks2 =c(-130,-145)) ## auto
# dat <- getGR(tempdf = full_data, breaksdf) ## assign breakpoints
# dat$selType <- ifelse(dat$REG == 'BC', 2, 1) ## 2 does truncation, 1 leaves as-is
# 
# ## add in selex for BC data, using DFO values
# selectivity <- function(Length){
#   selex <- 1/(1+exp( 52.976 - Length)) ## email from Sam Johnson/
#   return(selex)
# }
# dat$Sel <- ifelse(dat$REG == 'BC', selectivity(dat$Length_cm), 1)
# # dat$Sel <- ifelse(dat$Sel == 0, 1E-5, dat$Sel) ## to avoid div by zero
# 
# DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat)) 
# cat(phase," Built Dat \n")
# 
# ## sanity check
# # dat %>% group_by(gamREG) %>% summarise(mnlat = mean(Latitude_dd), mnlon = mean(Longitude_dd))
# # dat %>% group_by(gamREG,REG) %>% summarise(n = n())
# 
# dat$cREG <- paste0(dat$gamREG,"_",dat$Period,"_",dat$Sex)
# 
# ## manual override based on Linf ambiguity ----
# ## phase 2 
# ## males r1/r2
# if(phase == 'phase2'){
#   dat$cREG[dat$cREG == "R1_early_M"  | dat$cREG == "R1_late_M" ] <- "R1_pool_M"
#   dat$cREG[dat$cREG == "R2_early_M"  | dat$cREG == "R2_late_M" ] <- "R2_pool_M"
#   ## both sexes r3/R4//r5
#   dat$cREG[dat$cREG == "R3_early_M"  | dat$cREG == "R3_late_M" ] <- "R3_pool_M"
#   # dat$cREG[dat$cREG == "R3_early_F"  | dat$cREG == "R3_late_F" ] <- "R3_pool_F"
#   
#   dat$cREG[dat$cREG == "R4_early_M"  | dat$cREG == "R4_late_M" ] <- "R4_pool_M"
#   # dat$cREG[dat$cREG == "R4_early_F"  | dat$cREG == "R4_late_F" ] <- "R4_pool_F"
#   
#   dat$cREG[dat$cREG == "R5_early_M"  | dat$cREG == "R5_late_M" ] <- "R5_pool_M"
#   dat$cREG[dat$cREG == "R5_early_F"  | dat$cREG == "R5_late_F" ] <- "R5_pool_F"
#   ## females r4
#   length(unique(dat$cREG)) ## 12
# }

## MK START HERE
load("./sabdat_Oct2019_formatted.Rda")
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
    a2 = 30 ## wtf is this supposed to be
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
                    data.frame(c(rep(keybase, 5), rep("ALL", 1)))))

## reformat outputs ----
names(rep0) <- c('variable', 'value','sd', 'REG')
write.csv(rep0, file = paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),'_',phase,'.csv'),row.names = F)

ypreds0 <- cbind(dat0,dat) %>% data.frame()  
names(ypreds0)[1] <- c('Predicted')
write.csv(ypreds0,  paste0("./GAM_output/SAB_predicts_",Sys.Date(),'_',phase,".csv"),row.names = F)

cat(phase," Fit TMB model  & saved outputs \n")