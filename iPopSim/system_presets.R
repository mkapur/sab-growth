## Unchanging System Presets -- Universal for SpatioTemporal Regimes
## kapurm@uw.edu
## System Functions and Presets
## iPopSim MK
## 2017 Summer
## kapurm@uw.edu
options(warn = -1)
pow=function(x,y){x^y}

## GENERAL IBM PARAMETERS ----
# setup the calendar year
start_yr = 0

# num of M only year 
M_only_yr = 50

# Fishing mortality starting year
# F_start_yr = 60
F_start_yr = 50

simu_year = M_only_yr+F_start_yr

# ## GROWTH MODULE ---- MOVED TO REGION SPECIFIC FILE
# # length-weight
# lw_a = 1.3528e-006
# lw_b = 3.4297
# 
# # SS3 growth model (vonBert with L1&L2, SS3 3.34f mannual, Page 42)
# L1 = 62.69
# L2 = 216.72
# a1 = 0
# a2 = 15
# growth_lognormal_SD = 0.025
# VBGF_K = 0.258

Linf =L1+(L2-L1)/(1-exp(-VBGF_K*(a2-a1))) 

## Growth increment FUNCTION, NON TV: FROM V7

growth_incre  = function(fish_size){ 
  y=(Linf-fish_size)*(1-exp(-VBGF_K))
  y
}

# length-weight FUNCTION
lw  = function(fish_size){
  y=lw_a*fish_size^lw_b
  y
}  

## RECRUITMENT/REPRO MODULE ----
# Maturity ogive
r = -0.1034
L50 = 143.68

# Maturity ogive FUNCTION
prob_mature = function(fish_size){
  y = 1/(1+exp(r*(fish_size-L50)))
  y
}

# num of super individuals - max recruits per year
R0_super = 12
# R0_super = 3000

## Beverton-Holt SRR 
sigma_R = 0.1
steepBH = 0.9 ## bev holt steepness

# B-H S-R FUNCTION                                              
BH_SR = function(SSBy){                                         
  y = (4*steepBH*R0_super*SSBy/(SSB0*(1-steepBH)+SSBy*(5*steepBH-1)))   
  y                                                               
}                                                               

## LFSR

## Explicit derivation of LFSR parameters
beta = 2 
sfrac = 0.391
z0 = ## neg log of pre-recruit mortality at unfished equilibrium, log(r0/pups0)
  calc.LFSRh = function(beta,z0,zfrac){
    h = 0.2*exp(z0*zfrac*(1-pow(0.2,beta)))
    return(h)
  }
LF_SR = function(SSBy,beta,z0,zfrac,zmin){
  y = exp(-z0+(z0-zmin)*(1-pow(SBY/SSB0,rho)))*SBY*exp()
  y
}

## MORTALITY MODULE --==
# mortality setup 
Inst_M = 0.25

# catch data sampling rate
#catch_sampling_rate = 0.9 

# selectivity ogive
# SEL_50 = 147 
# SEL_95 = 176
## ages -- eyeballed from plot!
SEL_50 <- 2
SEL_95 <- 4

# fishery-dependent catchability
fishery_q = 0.00001
sigma_CPUE = 0.1

# Age composition
# Ageing error definition 
# Mean age' for each true age
# Caution: need to define a full age range
max_life_chance = 200
mean_age_true_age = c(0:max_life_chance)+0.5

# SD age' for each true age
# Caution: need to define a full age range
mean_SD_age_true_age = rep(0.001,length(mean_age_true_age))
Effect_N_Age_COM = 20

# Caution: "No need" to define a full age range 
OBS_age_bin = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)

# Mean-size-at-age
Nsamp_Mean_length_at_age = 20

# selectivity FUNCTION
# selectivity  = function(fish_size){
#   y=1/(1+exp(-log(19)*(fish_size-SEL_50)/(SEL_95-SEL_50)))
#   y
# }
selectivity  = function(true_age){
  y=1/(1+exp(-log(19)*(true_age-SEL_50)/(SEL_95-SEL_50)))
  y
}

test_fish_size = c(0:400)

## Effort Simulator ----
generate.effort = function(level, nyears) {
  getFhist<-function(nsim,Esd,nyears,dFmin,dFmax,bb,scale){
    ne<-nsim*10                                           # Number of simulated effort datasets
    dEfinal<-runif(ne,dFmin,dFmax)                        # Sample the final gradient in effort
    a<-(dEfinal-bb)/nyears                                # Derive slope to get there from intercept
    a<-array(a,dim=c(ne,nyears))                          # Slope array
    bb<-array(bb,dim=c(ne,nyears))                        # Intercept array
    x<-array(rep(1:nyears,each=ne),dim=c(ne,nyears))      # Year array
    dE<-a*x+bb                                            # Change in effort
    E<-array(NA,dim=c(ne,nyears))                         # Define total effort array
    E[,1]<-dE[,1]
    for(y in 2:nyears){
      E[,y]<-apply(dE[,1:y],1,sum)
    }
    
    E<-scale*E/array(apply(E,1,mean),dim=c(ne,nyears))          # Standardise Effort to average 1
    
    cond<-apply(E,1,min)>0
    pos<-(1:ne)[cond]
    pos<-pos[1:nsim]
    
    E<-E[pos,]                                            # Sample only those without negative effort
    Emu<--0.5*Esd^2
    Eerr<-array(exp(rnorm(nyears*nsim,rep(Emu,nyears),rep(Esd,nyears))),c(nsim,nyears))
    E*Eerr
  }
  
  ## just run a single vector and save it to an ongoing data frame
  scaleL <- ifelse(level == 'HIGH', 0.3/1.323, ifelse(level == 'MED', 0.2/1.323, 0.1/1.323))
  # effvec <- as.vector(getFhist(nsim=1,Esd=0.001,nyears=F_start_yr+1,dFmin=-0.07,dFmax=0.07,bb=-0.1,scale = scaleL))
  ## create a data frame based on your number of sims and desired fishing effort. Add in an extra year.

  if (level == 'HIGH') {
    effvec <-   runif(nyears, 0.3, 0.4)
  } else if (level == 'MED') {
    effvec <-     runif(nyears, 0.2, 0.3)
  } else{
    effvec <-     runif(nyears, 0.1, 0.2)
  }
  return(effvec)
}
#   df <- data.frame(
#     'SIM' = rep(1:nsim, each = F_start_yr+1),
#     'YEAR' = rep(seq(1, F_start_yr+1, 1), nsim),
#     'FMORT' = if (level == 'HIGH') {
#       as.vector(getFhist(nsim=1,Esd=0.001,nyears=F_start_yr+1,dFmin=-0.07,dFmax=0.07,bb=-0.1,scale = 0.3/1.323))
#     }
#     else if (level == 'MED') {
#       as.vector(getFhist(nsim=1,Esd=0.001,nyears=F_start_yr+1,dFmin=-0.07,dFmax=0.07,bb=-0.1,scale = 0.2/1.323))
#     } else {
#       as.vector(getFhist(nsim=1,Esd=0.001,nyears=F_start_yr+1,dFmin=-0.07,dFmax=0.07,bb=-0.1,scale = 0.1/1.323))
#       
#     }
#   )
#   write.csv(df,
#             paste0(getwd(),'/inputs/FSIM_', level, '.csv'),
#             row.names = F)
# }




## Bookkeeping - set up empty dataframes for output ----
maketables = function(path_name){
  SAA = data.frame(cohort=c(-999),uni_ind = c(-999),ind =c(-999),
                   fish_size=c(-999),Year=(-999),Age=(-999),AgeE = c(-999),
                   FMORT=c(-999),M=c(-999),Sel=c(-999), FSIM = c(-999), REG = c(-999))
  
  dinfo = data.frame(cohort=c(-999),ind =c(-999),fish_size=c(-999),
                     Year=c(-999),Age=c(-999),AgeE = c(-999),dtype="HEAD",
                     FMORT=c(-999),M=c(-999),Sel=c(-999), FSIM = c(-999), REG = c(-999))
  eff <- data.frame(Year = NA, FMORT = NA, SIMID = NA)
  
  EXP_info = SAA ## same structure as SAA
  
  ## Remove pre-existing MATRIX files
  do.call(file.remove, list(list.files(path_name, full.names = TRUE)))
  write.table(SAA[-1,],paste0(path_name,"/IBM_SAA_MATRIX.txt"),
              append = TRUE,row.names = FALSE, col.names = TRUE,sep=",")
  
  write.table(dinfo[-1,],paste0(path_name,"/IBM_DEAD_MATRIX.txt"),
              append = TRUE,row.names = FALSE,col.names = TRUE,sep=",")
  
  write.table(EXP_info[-1,],paste0(path_name,"/IBM_EXP_MATRIX.txt"),append = TRUE,
              row.names = FALSE,col.names = TRUE,sep=",")
  
  # for(level in c('LOW','MED','HIGH')){ ## make all three
  # write.table(SAA[-1,],paste0(path_name,"/IBM_SAA_MATRIX_",level,".txt"),
  #             append = TRUE,row.names = FALSE, col.names = TRUE,sep=",")
  # 
  # write.table(dinfo[-1,],paste0(path_name,"/IBM_DEAD_MATRIX_",level,".txt"),
  #             append = TRUE,row.names = FALSE,col.names = TRUE,sep=",")
  # }## end of level
  # 
  # write.table(EXP_info[-1,],paste0(path_name,"/IBM_EXP_MATRIX.txt"),append = TRUE,
  #             row.names = FALSE,col.names = TRUE,sep=",")
} ## end of function

## Deprecated

## old Fgen was having issues and didn't use scale ----
# getFhist<-function(nsim,Esd,nyears,dFmin,dFmax,bb){
#   
#   ne<-nsim*10                                           # Number of simulated effort datasets
#   dEfinal<-runif(ne,dFmin,dFmax)                        # Sample the final gradient in effort
#   a<-(dEfinal-bb)/nyears                                # Derive slope to get there from intercept
#   a<-array(a,dim=c(ne,nyears))                          # Slope array
#   bb<-array(bb,dim=c(ne,nyears))                        # Intercept array
#   x<-array(rep(1:nyears,each=ne),dim=c(ne,nyears))      # Year array
#   dE<-a*x+bb                                            # Change in effort
#   E<-array(NA,dim=c(ne,nyears))                         # Define total effort array
#   E[,1]<-dE[,1]
#   for(y in 2:nyears){
#     E[,y]<-apply(dE[,1:y],1,sum)
#   }
#   E<-E/array(apply(E,1,mean),dim=c(ne,nyears))          # Standardise Effort to average 1
#   cond<-apply(E,1,min)>0
#   pos<-(1:ne)[cond]
#   pos<-pos[1:nsim]
#   
#   E<-E[pos,]                                            # Sample only those without negative effort
#   Emu<--0.5*Esd^2
#   Eerr<-array(exp(rnorm(nyears*nsim,rep(Emu,nyears),rep(Esd,nyears))),c(nsim,nyears))
#   E*Eerr
# # }
# df <- data.frame(
#   'SIM' = rep(1:nsim, each = F_start_yr+1),
#   'YEAR' = rep(seq(1, F_start_yr+1, 1), nsim),
#   'FMORT' = if (level == 'HIGH') {
#     as.vector(getFhist(
#       nsim = nsim,
#       Esd = 0,
#       nyears = F_start_yr+1, 
#       dFmin = -0.1,
#       dFmax = 0.1,
#       bb = 0.15
#     ))
#   }
#   else if (level == 'MED') {
#     as.vector(getFhist(
#       nsim = nsim,
#       Esd = 0,
#       nyears = F_start_yr+1,
#       dFmin = -0.1,
#       dFmax = 0.1,
#       bb = 0.25
#     ))
#   } else {
#     as.vector(getFhist(
#       nsim = nsim,
#       Esd = 0,
#       nyears =  F_start_yr+1,
#       dFmin = -0.15,
#       dFmax = 0.15,
#       bb = 0.35
#     ))
#   }
# )
# 
# select.effort = function(level, nsim, location){
#   ## pick a random vector (given by column SIM)
#   df = read.csv(paste0(location,'FSIM_',level,'.csv'))
#   effort.vect = floor(runif(1, 0, nsim))
#   Inst_F_vector = c(rep(0, M_only_yr), df[df$SIM == effort.vect,'FMORT'])
#   # print(summary(df$FMORT))
#   # return(Inst_F_vector)
# }