## IBM Master REBOOT
## kapurm@uw.edu June 2019
## Code to generate fake population with spatial structure 
## from scenarios.csv
library(lattice)
require(dplyr)
require(purrr)
library(RColorBrewer)
require(tidyr)
require(ggplot2)
# devtools::install_github("mkapur/kaputils")
# install.packages( "c:/users/maia kapur/dropbox/kaputils",  repos = NULL, type = "source", dependencies = T)
library(kaputils)
rm(list = ls())

options(warn = -1)
nboot <- 100
scenarios <- read.csv('./input_data/scenarios.csv',na.strings = 'NA') ## manual file
source('./functions/IBM_wrangle_data.R')
# fdf <- data.frame(); idx <- 1
p <- proc.time() ## set timer

for(l in 1:(nrow(scenarios)-1)){
  
  ## if there's a spatial break
  if (!is.na(scenarios[l, 'SPATIAL'])){
    ## special VB params
    source(paste0(getwd(),"/functions/vb_params_", scenarios[l, 'REG'], '.R')); source("./functions/IBM_system_presets.R") ## system functions 
  } else {
    source("./functions/vb_params_R1.R"); source("./functions/IBM_system_presets.R") ## system functions
  }
  
  scenname <- paste0(getwd(),"/IBM_output/",scenarios[l,"DESC"]);   if(!exists(scenname)) {dir.create(scenname, showWarnings = F)}
  regID <- ifelse(!is.na(scenarios[l,"REG"]),paste(scenarios[l,"REG"]),"")
  out_file <- paste0(scenname,"/",regID,"/"); if(!exists(out_file)) {dir.create(out_file, showWarnings = F)}
  # nboot.temp <- ifelse(scenarios[l,"DESC"] == 'NoBreaks',nboot*2,nboot)
  
  for(boot_number in 1:nboot) { 
    cat(paste0("Generating population for ",scenarios[l,'DESC'],' boot ',boot_number, "\n"))
    ## Designate WRITE location
    # path_name <- paste0(out_file,"/boot_",boot_number) 
    # dir.create(path_name, showWarnings = F) ## make individual sim file
    # maketables(path_name) ## overwrite pre-existing tables, in system_presets
    
    ## Generate Equilibrium Population ----
    Neqn <- rep(0,a2)
    S <- exp(-Inst_M)
    Neqn[1] <- R0_super ## AEP set to 100
    for (Iage in 2:a2) Neqn[Iage] <- Neqn[Iage-1]*S ## survivorship
    Neqn[a2] <- Neqn[a2]/(1-S) ## plus group
    
    # Length-at-age 
    init_length = NULL
    error = exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2)
    init_length[1] = L1*error
    
    for(g in 1:(a2)){
      init_length[g+1] = init_length[g]+growth_incre(init_length[g])*error
    }
    init_length <- init_length[1:a2]
    # Weight-at-age
    init_weight = lw(init_length)
    # Maturity-at-age
    init_mat = round(prob_mature(init_length) ,1)
    SSB0 <- sum(Neqn*init_mat*init_weight) ## What we'd expect from a single population
    
    NIBM <- data.frame("Age" = NA, "Length_cm" = NA, "Mature" = NA, "Weight" = NA,
                       "Year" = NA, 'UniqueID' = NA)
    ## Populate initial individuals @ age. this number will not change, it is fixed, scaled by R0
    Ipnt <- 0
    for (Iage in 1:a2){ ## loop possible ages
      Nage <- round(Neqn[Iage]) # Integer numbers
      for (Inum in 1:Nage){ ## loop lifespan
        Ipnt <- Ipnt + 1 ## counting up individuals, new row for each
        NIBM[Ipnt,1] <- Iage-1 ## age
        # Length (initial) 
        NIBM[Ipnt,2] <- init_length[Iage]+growth_incre(init_length[Iage])*exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2) ## recalc error each time
        # Is this animal mature off the bat
        NIBM[Ipnt,3] <- runif(1,0,1) < prob_mature(NIBM[Ipnt,2])  ## logistic ogive, coinflip here gives x                     
        # Weight
        NIBM[Ipnt,4] <- lw( NIBM[Ipnt,2])
        NIBM[Ipnt,5] <- 0
        NIBM[Ipnt,6] <- Ipnt

      }  
    }  
    Ipnt0 <- Ipnt
    
    # Compute the Ps -- this is the conditional probability step, based on raw...expectation?
    P <- rep(1,a2); P[1] <- 0
    for (Iage in 2:a2) 
      if (init_mat[Iage-1] < 1)  
        P[Iage] <- (init_mat[Iage]-init_mat[Iage-1])/(1-init_mat[Iage-1])
    
    
    # Lognormal recruitment deviation module
    recruit_dev = rnorm(simu_year,0,sigma_R)
    recruit_dev_adj = (recruit_dev-0.5*sigma_R^2)
    # RECy = round(R0_super*exp(recruit_dev_adj)) ## generate estimated recruit vector of length simu_year
    RECy <- NULL
    # RECy[1] <- R0_super
    SSBy = NULL
    SSBy[1] = SSB0
    # with(NIBM, plot(Length_cm, Mature))
    # Now loop forward [this is where mine comes in]
    for(ry in 1:simu_year){ ## iterate sim years
      # Alive <- !is.na(NIBM[,'Age']) & NIBM[,'Age'] > -1
      # if(sum(is.na(NIBM[Alive,'Mature']))>0) stop(min(which(is.na(NIBM[Alive,'Mature'])))," is NA mature")
      
      # cat("Alive",sum(Alive),"\n")
      # Compute SSB
      # Compute recruits here, this is how many we need to add later
      if(ry == 1){
        SSBy[ry] <- SSB0 ## SSB from last yr informs recruits
        RECy[ry] <- round(BH_SR(SSB0) * exp(recruit_dev_adj[ry])) 
      }else{
        SSBy[ry] <- sum(NIBM$Mature[NIBM$Year == (ry-1)] * NIBM$Weight[NIBM$Year == (ry-1)]) ## Mature of living x weight of living [survived]
        RECy[ry] <- round(BH_SR(SSBy[ry]) * exp(recruit_dev_adj[ry])) 
      }
      
      if(SSBy == 0)stop(ry," SSBY = 0 \n")
      if(RECy[ry] < 1){ ## fail if less than one recruit
        RECy[ry]=1 ## coerce to 1
        print(paste0("!!! Recruitment Fail !!! year ",ry))
      }
      Ilast <- ifelse(ry == 1,1,min(which(NIBM$Year == (ry-1))))  ## where do individuals for prev year begin?
      Ipnt <- ifelse(ry == 1,Ipnt0, nrow(subset(NIBM, Year == (ry-1)))) ## number of individuals
      idx <- nrow(NIBM) ## last filled in row
      for (II in Ilast:(Ilast+Ipnt-1)){ ## loop individuals for this year [given by Ipnt]
        # if (Alive[II]){
        MyAge <- NIBM[II,1] ## check age
        # cat(NIBM[II,'UniqueID'],"\n") ## should not duplicate
        NIBM[idx+1,'UniqueID'] <- NIBM[II,'UniqueID'] ## retain ID
        NIBM[idx+1,'Year'] <- ry ## update year
        
        # Died this year?
        NIBM[idx+1,'Age'] <- ifelse(runif(1,0,1) > S, -1, ifelse(NIBM[II,'Age']+1 <= 15,NIBM[II,'Age']+1,15 )) ## Update dead, if not increase age by 1 year, coerce to 15 if older
        # Now grow the animal and update its weight
        NIBM[idx + 1, 'Length_cm'] <- NIBM[II, 'Length_cm'] +growth_incre(NIBM[II, 'Length_cm'])*exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2)
        NIBM[idx + 1, 'Weight'] <-  lw( NIBM[II,'Length_cm'])
        # Mature this year?
        if (is.na(NIBM[II,'Mature'])){NIBM[II,'Mature'] <- 0} ## replace NAs
        if (NIBM[II,'Mature'] == 1){NIBM[idx + 1,'Mature'] <- 1} ## cannot un-mature
        if(MyAge < 15 & NIBM[II,'Mature'] == 0){
          NIBM[idx + 1,'Mature'] <- ifelse(runif(1,0,1) < P[MyAge+1],1,0)
        } else if(MyAge >= 15 & NIBM[II,'Mature'] == 0){
          NIBM[idx + 1,'Mature'] <- ifelse(runif(1,0,1) < P[MyAge],1,0)
        }
        idx <- idx + 1 ## bump down a row
        # Add recruits starting from slots Ipnt + 1
        
        # Ipnt <- Ipnt + 1 ## counting up individuals, new row for each
        # } ## end if living
      } ## end individuals
      for(ir in 1:RECy[ry]){ ## now add individuals
        NIBM[idx + 1,'UniqueID'] <- max(NIBM$UniqueID) + 1 ## new id; DO FIRST OTHERWISE FILLS NA
        NIBM[idx + 1,'Year'] <- ry
        NIBM[idx + 1,'Age'] <- 0
        NIBM[idx + 1, 'Length_cm'] <- L1*exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2) ## recalc error each time
        NIBM[idx + 1, 'Weight'] <-  lw(NIBM[II,'Length_cm'])
        NIBM[idx + 1,'Mature'] <- ifelse(runif(1,0,1)< P[1],1,0)
        idx <- idx + 1
      } ## end adding recruits
      
      NIBM <- NIBM %>% filter(Age != -1) # To speed up remove all Age -1 animals at each time-step. idx is invalid now.
    } ## end simu_year
    ## save NIBM
    NIBM <- NIBM %>% filter(Age != -1 &  Year > 24) # To speed up remove all Age -1 animals at each time-step. idx is invalid now.
    NIBM$Year <- NIBM$Year - 25 ## remove burn in
    
    write.csv(NIBM %>% mutate('boot' = boot_number, 
                              'REG' = as.factor(scenarios[l, "REG"])),
              paste0(out_file,"/",boot_number,"_sim_Comp.csv"),row.names = F)
    
    
    
  } ## end of boots
  } ## end make NIBM

#   out_file <- paste0("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/IBM_output/tempvar_R1R2")
#   # out_file <- paste0("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/IBM_output/NoBreaks")
 
# ## check length vs time -- should be no swoops
#   list.files(out_file, pattern = "*.csv",recursive = F,full.names = T) %>%
#     lapply(., read.csv) %>%
#     bind_rows() %>%
#     filter(Age == 6) %>%
#     group_by(boot) %>%
#     # dplyr::summarise(n = n()) %>%
#     ggplot(.,aes(x = Year, y = Length_cm, group = boot, color = factor(boot))) +
#     geom_point() +
#     labs(y = 'Length_cm', color = 'replicate')
#   ## check n individuals vs time -- should be no swoops
#   list.files(out_file, pattern = '*.csv',recursive = F,full.names = T) %>%
#     lapply(., read.csv) %>%
#     bind_rows() %>%
#     # filter(Year != 0) %>%
#     group_by(Year, boot) %>%
#     dplyr::summarise(n = n()) %>%
#     ggplot(.,aes(x = Year, y = n, group = boot, color = factor(boot))) +
#     geom_line() +
#     labs(y = '# Individuals', color = 'replicate')


## now make latitude
# sptl <- scenarios[l,"SPATIAL"]
# scenname <- paste0(getwd(),"/IBM_output/",scenarios[l,"DESC"])
# regID <- ifelse(!is.na(scenarios[l,"REG"]),paste(scenarios[l,"REG"]),"")
# out_file <-  paste0(scenname,"/",regID)
# 
# dat <- list.files(scenname, full.names = T,recursive = T)[grep(paste0("/",boot_number,'_sim_Comp'),list.files(scenname, full.names = T, recursive = T))] %>%
#   lapply(read.csv,sep=",",header=T) %>%
#   reduce(bind_rows)
# makeLat(dat) %>%
#   write.csv(.,  paste0(getwd(),"/IBM_output/datasets/",basename(scenname),"_",boot_number,".csv"),
#             row.names = F)

## Add in Spatial - needs to happen once R1 and R2 are built for this scenario ----
source("./functions/makeLat.R")

for(l in 1:(length(unique(scenarios$DESC)))){
  sptl <- subset(scenarios, DESC == unique(scenarios$DESC)[l])$SPATIAL[1] 
  scenname <- paste0(getwd(),"/IBM_output/",unique(scenarios$DESC)[l])
  # regID <- ifelse(!is.na(subset(scenarios, DESC == unique(scenarios$DESC)[l])$REG[1] ),
  #                 paste(subset(scenarios, DESC == unique(scenarios$DESC)[l])$REG[1]),"")
  # # out_file <-  paste0(scenname,"/",regID)
  
  # if(fLevs[l,'CAT'] %in% c('l = 6NONE','SPATIAL')){ ## if spatial only, aggregate all and split regionally
  # nboot.temp <- ifelse(scenarios[l,"DESC"] == 'NoBreaks',nboot*2,nboot)
  
  for(b in 1:nboot){
    cat(paste0("Making Spatial Breaks for ",unique(scenarios$DESC)[l]," boot ",b,"\n"))
    
    # cat(basename(scenname)," boot ",b, " MakeLat ", "\n")
    dat <- list.files(scenname, full.names = T,recursive = T)[grep(paste0("/",b,'_sim_Comp'),
                                                                   list.files(scenname, full.names = T, recursive = T))] %>%
      lapply(read.csv,sep=",",header=T) %>%
      reduce(bind_rows)
    makeLat(dat) %>%
      write.csv(.,  paste0(getwd(),"/IBM_output/datasets/",basename(scenname),"_",b,".csv"),
                row.names = F)
    
    
  }
} ## end testrows


# tt <- read_csv("IBM_output/datasets/F0L1S_25_1.csv")
# unique(tt$REG)
# with(tt, plot(Length_cm~ Latitude_dd))

## Post Hoc -- Creation of temp var via stitching
for(b in 1:nboot){
  cat(paste0("Making temporal Breaks boot ", b,"\n"))
  
  p1 <- read.csv(paste0(getwd(),"/IBM_output/datasets/NoBreaks_",b,".csv")) %>% filter(Year < 50)
  # for(n in c('F0L1S_25_',"F0L1S_R3_")){
  # for(n in c('F0L1S_25_',"F0L1S_R3_")){

  p2 <- read.csv(paste0(getwd(),"/IBM_output/datasets/F0L1S_25_",b,".csv")) %>% filter(Year >= 50 & REG != 'R1')
  tempreg <- unique(p2$REG)
  tempdf <- rbind(p1,p2)
  for(i in 1:nrow(tempdf)){
    tempdf$Latitude_dd[i] <- runif(1, 0.0, 50.0) ## uniform range
    tempdf$Longitude_dd[i] <- runif(1, 0.0, 50.0)## uniform range
  } ## end rows
  levels(tempdf$REG) <- c('R0','R1','R2','R3')
  tempdf$REG[tempdf$REG == 'R0'] <- "R1"
  write.csv(tempdf, paste0(getwd(),"/IBM_output/datasets/tempvar_R1",tempreg,"_",b,".csv"), row.names=F)
  # cat(n," ",b,"\n")
  # } ## end R2/R3
} ## end boots
(proc.time() - p)[1] / 60 ## minutes

