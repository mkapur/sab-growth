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
testrows <- 1:nrow(scenarios)
fdf <- data.frame(); idx <- 1
p <- proc.time() ## set timer

for(l in testrows){
  
  ## if there's a spatial break
  if (!is.na(scenarios[l, 'SPATIAL'])){
    ## special VB params
    source(paste0(getwd(),"/functions/vb_params_", scenarios[l, 'REG'], '.R')); source("./functions/IBM_system_presets.R") ## system functions 
  } else {
    source("./functions/vb_params_R1.R"); source("./functions/IBM_system_presets.R") ## system functions
  }
  
  scenname <- paste0(getwd(),"/IBM_output/",scenarios[l,"DESC"]);   if(!exists(scenname)) {dir.create(scenname, showWarnings = F)}
  regID <- ifelse(!is.na(scenarios[l,"REG"]),paste(scenarios[l,"REG"]),"")
  out_file <- paste0(scenname,"/",regID); if(!exists(out_file)) {dir.create(out_file, showWarnings = F)}
  nboot.temp <- ifelse(scenarios[l,"DESC"] == 'NoBreaks',nboot*2,nboot)
  
  for(boot_number in 1:nboot.temp) { ## full nsims, line 20
    cat(paste0(paste(scenarios[l,'DESC']," ",scenarios[l,'ID']),' boot ',boot_number, "\n"))
    ## Designate WRITE location
    path_name <- paste0(out_file,"/boot_",boot_number) 
    dir.create(path_name, showWarnings = F) ## make individual sim file
    maketables(path_name) ## overwrite pre-existing tables, in system_presets
    
    ## Generate Equilibrium Population
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
    # Weight-at-age
    init_weight = lw(init_length)
    # Maturity-at-age
    init_mat = prob_mature(init_length)    
    SSB0 <- sum(Neqn*init_mat*init_weight) ## What we'd expect from a single population

    # Nmature (for checking)
    # print(sum(Neqn*Matt))
    # NCnt <-20*sum(Neqn) ## never used again
    # NIBM <- matrix(NA,nrow=20*sum(Neqn),ncol=4) ## empty matrix
    NIBM <- data.frame("Age" = NA, "Length_cm" = NA, "Mature" = NA, "Weight" = NA, "Year" = NA, 'UniqueID' = NA)
    ## Populate initial individuals @ age
    Ipnt <- 0
    for (Iage in 1:a2){ ## loop possible ages
      Nage <- round(Neqn[Iage]) # Integer numbers
      for (Inum in 1:Nage){ ## loop lifespan
        Ipnt <- Ipnt + 1 ## counting up individuals, new row for each
        NIBM[Ipnt,1] <- Iage-1 ## age
        # Length (initial) - don't worry this will burn out quickly
        NIBM[Ipnt,2] <- init_length[Iage]+growth_incre(init_length[Iage])*exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2) ## recalc error each time
        # NIBM[Ipnt,5] <-100*(1.0-exp(-0.2*(Iage-1)))*exp(rnorm(1,0,0.2)-0.2^2/2.0) ## bias corr length
        # Is this animal mature off the bat
        NIBM[Ipnt,3] <- runif(1,0,1) < prob_mature(NIBM[Ipnt,2])  ## logistic ogive, coinflip here gives x                     
        # NIBM[Ipnt,3] <- runif(1,0,1) < Matt[Iage]   ## logistic ogive                       
        # Weight
        NIBM[Ipnt,4] <- lw( NIBM[Ipnt,2])
        NIBM[Ipnt,5] <- 'EQ'
        NIBM[Ipnt,6] <- Ipnt
      }  
    }  
    # head(subset(NIBM, Age == 10))
    # How many animals
    Ipnt0 <- Ipnt
    # How many mature animals
    # print(sum(NIBM[1:Ipnt,3]))
    
    # Compute the Ps -- this is the conditional probability step, based on raw...expectation?
    P <- rep(1,a2); P[1] <- 0
    for (Iage in 2:a2) 
      if (init_mat[Iage-1] < 1)  
        P[Iage] <- (init_mat[Iage]-init_mat[Iage-1])/(1-init_mat[Iage-1])
    # print(P)
    
    # Lognormal recruitment deviation module
    recruit_dev = rnorm(simu_year,0,sigma_R)
    recruit_dev_adj = (recruit_dev-0.5*sigma_R^2)
    RECy = round(R0_super*exp(recruit_dev_adj)) ## generate estimated recruit vector of length simu_year
    
    SSBy = NULL
    SSBy[1] = SSB0
    uni_ind = 0
    
    # Now loop forward [this is where mine comes in]
    for(ry in 1:3){ ## iterate sim years
      
      Alive <- !is.na(NIBM[,'Age']) & NIBM[,'Age'] > -1
      if(sum(is.na(NIBM[Alive,'Mature']))>0) stop(min(which(is.na(NIBM[Alive,'Mature'])))," is NA mature")
      
      # cat("Alive",sum(Alive),"\n")
      # Compute SSB
      SSBy[ry] <- sum(NIBM[Alive,'Mature']*NIBM[Alive,'Weight']) ## Mature of living x weight of living [survived]
      # Compute recruits here, this is how many we need to add later
      if(ry == 1){
        RECy[ry] <- round(BH_SR(SSB0) * exp(recruit_dev_adj[ry])) 
      } else{
        RECy[ry] <- round(BH_SR(SSBy[ry]) * exp(recruit_dev_adj[ry])) 
      }
      Ilast <- ifelse(ry == 1,1,min(which(NIBM$Year == (ry-1))))  ## where do individuals for prev year begin?
      Ipnt <- ifelse(ry == 1,Ipnt0, nrow(subset(NIBM, Year == (ry-1)))) ## number of individuals
      idx <- nrow(NIBM) ## last filled in row
      for (II in Ilast:(Ilast+Ipnt-1)){ ## loop individuals for this year [given by Ipnt]
        if (Alive[II]){
          MyAge <- NIBM[II,1] ## check age
          cat(MyAge," ",II,"\n")
          NIBM[idx+1,'UniqueID'] <- NIBM[II,'UniqueID'] ## retain ID
          NIBM[idx+1,'Year'] <- ry ## update year
          
          # Died this year?
          NIBM[idx+1,'Age'] <- ifelse(runif(1,0,1) > S, -1, ifelse(NIBM[II,'Age']+1 <= 15,NIBM[II,'Age']+1,15 )) ## Update if dead, coerce to 15 if older
          # Now grow the animal and update its weight
          NIBM[idx + 1, 'Length_cm'] <- NIBM[II, 'Length_cm'] +growth_incre(NIBM[II, 'Length_cm'])*exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2)
          NIBM[idx + 1, 'Weight'] <-  lw( NIBM[II,'Length_cm'])
          # Mature this year?
          if (is.na(NIBM[II,'Mature'])){NIBM[II,'Mature'] <- 0}
          if(MyAge <= 15){
            NIBM[idx + 1,'Mature'] <- ifelse(runif(1,0,1) < P[MyAge+1],1,0)
          } else{
            NIBM[idx + 1,'Mature'] <- ifelse(runif(1,0,1) < P[MyAge],1,0)
            
          }
       ifelse(MyAge <= 15 &  second coinflip
          idx <- idx + 1 ## bump down a row
          # Add recruits starting from slots Ipnt + 1
            
        # Ipnt <- Ipnt + 1 ## counting up individuals, new row for each
        } ## end if living
      } ## end individuals
      for(r in 1:RECy[ry]){ ## now add individuals
        NIBM[idx + 1,'Year'] <- ry
        NIBM[idx + 1,'Age'] <- 0
        NIBM[idx + 1,'UniqueID'] <- max(NIBM$UniqueID) + 1 ## new id
        NIBM[idx + 1, 'Length_cm'] <- L1*exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2) ## recalc error each time
        NIBM[idx + 1, 'Weight'] <-  lw(NIBM[II,'Length_cm'])
        NIBM[idx + 1,'Mature'] <- ifelse(runif(1,0,1)< P[1],1,0)
        idx <- idx + 1
      } ## end adding recruits

      NIBM <- NIBM %>% filter(Age != -1 ) # To speed up remove all Age -1 animals at each time-step. idx is invalid now.
      # NIBM[50:75,]
      cat(idx," ",ry,"\n")    } ## end simu_year

    
    
    
      if(RECy[ry] < 1){ ## fail if less than one recruit
        RECy[ry]=1 ## coerce to 1
        print(paste0("!!! Recruitment Fail !!! year ",ry))
      }
      
      # current_year_by individual
      RECRUITs_year =  rep(1*(ry-1),RECy[ry]) ## vector that lists current year once for each recruit 
      
      # a loop for running each by ind  
      for(i in 1:length(RECRUITs_year)){   ## iterate through each new recruit
        
        uni_ind <- uni_ind + 1
        accumu_yr_1 <- RECRUITs_year[i]+start_yr   
        
        ## GROWTH MODULE: 
        num_of_chance <- length(Inst_F_vector[ry:length(Inst_F_vector)])
        # num_of_chance =  length(ry:length(Inst_F_vector)) ## find length from current recruitment year to end of fishing mortality vector (includes burn in); number of years left
        
        fish_sizes = NULL
        ind_growth_error = exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2) ## set individual growth error
        fish_sizes[1]=L1*ind_growth_error ## establish first size
        
        for(g in 1:(num_of_chance-1)){          # next mid year (interval 1 year)
          # fish_sizes[g+1] <- (fish_sizes[g]+growth_incre(fish_sizes[g]))*ind_growth_error
          fish_sizes[g+1] <- fish_sizes[g]+growth_incre(fish_sizes[g])*ind_growth_error
          
        }
        
        
        # Generate a vector of inst_total mortality for this individual;
        ## inst_M is constant, select F vector from given year to end of F vector and multiply by selectivity across fish's lifetime
        inst_Z_vector = Inst_M + Inst_F_vector[ry:length(Inst_F_vector)]*selectivity(fish_sizes)  
        # subset the fish size from the full size-at-age vector (only those that survive grow)
        
        inst_Z_vector = inst_Z_vector[1:a2]
        
        # random number r1 - bernoulli trial
        r1 = runif(length(inst_Z_vector), 0, 1) ## bernoulli trial for each ind/yr/age
        Zdeath = which.max(r1 <= (1 - exp(-inst_Z_vector))) ## index where r1 is first less than probability of survival -- first death
        fish_sizes = fish_sizes[1:Zdeath] ## trim fish_sizes to time at death
        ## as whichever age the random number exceeds e^-m, survivorship goes to zero ()
        if(length(fish_sizes) >a2){print('length fish size > a2 line 179')}
        
        ## SIZE-AT-TIME BOOKKEEPING 
        true_age = 1:Zdeath-1 ## actual lifespan before zdeath
        AgeE = rnorm(length(true_age),mean_age_true_age[true_age+1],mean_SD_age_true_age[true_age+1]) 
        AgeE[AgeE<0]=0
        
        SAA <- data.frame(
          cohort = rep(ry, Zdeath),
          uni_ind = rep(uni_ind, Zdeath), ## changed this to rep length Zdeath so it'd match
          ind = rep(i, Zdeath),
          fish_size = fish_sizes[1:Zdeath],
          Year = accumu_yr_1 + 1:Zdeath - 1,
          Age = true_age,
          AgeE = AgeE,
          FMORT = c(Inst_F_vector[ry:(ry + Zdeath - 1)]),
          M = rep(Inst_M, Zdeath),  ## changed this to rep length Zdeath so it'd match
          Sel = c(selectivity(fish_sizes)),
          # Sel = c(selectivity(true_age)),
          
          # FSIM = paste(scenarios[l,"FMORT_1"],scenarios[l,"FMORT_2"],scenarios[l,"FMORT_3"]),
          REG = scenarios[l,"REG"]
        )
        
        ## MORTALITY MODULE:
        # seperate inst_total mortality into F and M
        # M_dead_por = (1 - exp(-Inst_M)) / (1 - exp(-inst_Z_vector[Zdeath])) ## nat survivorship over overall survivorship. Will be 1 for years with only nat mort in effect.
        # r2 = runif(1, 0, 1) ## bernoulli trial for survivorship
        true_age = Zdeath - 1
        AgeE = rnorm(1, mean_age_true_age[true_age + 1], mean_SD_age_true_age[true_age +
                                                                                1])
        AgeE[AgeE < 0] = 0
        
        
        ## write data to file - one for each F level
        write.table(SAA, paste(path_name,"/IBM_SAA_MATRIX.txt",sep=""),
                    append = TRUE,
                    row.names = FALSE,
                    col.names = FALSE,
                    sep=",")
        # write.table(dinfo,paste(path_name,"/IBM_DEAD_MATRIX.txt",sep=""),append = TRUE,row.names = FALSE,col.names = FALSE,sep=",")
        
        
      } ## end of individual recruit
      
      ## STOCK-RECRUITMENT MODULE
      # Update number of recruits by B-H stock recruitment relationship
      # Recruitment derived from B-H S-R start from the second year of fishing mortality
      xIY <-  NULL; idx = 1
      if (ry > M_only_yr & ry <= (simu_year - 1)) { ## only perform SRR for years after M only and before last year
        
        tempSAA = read.table(paste0(path_name, "/IBM_SAA_MATRIX.txt"),header = T, sep = ',') ## read current size-at-age-table
        tempN_prev = subset(tempSAA, Year == (ry + start_yr - 1)) ## extract PREVIOUS year
        tempN_curr = subset(tempSAA, Year == (ry + start_yr )) ## extract CURRENT year
        
        tempW = lw(tempN_curr$fish_size) ## calc weights
        for(i in unique(tempN_curr$uni_ind)){
          tempN_prev_i <- subset(tempN_prev, uni_ind == i); if(nrow(tempN_prev_i) == 0) next(); ## skip if it died last year
          tempN_curr_i <- subset(tempN_curr, uni_ind == i)
          
          tempMat_prev = prob_mature(tempN_prev_i$fish_size) 
          tempMat_curr = prob_mature(tempN_curr_i$fish_size) ## calc prob of maturity at size, year for individual
          # r1 <- runif(1,0,1)
          # mat01_prev <- ifelse(tempMat_prev > r1,1,0)
          # mat01_curr <- ifelse(tempMat_curr > r1,1,0)
          xIY[idx] <-       (tempMat_curr - tempMat_prev)/(1-tempMat_prev); idx = idx+1 ## this is the prob mat for each individual in this year
        }
        
        
        ## bernoulli trial for maturity
        # r1 <- runif(length(tempMat), 0, 1)
        # mat01 <- (r1 <= tempMat) ## index where r1 is first less than probability of survival (incorporates change)
        # mat01[mat01 == TRUE] <- 1; mat01[mat01==FALSE] <- 0
        # tempE = 1 * tempW ^ 0 ## ! unsure..this is supposed to be eggs
        # SSBy[ry] = sum(tempW * tempMat * tempE) / 2 ## SSB for this year (half step) is the sum of all biomasses & prob maturity
        SSBy[ry] = sum(tempW * xIY ) / 2 ## SSB for this year is the sum of all biomasses & prob maturity
        RECy[ry + 1] = round(BH_SR(SSBy[ry]) * exp(recruit_dev_adj[ry])) ## next year's recruits are the bev-holt of the SSB with error
      } ## end of recruitment generation for that year
    } ## end of simu_year (ry)
   
    # 
    # init_comp <- NULL
    # init_comp[1] <- 1
    # for(g in 2:(a2+1)){
    #   init_comp[g]=init_comp[g-1]*exp(-Inst_M) ## instantaneous survivorship at age
    # }
    # init_comp[a2+1] = init_comp[a2+1]/(1-exp(-Inst_M)) #plus group
    # 

    # init_length=init_length[1:a2]
    

    # Eggs-body mass relationship
    # init_eggs = 1*init_weight^0
    # Calculate SSB per recruit
    ## bernoulli trial
    # r1 <- runif(length(init_mat), 0, 1)
    # mat01 <- (r1 <= init_mat) ## index where r1 is first less than probability of survival (incorporates change)
    # mat01[mat01 == TRUE] <- 1; mat01[mat01==FALSE] <- 0
    # # SSB_per_recruit = sum(0.5*init_comp*init_weight*init_mat*init_eggs)
    # SSB_per_recruit = sum(0.5*init_weight*mat01) ## sum weight * maturity * sex ratio
    # # initial SSB
    # SSB0 = SSB_per_recruit*R0_super; cat(SSB0,"\n")
    
 
    ## Simple gen F vector using original etup
    Inst_F_vector <- rep(0.000,simu_year+1) ## to make them longer.
    
    ## save simulated F
    # fdf <- rbind(fdf, data.frame(cbind(Inst_F_vector,
    #                                    1:101, 
    #                                    rep(idx, 101), 
    #                                    rep(boot_number,101))))
    
   
  } ## end of boots
  idx <- idx + 1
  genDat <- build_simComp(out_file,l=l)
  rm(genDat)
  
}