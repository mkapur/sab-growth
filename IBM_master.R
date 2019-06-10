## IBM Master
## kapurm@uw.edu Winter 2019
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
    
    init_comp <- NULL
    init_comp[1] <- 1
    for(g in 2:(a2+1)){
      init_comp[g]=init_comp[g-1]*exp(-Inst_M) ## instantaneous survivorship at age
    }
    init_comp[a2+1] = init_comp[a2+1]/(1-exp(-Inst_M)) #plus group
    
    # Length-at-age 
    init_length = NULL
    error = exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2)
    init_length[1] = L1*error
    
    for(g in 1:(a2)){
      # init_length[g+1] = (init_length[g]+growth_incre(init_length[g]))*error
      init_length[g+1] = init_length[g]+growth_incre(init_length[g])*error
      
    }
    # init_length=init_length[1:a2]
    
    # Weight-at-age
    init_weight = lw(init_length)
    # Maturity-at-age
    init_mat = prob_mature(init_length)
    # Eggs-body mass relationship
    init_eggs = 1*init_weight^0
    # Calculate SSB per recruit
    ## bernoulli trial
    r1 <- runif(length(init_mat), 0, 1)
    mat01 <- (r1 <= init_mat) ## index where r1 is first less than probability of survival (incorporates change)
    mat01[mat01 == TRUE] <- 1; mat01[mat01==FALSE] <- 0
    # SSB_per_recruit = sum(0.5*init_comp*init_weight*init_mat*init_eggs)
    
    SSB_per_recruit = sum(0.5*init_comp*init_weight*mat01*init_eggs)
    # initial SSB
    SSB0 = SSB_per_recruit*R0_super; cat(SSB0,"\n")
    
    # Lognormal recruitment deviation module
    recruit_dev = rnorm(simu_year,0,sigma_R)
    recruit_dev_adj = (recruit_dev-0.5*sigma_R^2)
    RECRUITs = round(R0_super*exp(recruit_dev_adj)) ## generate estimated recruit vector of length simu_year
    
    SSBy = NULL
    SSBy[1] = SSB0
    
    uni_ind = 0
    ## Simple gen F vector using original etup
    Inst_F_vector <- rep(0.000,simu_year+1) ## to make them longer.
    
    ## save simulated F
    fdf <- rbind(fdf, data.frame(cbind(Inst_F_vector,
                                       1:101, 
                                       rep(idx, 101), 
                                       rep(boot_number,101))))
    
    for(ry in 1:simu_year){ ## iterate sim years
      if(RECRUITs[ry] < 1){ ## fail if less than one recruit
        RECRUITs[ry]=1 ## coerce to 1
        print(paste0("!!! Recruitment Fail !!! year ",ry))
      }
      
      # current_year_by individual
      RECRUITs_year =  rep(1*(ry-1),RECRUITs[ry]) ## vector that lists current year once for each recruit 
      
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
        r1 = runif(length(inst_Z_vector), 0, 1)
        Zdeath = which.max(r1 <= (1 - exp(-inst_Z_vector))) ## index where r1 is first less than probability of survival (incorporates change)
        fish_sizes = fish_sizes[1:Zdeath] ## trim fish_sizes to time at death
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
        M_dead_por = (1 - exp(-Inst_M)) / (1 - exp(-inst_Z_vector[Zdeath])) ## nat survivorship over overall survivorship. Will be 1 for years with only nat mort in effect.
        r2 = runif(1, 0, 1) ## bernoulli trial
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
      
      if (ry > M_only_yr & ry <= (simu_year - 1)) { ## only perform SRR for years after M only and before last year
        tempSAA = read.table(paste0(path_name, "/IBM_SAA_MATRIX.txt"),header = T, sep = ',') ## read current size-at-age-table
        tempN = subset(tempSAA, Year == (ry + start_yr - 1)) ## extract PREVIOUS year
        tempW = lw(tempN$fish_size) ## calc weights
        tempMat = prob_mature(tempN$fish_size) ## calc prob of maturity at size

        ## bernoulli trial
        r1 <- runif(length(tempMat), 0, 1)
        mat01 <- (r1 <= tempMat) ## index where r1 is first less than probability of survival (incorporates change)
        mat01[mat01 == TRUE] <- 1; mat01[mat01==FALSE] <- 0
        tempE = 1 * tempW ^ 0 ## ! unsure
        # SSBy[ry] = sum(tempW * tempMat * tempE) / 2 ## SSB for this year (half step) is the sum of all biomasses & prob maturity
        SSBy[ry] = sum(tempW * mat01 * tempE) / 2 ## SSB for this year (half step) is the sum of all biomasses & prob maturity
        
        RECRUITs[ry + 1] = round(BH_SR(SSBy[ry]) * exp(recruit_dev_adj[ry])) ## next year's recruits are the bev-holt or LFSR of the SSB with error
      } ## end of recruitment generation for that year
    } ## end of simu_year (ry)
  } ## end of boots
  idx <- idx + 1
  genDat <- build_simComp(out_file,l=l)
  rm(genDat)
} ## end of fLev[l]
print(paste("Full simulation time ",((proc.time()-p)/60)[1],"min(s)"))

## Add in Spatial
for(l in testrows){
  
  sptl <- scenarios[l,"SPATIAL"]
  scenname <- paste0(getwd(),"/IBM_output/",scenarios[l,"DESC"])
  regID <- ifelse(!is.na(scenarios[l,"REG"]),paste(scenarios[l,"REG"]),"")
  out_file <-  paste0(scenname,"/",regID)
  
  # if(fLevs[l,'CAT'] %in% c('l = 6NONE','SPATIAL')){ ## if spatial only, aggregate all and split regionally
  nboot.temp <- ifelse(scenarios[l,"DESC"] == 'NoBreaks',nboot*2,nboot)
  
  for(b in 1:nboot.temp){
    cat(basename(scenname)," boot ",b, " MakeLat ", "\n")
    dat <- list.files(scenname, full.names = T,recursive = T)[grep(paste0("/",b,'_sim_Comp'),list.files(scenname, full.names = T, recursive = T))] %>%
      lapply(read.csv,sep=",",header=T) %>%
      reduce(bind_rows)
    makeLat(dat) %>%
      write.csv(.,  paste0(getwd(),"/IBM_output/datasets/",basename(scenname),"_",b,".csv"),
                row.names = F)
    
    
  }
} ## end testrows


## Post Hoc -- Creation of temp var via stitching
for(b in 1:nboot){
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

