## routine to initialize the population
simu_year = 100
a1 <- 1
a2 <- 15 ## max age
Inst_M <- 0.25
growth_lognormal_SD <- 0.025
L1 <- 65
L2 <- 150
VBGF_K <- 0.3
Linf <-L1+(L2-L1)/(1-exp(-VBGF_K*(a2-a1)))
lw_a = 1.3528e-006
lw_b = 3.4297
r = -0.1034
L50 = 143.68

## some functions
# source("./functions/IBM_system_presets.R")
R0_super = 36

## just to initialize the population

init_comp <- NULL
init_comp[1] <- 1
for(g in 2:(a2+1)){
  init_comp[g]=init_comp[g-1]*exp(-Inst_M) ## instantaneous survivorship at age
}
init_comp[a2] = init_comp[a2+1]/(1-exp(-Inst_M)) #plus group
init_comp[1] <- R0_super

init_comp <- init_comp[1:a2]
# Length-at-age
init_length = NULL
error = exp(rnorm(1,0,growth_lognormal_SD)-growth_lognormal_SD^2/2)
init_length[1] = L1*error

for(g in 1:(a2)){
  init_length[g+1] = init_length[g]+growth_incre(init_length[g])*error

}
init_length=init_length[1:a2]

# Weight-at-age
init_weight = lw(init_length)
# Maturity-at-age
init_mat = prob_mature(init_length)

# Calculate SSB per recruit for initial pop
## bernoulli trial
r1 <- runif(length(init_mat), 0, 1)
# mat01 <- r1 >= init_mat ## index where r1 is first less than probability of survival (incorporates change)
# mat01[mat01 == TRUE] <- 0; mat01[mat01==FALSE] <- 1
# SSB_per_recruit = sum(0.5*init_weight*init_mat*init_comp)
# initial SSB
# SSB0 = SSB_per_recruit*R0_super; cat(SSB0,"\n") ## MAKE THIS EQUAL 2000
SSB0 = sum(0.5*init_weight*init_mat*init_comp)
# Lognormal recruitment deviation module
recruit_dev = rnorm(simu_year,0,sigma_R)
recruit_dev_adj = (recruit_dev-0.5*sigma_R^2)
RECRUITs = round(R0_super*exp(recruit_dev_adj)) ## generate estimated recruit vector of length simu_year

SSBy = NULL
SSBy[1] = SSB0

uni_ind = 0
## Simple gen F vector using original setup
Inst_F_vector <- rep(0.000,simu_year+1) ## to make them longer.

for(ry in 1:simu_year){ ## iterate sim years
  if(RECRUITs[ry] < 1){ ## fail if less than one recruit
    # RECRUITs[ry]=1 ## coerce to 1
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
  # if (ry > M_only_yr & ry <= (simu_year - 1)) { ## only perform SRR for years after M only and before last year
  
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
  r1 <- runif(length(xIY), 0, 1)
  mat01 <- r1 >= xIY
  mat01[mat01 == TRUE] <- 1; mat01[mat01==FALSE] <- 0
  # tempE = 1 * tempW ^ 0 ## ! unsure..this is supposed to be eggs
  # SSBy[ry] = sum(tempW * tempMat * tempE) / 2 ## SSB for this year (half step) is the sum of all biomasses & prob maturity
  SSBy[ry] = sum(tempW * mat01 ) / 2 ## SSB for this year is the sum of all biomasses & prob maturity; the survivorship has already been acct for
  RECRUITs[ry + 1] = round(BH_SR(SSBy[ry]) * exp(recruit_dev_adj[ry])) ## next year's recruits are the bev-holt of the SSB with error
  # } ## end of recruitment generation for that year
} ## end of simu_year (ry)
dat1 <- dat2 <- datN <- list()

# for(b in 1:nboot.temp){
  dat0 <- read.table(paste0(out_file,"/boot_",b,"/IBM_SAA_MATRIX.txt"),sep=",",header=T)
  
  dat0 %>% 
    mutate(Length_cm = fish_size) %>% 
    select(Year, Age, AgeE, Length_cm, FMORT, REG) %>% 
    write.csv(., paste0(out_file,"/boot_",b,"/",b,"_sim_Comp.csv"),row.names = F)
  
  
  dat1[[b]] <- dat0 %>%
    group_by(Year, Age) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n / sum(n)) %>% 
    filter(Year %in% c(50:99))
  
  dat1[[b]]$Flev <- 'No Fishing'
  
  
  ## length comps, use rounded vals for binning (only for illustration)
  dat2[[b]] <- dat0 %>%
    mutate(fish_size = round(fish_size)) %>%
    group_by(Year, fish_size) %>% 
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    filter(Year %in% c(50:99)) #%>%
  
  dat2[[b]]$Flev <- 'No Fishing'
  
  
  datN[[b]] <- dat1[[b]] %>% group_by(Year) %>% summarise(totN = sum(n))
  sptl <- scenarios[l,"SPATIAL"]
  scenname <- paste0(getwd(),"/IBM_output/",scenarios[l,"DESC"])
  regID <- ifelse(!is.na(scenarios[l,"REG"]),paste(scenarios[l,"REG"]),"")
  out_file <-  paste0(scenname,"/",regID)
  
  # if(fLevs[l,'CAT'] %in% c('l = 6NONE','SPATIAL')){ ## if spatial only, aggregate all and split regionally
  # nboot.temp <- ifelse(scenarios[l,"DESC"] == 'NoBreaks',nboot*2,nboot)
  
  # for(b in 1:nboot.temp){
    cat(basename(scenname)," boot ",b, " MakeLat ", "\n")
    dat <- list.files(scenname, full.names = T,recursive = T)[grep(paste0("/",b,'_sim_Comp'),list.files(scenname, full.names = T, recursive = T))] %>%
      lapply(read.csv,sep=",",header=T) %>%
      reduce(bind_rows)
    makeLat(dat) %>%
      write.csv(.,  paste0(getwd(),"/IBM_output/datasets/",basename(scenname),"_",b,".csv"),
                row.names = F)
  # }
# } ## end boots
## IGNORE was to make the plot
# require(dplyr);require(ggplot2)
read.csv("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/IBM_output/datasets/NoBreaks_44.csv") %>%
  filter(Year < 51) %>%
  group_by(Year) %>%
  dplyr::summarise(n = n()) %>%
  ggplot(.,aes(x = Year, y = n)) +geom_line(aes(color = "number alive"))+
  labs(y = '# Individuals') +
  geom_line(aes(y = RECRUITs[1:51],color = 'num recruits')) # +
  geom_line(aes(y = SSBy[1:51],color = 'SSB')) 
