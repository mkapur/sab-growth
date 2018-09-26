require(TMB)
options(scipen=999)
## load data
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")
load( paste0(getwd(),"/filtered_SAB_WC.rda")); load(paste0(getwd(),"/lenmat_WC.rda")); load(paste0(getwd(),"/agemat_WC.rda")) ## made in dataprep
nStrata <- length(unique(len$st_f))

## Estimation
s <- c(1:4)[1] ## choose between uniform, min, max, or dome selectivity
minSel <- c(42, 37, 39, 27, 23, 21) ## cutoff values for selectivity
maxSel <- c(74, 70, 77, 71, 71, 63)


## a function to get 500 samples within selectivity bounds
data.sub <- function(len, s){
  ## reshape raw length and age data
  agemat0 <-
    len %>% select(Age, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID'),'Age') %>% 
    reshape::cast(., ID + value ~ st, fill = NA, drop= T) %>% select(-ID, -value) 
  lenmat0 <- 
    len %>% select(Length_cm, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID'),'Length_cm') %>% 
    reshape::cast(., ID + value ~ st, fill = NA, drop= T) %>% select(-ID, -value) 
  
  agemat <- lenmat <-  matrix(NA, nrow = 500, ncol = length(unique(len$st)))
  
  for(i in 1:nStrata){
    temp0 <- matrix(na.omit(lenmat0[,i]),ncol = 1) ##if there is a value in one, should be in the other as well 
    if(s == 1) unif <- temp <- temp0 ## used in all sims
    if(s == 2) { temp <- temp0[temp0 >= minSel[i]]}
    # if(s == 3) { temp <- temp[temp <= maxSel[i]]}  ## sample size issues
    # if(s == 4) { temp <- temp[temp >= minSel[i] & temp <= maxSel[i]]} ## sample size issues
    idx <- sample(1:length(temp),500) ## pick 500 rows and store index
    lenmat[,i] <- temp[idx] %>% as.matrix()
    agemat[,i] <- data.frame(na.omit(agemat0[,i]))[idx,] %>% as.matrix() ## select same rows
    
    lenmat[,i] <- unif[idx] %>% as.matrix()
    agemat[,i] <- data.frame(na.omit(agemat0[,i]))[idx,] %>% as.matrix() ## select same rows
    
    rm(temp0)
  }
  return(list(lenmat,agemat))
}



for(s in 1:2){ ## ultimately loop over 4 selectivities
  matrices <- data.sub(len, s)
  
  data <-
    list(
      Length_cm = matrices[[1]], 
      Age = matrices[[2]],
      st = len[, 'st_f'],
      nStrata = nStrata,
      minSel = minSel,
      maxSel = maxSel,
      selType = s
    )
  parameters <-
    list(
      dummy = 0,
      log_Linf = rep(log(70), nStrata),
      log_k = rep(0, nStrata),
      t0 = rep(0, nStrata),
      log_Sigma = 0.1
    )
  
  compile("sptlvb.cpp")
  dyn.load(dynlib("sptlvb"))
  # Now estimate everything
  map <- NULL
  model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)
  
  fit <- nlminb(
    model$par,
    model$fn,
    model$gr,
    control = list(
      rel.tol = 1e-12,
      eval.max = 100000,
      iter.max = 1000
    )
  )
  best <- model$env$last.par.best
  print(c(best[1:6],exp(best[7:length(best)])))
  rep <- sdreport(model)
}
## plotting
require(ggplot2);require(dplyr);require(reshape)
dat0 <- model$report()$ypreds %>% data.frame() 
names(dat0) <- paste0(unique(len$st_f))
dat0 <- dat0 %>% reshape::melt() %>% plyr::rename(c('value' = 'Predicted', 'variable' = 'st'))
dat0$st_f <- factor(dat0$st, levels=c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s")) 

agedf <- agemat %>% data.frame() %>% reshape::melt() %>% plyr::rename(c('value' = 'Age'))
lendf <- lenmat %>% data.frame() %>% reshape::melt() %>% plyr::rename(c('value' = 'Length'))

mdf <- bind_cols(dat0,agedf,lendf) %>% select(-variable, -variable1, -st)
# mdf$st_f = factor(mdf$variable, levels=c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s"))
# save(mdf,file =  paste0(getwd(),"/uniform_predicts.rda"))

ggplot(mdf, aes(x = Age, y = Length)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = c(0.9,0.15))+
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  geom_point(alpha = 0.2) +
  geom_line(aes(y = Predicted), color = 'blue', lwd = 1.1)+
  facet_wrap(~ st_f) +
  labs(title = "Predicted Model Fits with Uniform Selectivity", y = 'Length (cm)', x= 'Age (yr)', subtitle = 'Points are actual subsampled data used in parameter fitting')

# ggsave(file = paste0(getwd(),"/plots/uniform_fits_sub.png"), plot = last_plot(), height = 5, width = 7, unit = 'in', dpi = 520)


ggplot(len, aes(x = Age, y = Length_cm, color = Sex)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = c(0.9,0.15))+
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  scale_color_brewer(palette = 'Dark2') +
  geom_point(alpha = 0.2) +
  geom_line(data = mdf, aes(x = Age, y = Predicted), col = 'blue', lwd = 1.1) +
  facet_wrap(~ st_f) +
  labs(title = "Predicted Model Fits with Uniform Selectivity", y = 'Length (cm)', x= 'Age (yr)', subtitle = 'Parameters estimated using subset of each strata')

# ggsave(file = paste0(getwd(),"/plots/uniform_fits.png"), plot = last_plot(), height = 5, width = 7, unit = 'in', dpi = 520)
