require(TMB); require(dplyr); require(ggplot2); require(reshape)
options(scipen=999)
rm(list = ls())


## load objects made in dataprep
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")
load( paste0(getwd(),"/data/all_data.rda")) ## all_data
load( paste0(getwd(),"/data/DES.rda")) ## all_data



# load(paste0(getwd(),"/lenarray.rda")) ## length data
# load(paste0(getwd(),"/agearray.rda")) ## age data
# load( paste0(getwd(),"/selarray.rda")) ## selectivity values
# nStrata <- nrow(nmat)
# stnames <- c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s",
#              "AK 0","AK 1","AK 2","AK 3","AK 4","AK 5","AK 6","AK 7","AK 8")

## Estimation
# s <- c(1:4)[1] ## choose between uniform, min, maax, or dome selectivity


DES <-
compile("sptlvb.cpp")
dyn.load(dynlib("sptlvb"))

dat0 <- rep0 <- NULL ## later storage for predicts
for(s in c(1:ncol(DES))){ ## 4 hypotheses for spatial groupings

  nStrata <- length(unique(DES[,s]))
  data <-
    list(
      Length_cm = all_data[,"Length_cm"],
      Age = all_data[,"Age"],
      DES = as.vector(DES[,s]),
      nStrata = nStrata
    )

  parameters <-
    list(
      dummy = 0,
      log_Linf = rep(log(70), nStrata),
      log_k = rep(0, nStrata),
      t0 = rep(0, nStrata),
      sigma0 = 0.1,
      sigma1 = 1
    )
  
  # Now estimate everything
  map <- NULL
  model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)
  fit <- nlminb(model$par, model$fn, model$gr, control = list(rel.tol = 1e-12, eval.max = 100000,iter.max = 10000), lower = 0.01)
  for (k in 1:3)  fit <- nlminb(model$env$last.par.best, model$fn, model$gr) ## start at last-best call, for stability
  
  best <- model$env$last.par.best
  print(c(best[1:6],exp(best[7:length(best)])))
  print(fit$objective)
  rep <- sdreport(model)
  dat0 <- bind_cols(dat0, model$report()$ypreds %>% data.frame()) ## each 6 cols is new sim
  rep0 <- bind_rows(rep0, bind_cols(data.frame(names(rep$value)),data.frame(rep$value)))
  
}
# model$report()$trunc_fun
AAA


rep0 <- rep0 %>% mutate('model' = c(rep('uniform',nStrata*2*3),rep('corrected',nStrata*2*3) ),
                        'sex' = rep(c('M','F'),nStrata*6),
                        'strata' = rep(rep(stnames, each = 2),6))
# write.csv(rep0, paste0(getwd(),"/results/parEst_",Sys.Date(),'.csv'),row.names = F)              



## plotting ----
names(dat0) <- paste0(c(rep('uniform_',length(stnames)*2),
                        rep('corrected_',length(stnames)*2)),
                      stnames,c(rep("_M",length(stnames)),rep("_F",length(stnames))))


## fill precise 0 with NAs (these are autofilled)
dat0[dat0==0] <- NA

agedf <- agemat  %>% data.frame() %>% reshape::melt()  %>% plyr::rename(c("value" = "Age")) %>% select(-variable) 
lendf <- lenmat %>% data.frame() #%>% reshape::melt()  %>% plyr::rename(c("value" = "Length")) %>% select(-variable) 
# names(lendf) <-  paste0(rep('Length_',6),paste0(unique(len$st_f)))
names(lendf) <- paste0(rep('Length_',length(stnames)), stnames,c(rep("_M",length(stnames)),rep("_F",length(stnames))))



mdf <-
  bind_cols(dat0,  lendf) %>% 
  reshape::melt(id = c())  %>% 
  plyr::rename(c('variable' = 'ID', 'value' = 'Length')) %>%
  mutate(
    model = sub('_.*$', '', ID),
    st = sub("_[^_]+$", "",  sub("^(?:[^_]+_){1}", "\\1",ID)),
    Sex = sub(".*_ *(._?)", "\\1", ID),
    age = rep(agedf$Age, 3)
  )  %>%
  select(-ID) %>%
  filter(!is.na(age) & !is.na(Length))
## for plotting
mdf$st_f = factor(mdf$st, levels=stnames)
# save(mdf,file =  paste0(getwd(),"/uniform_predicts.rda"))


ggplot(subset(mdf, model == 'Length'), aes(x = age, y = Length, col = Sex, group = model)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = c(0.9,0.1))+
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  scale_color_brewer(palette = 'Dark2') +
  scale_alpha(guide = 'none') +
  geom_point(alpha = 0.2) +
  geom_line(data = subset(mdf, model != 'Length'),
            aes(x = age, y = Length, color = model, group = model), lwd = 1.1) +
  facet_wrap(~ st_f + Sex) +
  labs(title = "Predicted Model Fits and Raw Data", 
       y = 'Length (cm)', x= 'Age (yr)', 
       color = 'selectivity model')

# ggsave(file = paste0(getwd(),"/plots/dome_uni_fits2.png"), plot = last_plot(), height = 5, width = 7, unit = 'in', dpi = 520)
