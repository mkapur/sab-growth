require(TMB); require(dplyr); require(ggplot2); require(reshape)
options(scipen=999)
rm(list = ls())
## load data
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")
load( paste0(getwd(),"/filtered_SAB_WC.rda"));load( paste0(getwd(),"/nmat.rda")); load(paste0(getwd(),"/lenarray.rda")); load(paste0(getwd(),"/agearray.rda")) ## made in dataprep
nStrata <- length(unique(len$st_f))

## Estimation
s <- c(1:4)[1] ## choose between uniform, min, maax, or dome selectivity
minSel <- c(42, 37, 39, 27, 23, 21) ## cutoff values for selectivity
maxSel <- c(74, 70, 77, 71, 71, 63)

parameters <-
  list(
    dummy = 0,
    log_Linf = rep(log(70), nStrata),
    log_k = rep(0, nStrata),
    t0 = rep(0, nStrata),
    log_Sigma = 1
  )

compile("sptlvb.cpp")
dyn.load(dynlib("sptlvb"))
dat0 <- NULL ## later storage for preditcs
for(s in c(1)){ ## ultimately loop over 4 selectivities
  data <-
    list(
      Length_cm = lenmat,
      Age = agemat,
      st = len[, 'st_f'],
      nStrata = nStrata,
      minSel = minSel,
      maxSel = maxSel,
      selType = s,
      nmat = nmat
    )

  # Now estimate everything
  map <- NULL
  model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)
  fit <- nlminb(model$par, model$fn, model$gr, control = list(rel.tol = 1e-12, eval.max = 100000,iter.max = 10000), lower = 0.01)
  # for (k in 1:3)  fit <- nlminb(model$env$last.par.best, model$fn, model$gr) ## start at last-best call, for stability
  
  best <- model$env$last.par.best
  print(c(best[1:6],exp(best[7:length(best)])))
  print(fit$objective)
  rep <- sdreport(model)
  dat0 <- bind_cols(dat0, model$report()$ypreds %>% data.frame()) ## each 6 cols is new sim
  
}
AAA
## plotting ----
# names(dat0) <- paste0(unique(len$st_f))
# names(dat0) <- paste0(c(rep('uniform_',6),rep('min_',6),rep('max_',6)),paste0(unique(len$st_f)))
names(dat0) <- paste0(c(rep('uniform_',6),rep('dome_',6)),paste0(unique(len$st_f)))
agedf <- agemat  %>% data.frame() %>% reshape::melt()  %>% plyr::rename(c("value" = "Age")) %>% select(-variable) 
lendf <- lenmat %>% data.frame() 
# names(lendf) <-  paste0(rep('Length_',6),paste0(unique(len$st_f)))
names(lendf) <- paste0(rep('Length_',6),c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s"))

mdf<- bind_cols(dat0,  lendf) %>% reshape::melt(id = c())  %>% plyr::rename(c( 'variable' = 'ID', 'value' = 'Length')) %>% 
  mutate(model = sub('_.*$','', ID), st = sub(".*_ *(.*?)", "\\1", ID), age = rep(agedf$Age,3))  %>% select(-ID)
## for plotting
mdf$st_f = factor(mdf$st, levels=c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s"))
# save(mdf,file =  paste0(getwd(),"/uniform_predicts.rda"))


ggplot(len, aes(x = Age, y = Length_cm)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = c(0.9,0.15))+
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  # scale_color_brewer(palette = 'Dark2') +
  geom_point(alpha = 0.2) +
  geom_line(data = subset(mdf, model != 'Length'), aes(x = age, y = Length, color = model), lwd = 1.1) +
  facet_wrap(~ st_f) +
  labs(title = "Predicted Model Fits and Raw Data", 
       y = 'Length (cm)', x= 'Age (yr)', 
       subtitle = 'Parameters estimated using subset of each strata, not sex-specific',
       color = 'selectivity model')

# ggsave(file = paste0(getwd(),"/plots/dome_uni_fits2.png"), plot = last_plot(), height = 5, width = 7, unit = 'in', dpi = 520)

# ggplot(subset(mdf, model == 'Length'), aes(x = age, y = Length)) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(), legend.position = c(0.9,0.15))+
#   # scale_y_continuous(limits = c(0,100)) +
#   # scale_x_continuous(limits = c(0,50)) +
#   geom_point(alpha = 0.2) +
#   geom_line(data = subset(mdf, model != 'Length'), aes(x = age, y = Length, color = model), lwd = 1.1)+
#   facet_wrap(~ st_f) +
#   labs(
#     title = "Predicted Model Fits and Subsampled Data",
#     y = 'Length (cm)',
#     x = 'Age (yr)',
#     subtitle = 'Points are actual subsampled data used in parameter fitting',
#     color = 'selectivity model'
#   )
# ggsave(file = paste0(getwd(),"/plots/dome_uni_fits1.png"), plot = last_plot(), height = 5, width = 7, unit = 'in', dpi = 520)
