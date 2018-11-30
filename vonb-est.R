require(TMB); require(dplyr); require(ggplot2); require(reshape)
options(scipen=999)
rm(list = ls())

## load objects made in dataprep
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")
load( paste0(getwd(),"/data/all_data.rda")) ## all_data
load( paste0(getwd(),"/data/DES.rda")) ## DES
# load( paste0(getwd(),"/data/REG.rda")) ## REG
# load( paste0(getwd(),"/data/SEL.rda")) ## SEL
load( paste0(getwd(),"/data/KEY.rda")) ## KEY

compile("sptlvb.cpp")
dyn.load(dynlib("sptlvb"))

dat0 <- rep0 <- aic0 <- NULL ## later storage
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
        log_Linf = rep(log(70), nStrata),
        log_k = rep(log(0.5), nStrata),
        t0 = rep(0, nStrata),
        log_Sigma = 0
      )
    
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
        iter.max = 10000
      )
    )
    # for (k in 1:3)  fit <- nlminb(model$env$last.par.best, model$fn, model$gr) ## start at last-best call, for stability
    
  best <- model$env$last.par.best

  rep <- sdreport(model)
  dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) ## each 6 cols is new sim
  rep0 <- bind_rows(rep0, bind_cols(data.frame(names(rep$value)),data.frame(rep$value),data.frame(rep$sd),data.frame(c(rep(unique(KEY[,s]),3), rep("ALL",1)))))
  aic0 <- c(aic0, model$report()$aic %>% data.frame())
}

## reformat outputs ----
aic <- aic0 %>% data.frame()
names(aic) <- c('Pooled','3 Zone', 'Lit', 'Survey Strata')
aic[which.min(aic)]

names(rep0) <- c('variable', 'value','sd', 'ID')
rep0 <- rep0 %>% mutate(
  model = c(
    rep('Pooled',length(unique(DES[,1]))*3 + 1),
    rep('3 Zone', length(unique(DES[,2]))*3 + 1),
    rep('Lit', length(unique(DES[,3]))*3 + 1),
    rep('Survey Strata', length(unique(DES[,4]))*3 + 1)
  ),
  Sex = sub('_.*$', '', ID),
  st = sub(".*_ *(._?)", "\\1", ID)
) %>% select(-ID)

write.csv(rep0, paste0(getwd(),"/results/parEst_",Sys.Date(),'.csv'),row.names = F)

ypreds0 <- cbind(dat0,all_data) %>% data.frame()  
names(ypreds0)[1:6] <- c('Pooled','3 Zone', 'Lit', 'Survey Strata', 'Age', 'Observed_Length')
ypreds <- ypreds0 %>% melt(id = c("Age","Sex","Observed_Length", "H2","H3","H4"))  %>% 
  mutate('REGION' = ifelse(variable == 'Pooled',"ALL",
                           ifelse(variable == '3 Zone',H2, 
                                  ifelse(variable == 'Lit',H3,
                                         H4)))) %>% select(-H2,-H3,-H4)
write.csv(ypreds, paste0(getwd(),"/results/predicts_",Sys.Date(),'.csv'),row.names = F)

## plotting ----

## fits
for(vr in unique(ypreds$variable)){
ggplot(subset(ypreds, variable == vr), aes(x = Age, y = value, col = Sex )) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9,0.1),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  scale_color_manual(values = c("#d8b365","#5ab4ac"))+
    scale_alpha(guide = 'none') +
  geom_point(alpha = 0.2, aes(y = Observed_Length)) +
  geom_line(lwd = 1.1) +
  labs(y = 'Length (cm)', title = vr, col = "") +
  facet_wrap(~ REGION) 
  ggsave(file = paste0(getwd(),"/plots/fits_",vr,".png"), 
         plot = last_plot(), height = 11, width = 17, unit = 'in', dpi = 520)
}
  
  
## parest by sex and region
parest <- read.csv("C:/Users/mkapur/Dropbox/UW/sab-growth/results/parEst_2018-11-21.csv") %>% filter(variable != "Sigma")
## exponentiate logk
parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))

for(vr in unique(ypreds$variable)){
  ggplot(subset(parest, model == vr), aes(x = st, y = value, col = Sex))+
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = c(0.25,0.9),
          legend.text = element_text(size = 14),
          strip.text = element_text(size=14))+
    scale_color_manual(values = c("#d8b365","#5ab4ac"))+
    geom_point() +
    geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    labs(x = 'Spatial Stratum', y = "") +
    facet_wrap(~model+variable, scales = "free_y") 
  
  ggsave(file = paste0(getwd(),"/plots/parplot_",vr,".png"), 
         plot = last_plot(), height = 6, width = 8, unit = 'in', dpi = 520)
}
