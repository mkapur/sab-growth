makeStrata <- function(mod)

## now re-aggregate the data ----
## generate DES matrix of vectors and a KEY for later comparison
DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(all_data)) ## H1 is all zeros
## get region
# temp0 <- with(all_data, ifelse(Latitude_dd >= 49 & Year < 2005 ,0,
#                                ifelse(Latitude_dd >= 49 & Year >= 2005,1,
#                                       ifelse(Latitude_dd < 49 & Year < 2005,2,3))))
# ## get region and sex combos as factors
# DES <- as.numeric(factor(paste0(temp0,all_data[,'Sex'])))-1
# KEY <- paste0(c("North_Early","North_Late","South_Early","South_Late")[(temp0+1)],"_",all_data[,'Sex'])
keybase <- apply(expand.grid(c("North_Early","North_Late","South_Early","South_Late"), c('F','M')), 1,
                 paste, collapse="_")

# save(DES, file = paste0(getwd(),"/data/DES_gam.rda"))
# save(KEY, file = paste0(getwd(),"/data/KEY_gam.rda"))
load( paste0(getwd(),"/data/DES_gam.rda")) ## DES
load( paste0(getwd(),"/data/KEY_gam.rda")) ## KEY

## run TMB parest----
compile("sptlvb.cpp")
dyn.load(dynlib("sptlvb"))

dat0 <- rep0 <- NULL ## later storage
nStrata <- length(unique(DES))
data <-
  list(
    Length_cm = all_data[,"Length_cm"],
    Age = all_data[,"Age"],
    DES = as.vector(DES),
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
rep0 <- bind_rows(rep0,
                  bind_cols(
                    data.frame(names(rep$value)),
                    data.frame(rep$value),
                    data.frame(rep$sd),
                    # data.frame(c(rep(unique(
                    #   KEY
                    # ), 3), rep("ALL", 1)))
                    data.frame(c(rep(keybase
                                     , 3), rep("ALL", 1)))
                    
                  ))

## reformat outputs ----
names(rep0) <- c('variable', 'value','sd', 'ID')
rep0 <- rep0 %>% mutate(
  st =sub('_([^_]*)$', '',  ID),
  Sex = sub(".*_ *(._?)", "\\1", ID)
) %>% select(-ID)

write.csv(rep0, file = paste0(paste0(getwd(),"/results/parEst_gam_",Sys.Date(),'.csv')),row.names = F)

ypreds0 <- cbind(dat0,all_data) %>% data.frame()  
names(ypreds0)[1] <- c('Predicted')
ypreds0$REG <- as.factor(sub('_([^_]*)$', '', KEY))

write.csv(ypreds0,  paste0(getwd(),"/results/predicts",Sys.Date(),".csv"))
ypreds <- read.csv( paste0(getwd(),"/results/predicts",Sys.Date(),".csv"))
## plotting ----

## fits

ggplot(ypreds, aes(x = Age, y = Predicted, col = Sex )) +
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
  geom_point(alpha = 0.2, aes(y = Length_cm)) +
  geom_line(lwd = 1.1) +
  labs(y = 'Length (cm)', col = "") +
  facet_wrap(~REG)
ggsave(file = paste0(getwd(),"/plots/fits_gam",Sys.Date(),".png"), 
       plot = last_plot(), height = 8, width = 10, unit = 'in', dpi = 520)


## parest by sex and region
parest <- read.csv( paste0(getwd(),"/results/parEst_gam_",Sys.Date(),".csv")) %>%
  filter(variable != "Sigma")
## exponentiate logk
parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))

ggplot(parest, aes(x = st, y = value, col = Sex))+
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.25,0.9),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14),
        axis.text = element_text(angle = 45))+
  scale_color_manual(values = c("#d8b365","#5ab4ac"))+
  geom_point() +
  geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
  labs(x = 'Spatial Stratum', y = "") +
  facet_wrap(~variable, scales = "free_y") 

ggsave(file = paste0(getwd(),"/plots/parplot_gam",Sys.Date(),".png"), 
       plot = last_plot(), height = 6, width = 8, unit = 'in', dpi = 520)

