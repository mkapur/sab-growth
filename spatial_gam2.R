require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)
library(ggsidekick);require(TMB);require(MASS)
## load and cbind data
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")
load( paste0(getwd(),"/data/all_data.rda")) ## all_data
load( paste0(getwd(),"/data/DES.rda")) ## DES
load( paste0(getwd(),"/data/KEY.rda")) ## KEY

compile("sptlvb.cpp")
dyn.load(dynlib("sptlvb"))

dat0 <- rep0 <- aic0 <- NULL ## later storage
## get K, Linf estimates for ALL data
nStrata <- length(unique(DES[,1])) ## pooled hypothesis
data <-
  list(
    Length_cm = all_data[,"Length_cm"],
    Age = all_data[,"Age"],
    DES = as.vector(DES[,1]),
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
rep0 <- bind_rows(rep0, bind_cols(data.frame(names(rep$value)),data.frame(rep$value),
                                  data.frame(rep$sd),data.frame(c(rep(unique(KEY[,1]),3), rep("ALL",1)))))
aic0 <- c(aic0, model$report()$aic %>% data.frame())
names(rep0) <- c('variable', 'value','sd', 'ID')
rep0 <- rep0 %>% mutate(
  Sex = sub('_.*$', '', ID),
  st = sub(".*_ *(._?)", "\\1", ID)
) %>% select(-ID)


parest <- rep0
## exponentiate logk
parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))
write.csv(parest, paste0(getwd(),"/results/parEstPOOLED_",Sys.Date(),'.csv'),row.names = F)
parest <- read.csv("C:/Users/mkapur/Dropbox/UW/sab-growth/results/parEstPOOLED_2018-12-09.csv") %>% filter(variable != "Sigma")

## extract covariance matrix for parests 
covMat <- cov2cor(rep$cov.fixed)[3:6,3:6] %>% data.frame()
rownames(covMat) <- NULL
## use mvrnorm to simulate 10k values of Linf and K
simd <- mvrnorm(n = 10E4, mu =  c(parest$value[3:6]),
                Sigma = covMat, empirical = FALSE) %>% data.frame()
names(simd) <-  c(paste(parest$variable[3:6],
                        parest$Sex[3:6]))
