require(TMB)

## load data
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")
load( paste0(getwd(),"/filtered_SAB_WC.rda")); load(paste0(getwd(),"/lenmat_WC.rda")); load(paste0(getwd(),"/agemat_WC.rda")) ## made in dataprep
nStrata <- length(unique(len$st_f))

## Estimation
data <- list(Length_cm = lenmat, Age = agemat, st = len[,'st_f'], nStrata = nStrata )
parameters <- list(dummy = 0, log_Linf = rep(log(70), nStrata), log_k = rep(0, nStrata), log_t0 = rep(0, nStrata), log_Sigma = 0.1)

compile("sptlvb.cpp")
dyn.load(dynlib("sptlvb"))
# Now estimate everything
map <- NULL
model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)

# Bounds on the parameters MAKE SURE THESE MATCH THE ORDER AND SIZE OF EACH PARAM
# lowbnd = c( -10, -7, rep(0.1,m),rep(-7,m), rep(-Inf,m))
# uppbnd = c( 10, 4, rep(1000,m),rep(4,m), rep(Inf,m))

fit <- nlminb(
  model$par,
  model$fn,
  model$gr,
  control = list(
    rel.tol = 1e-12,
    eval.max = 100000,
    iter.max = 1000
  )
  # ,
  # lower = lowbnd,
  # upper = uppbnd
)
best <- model$env$last.par.best
print(best)
rep <- sdreport(model)
print(summary(rep))
model$report()
