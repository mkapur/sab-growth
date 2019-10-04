require(TMB)
require(dplyr)

compile("sptlvbSel.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvbSel"))

compile("integTest0.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("integTest0"))

# compile("integTest.cpp") ## will throw 'incomplete line' msg, ignore
# dyn.load(dynlib("integTest0"))

data <- list( n = 2:4, A = matrix(1:9, ncol =3))
parameters <- list(b = rep(0.1,3), logSigma = 0.1, 
                   log_Linf = rep(80,3), t0 = 11:13)
map <- NULL
model <- MakeADFun(data, parameters,  DLL="integTest0",silent=T,map=map)
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

best <- model$env$last.par.best

# rep <- sdreport(model, bias.correct = TRUE)
rep <- sdreport(model)

# // Type operator()           // Evaluate integrand
# //   (vector<Type> x){
  # //   return 1/(2*3.14*Sigma *Age(i))*exp(-(yfit - theta)/2*(Age(i)*Sigma)) * selPred

## this will assign a unique DES depending on period X sex X region -- whatever is in DES
dat <- sample_n(dat, size = nrow(dat)*0.25)
data <-
  list(
    Length_cm = dat[,"Length_cm"],
    Age = dat[,"Age"],
    DES = as.vector(DES), ## keep this for master iterations
    selType = dat[,'selType'],
    Sel = dat[,'Sel'],
    nStrata = nStrata,
    a2 = 15
  )

parameters <-
  list(
    log_Linf = rep(log(70), nStrata),
    log_k = rep(log(0.5), nStrata),
    t0 = rep(0.1, nStrata),
    log_Sigma = 0.1
  )

# Now estimate everything
map <- NULL
model <- MakeADFun(data, parameters,  DLL="integTest",silent=T,map=map)
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


dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
rep0 <- bind_rows(rep0,
                  bind_cols(
                    data.frame(names(rep$value)),
                    data.frame(rep$value),
                    data.frame(rep$sd),
                    data.frame(c(rep(keybase, 5), rep("ALL", 1)))))