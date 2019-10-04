require(TMB)
require(dplyr)


compile("it10.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("it10"))

# compile("integTest.cpp") ## will throw 'incomplete line' msg, ignore
# dyn.load(dynlib("integTest0"))

data <-
  list(
    Length_cm = runif(9, 10, 15),
    Sel = runif(9, 0, 1),
    Age = 1:9,
    nStrata = 3,
    selType = sample(0:1, 9, replace = T),
    DES = runif(9, 0, 1),
    a2 = 30
  )

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


parameters <- list( log_Sigma = 0.1, log_Linf = rep(80,3), t0 = 2:5, log_k =  log(2:4))
map <- NULL
model <- MakeADFun(data, parameters,  DLL="it10",silent=T,map=map)
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
