
fitMod <- function(data, parameters, modversion = "sptlschnute", map){
  dat0 <- rep0 <- NULL ## later storage
  
  model <- MakeADFun(data, parameters,  DLL=modversion,silent=T,map=map)
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
  rep <- sdreport(model)
  
  dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
  
  rep0 <- bind_rows(rep0,
                    bind_cols(
                      data.frame(names(rep$value)),
                      data.frame(rep$value),
                      data.frame(rep$sd),
                      
                      data.frame(c(rep(keybase
                                       , 4), rep("ALL", 1)))))
  
  
  names(rep0) <- c('variable', 'value','sd', 'REG')
  
  rep0$value[rep0$variable %in% c('log_k','log_Lone','log_Ltwo')] <- exp(rep0$value[rep0$variable %in% c('log_k','log_Lone','log_Ltwo')])
  rep0$sd[rep0$variable %in% c('log_k','log_Lone','log_Ltwo')] <- exp(rep0$sd[rep0$variable %in% c('log_k','log_Lone','log_Ltwo')])
  rep0$variable <- factor(rep0$variable, levels = c("k","log_Ltwo","log_Lone","Sigma","t0","log_k","log_Linf","Lone","Ltwo")) ## enable new levels
  rep0$variable[rep0$variable == 'log_k'] <- 'k'
  rep0$variable[rep0$variable == 'log_Lone'] <- 'Lone'
  rep0$variable[rep0$variable == 'log_Ltwo'] <- 'Ltwo'
  
  return(list(dat0,rep0))
}
