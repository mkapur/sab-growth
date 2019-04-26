
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
                                       , 5), rep("ALL", 1)))))
  
  
  names(rep0) <- c('variable', 'value','sd', 'cREG')
  rep0$REG <- gsub("_.*", "\\1", rep0$cREG)
  rep0$cREG <- gsub(".*._", "\\1", rep0$cREG)
  
  # rep0$value[rep0$variable %in% c('log_k','log_Linf')] <- exp(rep0$value[rep0$variable %in% c('log_k','log_Linf' )])
  # rep0$sd[rep0$variable %in% c('log_k','log_Linf')] <- exp(rep0$sd[rep0$variable %in% c('log_k','log_Linf')])
  ## add in bias correction to estimates
  rep0$value[rep0$variable %in% c('log_k','log_Linf')] <- exp(rep0$value[rep0$variable %in% c('log_k','log_Linf' )]-(rep0$sd[rep0$variable %in% c('log_k','log_Linf')]^2)/2)
  rep0$sd[rep0$variable %in% c('log_k','log_Linf')] <- exp(rep0$sd[rep0$variable %in% c('log_k','log_Linf')])
  
  rep0$variable <- factor(rep0$variable, levels = c("k","log_Ltwo","Linf","Sigma","t0","log_k","log_Linf","L1","L2")) ## enable new levels
  rep0$variable[rep0$variable == 'log_k'] <- 'k'
  rep0$variable[rep0$variable == 'log_Linf'] <- 'Linf'

  
  return(list(dat0,rep0))
}
