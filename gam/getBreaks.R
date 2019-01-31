## given a GAM model, identify breakpoints via Deriv


getBreaks <- function(gammod = mod, dat, scenario = scen ){
  
  
  ## get & eval derivatives ----
  llsmooth <- mod$smooth[2][[1]]
  llsmooth$knots
  ## calc first derivatives
  pdat <- sample_n(dat,1000)
  pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
  p2 <- predict(mod, newdata = pdat) ## raw predicts
  pdat <- transform(pdat, predLen = p2, se2_spt = pTerm$se.fit[,3], se2_yr = pTerm$se.fit[,2])
  df.res <- df.residual(mod)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  ## variances are additive just FYI
  pdat <- transform(pdat,
                    upper = predLen + (crit.t * (se2_spt+se2_yr)),
                    lower = predLen - (crit.t * (se2_spt+se2_yr)))
  
  ## source from gist
  source(paste0(getwd(),"/Deriv.R"))
  
  ## here I iterate sex as a sanity check but it is actually not necessary -- results will be the same since sex is a covariate
  # breaksdf <- data.frame(Sex = NA, YrBreaks = NA, LatBreaks = NA)
  breaksdf <- list()
  idx <- 1
  for(s in 1:2){
    # breaksdf[s,'Sex'] = c("M","F")[s]
    for(t in 1:2){
      
      newD <- data.frame(Age = mean(seq(0,94,length = 100)), ## use mean age?
                         Sex = c("M","F")[s], 
                         Year = seq(1981,2017,length = 100),
                         Longitude_dd = seq(-186,-116, length = 100),
                         Latitude_dd = seq(0,64,length = 100))
      
      dr <- Deriv(mod,newdata = newD)
      confint.Deriv(dr)
      
      Term <- c("Year","Latitude_dd")[t]
      
      m2.d <- Deriv(mod, newdata = newD)
      m2.dci <- confint(m2.d, term = Term)
      
      # crit.eval = mean(m2.d[[Term]]$deriv) ## use mean
      crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
      ## identify where CI does NOT include some crit value (return NA where it does)
      m2.dsig <- signifD(m2.d$eval[[Term]], 
                         d = m2.d[[Term]]$deriv,
                         m2.dci[[Term]]$upper, 
                         m2.dci[[Term]]$lower, eval = crit.eval)
      pix <- c(which(!is.na(m2.dsig$incr)), which(!is.na(m2.dsig$decr)))
      vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
      breaksdf[[idx]] <- sort(c(unique(round(vals)))) ## get rounded unique
      # cat(idx)
      idx <- idx+1
      
    }
  }
  
  png( file = paste0("C:/users/maia kapur/dropbox/uw/sab-growth/gam/plots/",scenario,"gam_smooths.png"), height = 6, width = 8, units = 'in', res = 500)
  layout(matrix(1:4, ncol = 2))
  plot(mod,  select  =1,  scheme  
       =2,  lwd  =2, main = 'Year Smoother', cex.axis = 2)
  plot(mod,  select  =2,  scheme  =2, xlim = c(35,60),
       lwd  =2, main = 'Latitude Smoother', cex.axis = 2)
  plot.Deriv(m2.d, term = 'Year', cex.axis = 2, main = 'derivative of year')
  abline(v = breaksdf[[1]], col = 'red')
  plot.Deriv(m2.d, term = 'Latitude_dd',
             cex.axis = 2, xlim = c(35,60), main = 'derivative of latitude')
  abline(v = breaksdf[[2]], col = 'red')
  
  dev.off()
  
  
}
