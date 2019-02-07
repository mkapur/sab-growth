## given a GAM model, identify breakpoints via Deriv


getBreaks <- function(gammod = mod, dat, scenario = scen ){
  
  
  ## get & eval derivatives ----
  llsmooth <- mod$smooth[2][[1]]
  # llsmooth$knots
  ## calc first derivatives
  pdat <- sample_n(dat,100)
  # pdat <- dat
  pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
  p2 <- predict(mod, newdata = pdat) ## raw predicts
  pdat <- transform(pdat, predLen = p2, se2_spt = pTerm$se.fit[,2], se2_yr = pTerm$se.fit[,1])
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
  # for(s in 1:2){
    # breaksdf[s,'Sex'] = c("M","F")[s]
    for(t in 1:2){
      
      newD <- data.frame(Year = seq(1,50,length = 100),
                         Latitude_dd = seq(0,50,length = 100))
      
      dr <- Deriv(mod,newdata = newD)
      confint.Deriv(dr)
      
      Term <- c("Year","Latitude_dd")[t]
      
      m2.d <- Deriv(mod, newdata = newD)
      m2.dci <- confint(m2.d, term = Term)
      
      # crit.eval = mean(m2.d[[Term]]$deriv) ## use mean
      crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
      crit.eval.se = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$se.deriv) ## use tails
      
      ## identify where falls outside CI (NAs where inside CI)
      # m2.dsig.outliers <- signifD(m2.d$eval[[Term]],
      #                    d = m2.d[[Term]]$deriv,
      #                    m2.dci[[Term]]$upper,
      #                    m2.dci[[Term]]$lower,
      #                    eval = crit.eval)

      ## identify where mean crosses zero or falls out of bounds (NAs where it does)
      m2.dsig.zeros <- signifMK(x = m2.d$eval[[Term]],
                                d = m2.d[[Term]]$deriv,
                                upper = m2.dci[[Term]]$upper,
                                lower = m2.dci[[Term]]$lower,
                                crit.eval = crit.eval,
                                eval = 0)
      
      ## of passing ones, find outliers
      
      ## identify where CI crosses zero or falls out of bounds (NAs where it does)
      # m2.dsig.CIs <- signifCI(x = m2.d$eval[[Term]],
      #                         upper = m2.dci[[Term]]$upper,
      #                         lower = m2.dci[[Term]]$lower,
      #                          eval = 0)
      
      
      ## only retain vals that pass both crit
      # outkeep <- !is.na(m2.dsig.outliers$incr)[!is.na(m2.dsig.outliers$decr)]
      # zkeep <- !is.na(m2.dsig.outliers$incr)[!is.na(m2.dsig.outliers$decr)]
      # ikeep <- !is.na(m2.dsig.zeros$incr)[!is.na(m2.dsig.CIs$incr)]
      # dkeep <- !is.na(m2.dsig.zeros$decr)[!is.na(m2.dsig.CIs$decr)]
      # 

      pix <- !is.na(m2.dsig.zeros)
     # pix <- outkeep[zkeep]
      # pix <- zkeep[cikeep]
     vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
     breaksdf[[idx]] <- sort(c(unique(round(vals)))) ## get rounded unique
     # cat(idx)
     idx <- idx+1
     
    }
  
  
  png(file = paste0(outdir,"/",b,"_gam_smooths.png"), height = 6, width = 8, units = 'in', res = 500)
  
  layout(matrix(1:4, ncol = 2))
  
  
  plot(mod,  select  = 1,  scheme    =2,  lwd  =2, main = 'Year Smoother', cex.axis = 2)
  plot(mod,   select = 2, scheme  =2, lwd  =2, main = 'Latitude Smoother', cex.axis = 2)
  
  
  plot.Deriv(m2.d, term = 'Year', cex.axis = 2, main = 'derivative of year')
  abline(v = breaksdf[[1]], col = 'red')
  
  
  plot.Deriv(m2.d, term = 'Latitude_dd',cex.axis = 2,  main = 'derivative of latitude')
  abline(v = breaksdf[[2]], col = 'red')
  
  graphics.off()
  
  return(breaksdf)
}
