## given a GAM model, identify breakpoints via Deriv. Only care about latitude breaks.


getBreaks <- function(gammod = mod, dat, scenario = scen ){
  
  
  ## get & eval derivatives ----
  llsmooth <- mod$smooth[2][[1]]
  # llsmooth$knots
  ## calc first derivatives
  # pdat <- sample_n(dat,0.9*nrow(dat))
  pdat <- dat
  pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
  p2 <- predict(mod, newdata = pdat) ## raw predicts
  pdat <- transform(pdat, predLen = p2, se2_spt = pTerm$se.fit[,2], se2_yr = pTerm$se.fit[,1])
  # pdat <- transform(pdat, predLen = p2, se2_lat = pTerm$se.fit[,2],se2_lon = pTerm$se.fit[,3], se2_yr = pTerm$se.fit[,1])
  
  
  df.res <- df.residual(mod)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  ## variances are additive just FYI
  pdat <- transform(pdat,
                    upper = predLen + (crit.t * (se2_spt+se2_yr)),
                    lower = predLen - (crit.t * (se2_spt+se2_yr)))
  # pdat <- transform(pdat,
  #                   upper = predLen + (crit.t * (se2_spt)),
  #                   lower = predLen - (crit.t * (se2_spt)))
  
  ## source from gist
  source(paste0(getwd(),"/Deriv.R"))
  
  ## here I iterate sex as a sanity check but it is actually not necessary -- results will be the same since sex is a covariate
  # breaksdf <- data.frame(Sex = NA, YrBreaks = NA, LatBreaks = NA)
  breaksdf <- list()
  # idx <- 1
  # for(s in 1:2){
    # for(t in 1:2){
      
      # newD <- data.frame(Age = rep(6,100),
      #   Year = seq(1,50,length = 100),
      #                    Latitude_dd = seq(0,50,length = 100))
  
     newD <- data.frame(Year = seq(1,50,length = 100),
                     Latitude_dd = seq(0,50,length = 100))
      
      dr <- Deriv(mod,newdata = newD)
      confint.Deriv(dr)
      
      Terms <- c("Year","Latitude_dd")[2]
      
      m2.d <- Deriv(mod, newdata = newD)
      m2.dci <- confint(m2.d, term = Term)
      
      # crit.eval = mean(m2.d[[Term]]$deriv) ## use mean
      crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
      crit.eval.se = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$se.deriv) ## use tails
      

      ## identify where mean crosses zero or falls out of bounds (NAs where it does)
      m2.dsig.zeros <- signifMK(x = m2.d$eval[[Term]],
                                d = m2.d[[Term]]$deriv,
                                upper = m2.dci[[Term]]$upper,
                                lower = m2.dci[[Term]]$lower,
                                crit.eval = crit.eval,
                                eval = 0)
      pix <- !is.na(m2.dsig.zeros)
      
      vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
      breaksdf[[1]] <- sort(c(unique(round(vals)))) ## get rounded unique
      
      
    
     plist <- list()
     idx <- 1
     # ## extract plotting stuff
     pd <- plot(mod,   select = t, scheme  =2, lwd  =2, main = paste0(Terms[t],' Smoother'), cex.axis = 2, ylim = c(-10,ifelse(t != 3,10,500)))
     for(t in 1:length(Terms)){
       temp0 <- pd[t] ## get this smoother
       temp <- data.frame(cbind(temp0[[1]]$x,temp0[[1]]$fit, temp0[[1]]$se)); names(temp) = c('x','fit','se')
       plist[[idx]] <-  ggplot(temp, aes(x = x, y = fit)) +
         theme_minimal() +
         theme(panel.grid = element_blank())+
         geom_line(lwd = 1.1) +
         geom_line(aes(y= fit-se), linetype = 'dashed') +
         geom_line(aes(y= fit+se), linetype = 'dashed') +
         geom_rug(sides = 'b') +
         labs(x = Terms[t], y = "smoother", title = paste0('Smoother for ', temp0[[1]]$xlab))
       idx <- idx+1
       CI <- confint(m2.d, term = Terms[t])
       m2.dtemp <- data.frame(cbind(m2.d$eval[,Terms[t]], m2.d[[Terms[t]]]$deriv, CI[[1]]$upper, CI[[1]]$lower)); names(m2.dtemp) = c('x','deriv','upper','lower')
       plist[[idx]] <-  ggplot(m2.dtemp,aes(x = x, y = deriv))    +
         theme_minimal() +
         theme(panel.grid = element_blank())+
         geom_line(lwd = 1.1) +
         geom_hline(yintercept = 0, col = 'grey22') +
         geom_line(aes(y= upper), linetype = 'dashed') +
         geom_line(aes(y= lower), linetype = 'dashed') +
         geom_vline(xintercept = breaksdf[[t]], col = 'red') +
         labs(x = Terms[t], y = "f'(x)", title = paste0('First Derivative for ', temp0[[1]]$xlab))
       idx <- idx+1
     }
     
     lay <- rbind(c(2,3,1,1))
     
     grid.arrange(grobs = plist, layout_matrix = lay)  %>% 
       ggsave(plot = .,  file = paste0(outdir,"/",b,"_gam_smooths.png"), width = 11, height = 8, units = 'in', dpi = 480)
     
  
  # png(file = paste0(outdir,"/",b,"_gam_smooths.png"), height = 6, width = 8, units = 'in', res = 500)
  
  # layout(matrix(1:2, ncol = 2))
  
  
  # plot(mod,  select  = 1,  scheme    =2,  lwd  =2, main = 'Year Smoother', cex.axis = 2)
  # plot(mod,   select = 2, scheme  =2, lwd  =2, main = 'Latitude Smoother', cex.axis = 2)
  
  
  # plot.Deriv(m2.d, term = 'Year', cex.axis = 2, main = 'derivative of year')
  # abline(v = breaksdf[[1]], col = 'red')
  
  # plot.Deriv(m2.d, term = 'Latitude_dd',cex.axis = 2,  main = 'derivative of latitude')
  # abline(v = breaksdf[[1]], col = 'red')
  # 
  # graphics.off()
  # 
  # png(file = paste0(outdir,"/",b,"_gamAll.png"), height = 6, width = 8, units = 'in', res = 500)
  # layout(matrix(1:6, ncol = 2,byrow = T))
  # gam.check(mod)
  # 
  # plot(mod,   select = 2, scheme  =2, lwd  =2, main = 'Latitude Smoother', cex.axis = 2)
  # plot.Deriv(m2.d, term = 'Latitude_dd',cex.axis = 2,  main = 'derivative of latitude')
  # abline(v = breaksdf[[1]], col = 'red')
  graphics.off()
  return(breaksdf)
}
