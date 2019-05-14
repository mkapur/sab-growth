## BURN IT DOWN COMPARE BREAKS
compareNA <- function(v1,v2) { ## ADD RELAXED CRITERA plus/minus one
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}


## tabulate whether the detection & estimation methods were accurate.
parms <- c('L1','L2','k','Linf')[1:2]

matchdf <- read.csv('./input_data/matchdf3.csv')

parest <-  rep0 %>% 
  filter(variable %in% parms) %>% 
  mutate(source = 'Estimated')

parest$lwr <- parest$value - 1.96*parest$sd
parest$upr <- parest$value + 1.96*parest$sd

## subset & compare based on match.df
for(n in unique(parest$REG)){
  
  # if(scen != 'tempvar_R1R2' & is.na(breaksdf$yr_breaks)){
  if(is.na(breaksdf$yr_breaks)){
    
    cdf[idx,'scen'] <- scen
    cdf[idx,'boot'] <- b
    # cdf[idx,'REG'] <- t.REG ## if anything other than R1, will compare to appropriate RX (2 or 3)
    # # cdf[idx,'REG'] <- ifelse(factor(gsub("_.*", "\\1", gamREGS[iq]),levels = Rlevs) %in% REGS, REGS[REGS == factor(gsub("_.*", "\\1", gamREGS[iq]),levels = Rlevs)],
    
    cdf[idx,'gamLAT'] <- breaksdf$lat_breaks
    cdf[idx,'gamLON'] <- breaksdf$lon_breaks
    cdf[idx,'gamYR'] <- breaksdf$yr_breaks
    # 
    cdf[idx,"LAT"] <- compareNA(cdf[idx,"gamLAT"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
    cdf[idx,"LON"] <- compareNA(cdf[idx,"gamLON"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
    
    # if(scen == "F0LMW")  { # won't match "overlap"
    #   cdf[idx,"LAT"] <- cdf[idx,"gamLAT"] >= 20 & cdf[idx,"gamLAT"] <= 25
    #   cdf[idx,"LON"] <- cdf[idx,"gamLON"] >= 20 & cdf[idx,"gamLON"] <= 25
    # }
    if(scen == "F0LMW")  { # won't match "overlap"
      cdf[idx,"LAT"] <-  ifelse(!is.na(cdf[idx,'gamLAT']), cdf[idx,"gamLAT"] >= 20 & cdf[idx,"gamLAT"] <= 25,FALSE)
      cdf[idx,"LON"] <-  ifelse(!is.na(cdf[idx,'gamLON']), cdf[idx,"gamLON"] >= 20 & cdf[idx,"gamLON"] <= 25,FALSE)
    }
    cdf[idx,"YEAR"] <- compareNA(cdf[idx,"gamYR"],scenarios$TEMPORAL[scenarios$DESC == cdf[idx,'scen']][1])
    
    ptf <- ptf.miss <- c(NA,NA)
    
    for(v in  1:length(parms)){ ## loop variables
      tmp.true <- subset(matchdf, nscen == scen & guessREG == n & variable == parms[v] & guessCREG == 'pooled')
      # tmp.par <- parest %>% filter(REG == n & variable == parms[v])
      tmp.par <- parest %>% filter(REG == n & variable == parms[v] & source == 'Estimated')
      bounds <- tmp.par %>% select(lwr,upr)
      ptf[v] <- tmp.true$value >= bounds[1] & tmp.true$value <= bounds[2]
      ptf.miss[v] <- missby(ruling = ptf[v], bounds = bounds, truth = tmp.true)
      
    } ## end parms
    cdf[idx,'REG'] <- tmp.true$trueREG
    cdf[idx,'gamREG'] <- n
    cdf[idx,parms[1]] <- ptf[1];  cdf[idx,parms[2]] <- ptf[2]
    cdf[idx,'L1_MISS'] <- ptf.miss[1];cdf[idx,'L2_MISS'] <- ptf.miss[2]; idx <- idx + 1
    
    # } else if(scen == 'tempvar_R1R2' | !is.na(breaksdf$yr_breaks)){
  }else if( !is.na(breaksdf$yr_breaks)){
    for(c in unique(parest$cREG)){ ## also loop periods
      
      cdf[idx,'scen'] <- scen
      cdf[idx,'boot'] <- b
      # cdf[idx,'REG'] <- t.REG ## if anything other than R1, will compare to appropriate RX (2 or 3)
      # # cdf[idx,'REG'] <- ifelse(factor(gsub("_.*", "\\1", gamREGS[iq]),levels = Rlevs) %in% REGS, REGS[REGS == factor(gsub("_.*", "\\1", gamREGS[iq]),levels = Rlevs)],
      
      cdf[idx,'gamLAT'] <- breaksdf$lat_breaks
      cdf[idx,'gamLON'] <- breaksdf$lon_breaks
      cdf[idx,'gamYR'] <- breaksdf$yr_breaks
      # 
      cdf[idx,"LAT"] <- compareNA(cdf[idx,"gamLAT"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
      cdf[idx,"LON"] <- compareNA(cdf[idx,"gamLON"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
      
      
      if(scen == "F0LMW")  { # won't match "overlap"
        cdf[idx,"LAT"] <-  ifelse(!is.na(cdf[idx,'gamLAT']), cdf[idx,"gamLAT"] >= 20 & cdf[idx,"gamLAT"] <= 25,FALSE)
        cdf[idx,"LON"] <-  ifelse(!is.na(cdf[idx,'gamLON']), cdf[idx,"gamLON"] >= 20 & cdf[idx,"gamLON"] <= 25,FALSE)
      }
      cdf[idx,"YEAR"] <- compareNA(cdf[idx,"gamYR"],scenarios$TEMPORAL[scenarios$DESC == cdf[idx,'scen']][1])
      
      ptf <- ptf.miss <- c(NA,NA)
      
      for(v in  1:length(parms)){ ## loop variables
        tmp.true <- subset(matchdf, nscen == scen & guessCREG == c & guessREG == n & variable == parms[v])
        tmp.par <- parest %>% filter(REG == n & cREG == c &  variable == parms[v] & source == 'Estimated')
        if(dim(tmp.par)[1]==0) next("NO DATA")
        bounds <- tmp.par %>% select(lwr,upr)
        ptf[v] <- tmp.true$value >= bounds[1] & tmp.true$value <= bounds[2]
        ptf.miss[v] <- missby(ruling = ptf[v], bounds = bounds, truth = tmp.true)
      } ## end parms
      cdf[idx,'REG'] <- tmp.true$trueREG
      cdf[idx,'gamREG'] <- paste(n,c)
      cdf[idx,parms[1]] <- ptf[1];  cdf[idx,parms[2]] <- ptf[2]
      cdf[idx,'L1_MISS'] <- ptf.miss[1]; cdf[idx,'L2_MISS'] <- ptf.miss[2]; idx <- idx + 1
      
    } ## end unique periods
  } ## end not tempvar
} ## end unique REG in dat




