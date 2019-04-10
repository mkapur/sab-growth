## compareBreaks
## function to record for each gamboot whether the right breakpoints were detected AND whether the coverage includes the truth

vars <- paste0('R',1:4,"_");   vis <- c("pooled",'early','late'); Rlevs <- c("R1", "R2","R3","R4",apply(expand.grid(vars, vis), 1, paste, collapse=""))
  compareNA <- function(v1,v2) {
    same <- (v1 == v2) | (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
  }

  parest <-  rep0 %>% 
    filter(variable %in% parms) %>% 
    mutate(source = 'Estimated')
  ## only extract regions from 'actual' that were present in original
  if(scen == 'NoBreaks'){
    ## for no spatial var only compare with R1
    parest <- rbind(parest, (read.csv("./input_data/true_ibm_vals.csv") %>% 
                               filter(REG == "R1" & variable %in% parms))) 
  }else if(scen == "F0L1S_R3"){
    parest <- rbind(parest, (read.csv("./input_data/true_ibm_vals.csv") %>% 
                               filter(REG %in% c("R1","R3") & variable %in% parms))) 
  } else{
    parest <- rbind(parest, (read.csv("./input_data/true_ibm_vals.csv") %>% 
                               filter(REG %in% c("R1","R2") & variable %in% parms))) 
  }
  
  parest$REG <- factor(parest$REG, levels=Rlevs)
  parest$lwr <- parest$value - 1.96*parest$sd
  parest$upr <- parest$value + 1.96*parest$sd
  
   if(scen != 'tempvar_R1R2'){
    gamREGS <-  factor(unique(dat$gamREG), levels=Rlevs)
 } else {
    gamREGS <- factor(unique(dat$cREG),levels=Rlevs)
  }

  REGS <-  factor(unique(dat$REG),levels=Rlevs) ## what was used to GENERATE
  parms <- c('L1','L2','k','Linf')[1:2]
  if(is.na(breaksdf$yr_breaks)){
    for(iq in 1:length(gamREGS)){ ## loop regions, this leads to varying # samples (denoms)
      ptf <- NULL
      cdf[idx,'scen'] <- scen
      cdf[idx,'boot'] <- b
      cdf[idx,'REG'] <- as.factor(if(!is.na(REGS[iq])){REGS[iq]}else{last(REGS)}) ## if anything other than R1, will compare to appropriate RX (2 or 3)
      
      cdf[idx,'gamREG'] <- if(!is.na(gamREGS[iq])){gamREGS[iq]}else{last(gamREGS)}
      cdf[idx,'gamLAT'] <- breaksdf$lat_breaks
      cdf[idx,'gamLON'] <- breaksdf$lon_breaks
      cdf[idx,'gamYR'] <- breaksdf$yr_breaks
      
      cdf[idx,"LAT"] <- compareNA(cdf[idx,"gamLAT"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
      cdf[idx,"LON"] <- compareNA(cdf[idx,"gamLON"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
      cdf[idx,"YEAR"] <- compareNA(cdf[idx,"gamYR"],scenarios$TEMPORAL[scenarios$DESC == cdf[idx,'scen']][1])
      
      if(scen == "F0LMW")  { # won't match "overlap"
        cdf[idx,"LAT"] <- cdf[idx,"gamLAT"] >= 20 & cdf[idx,"gamLAT"] <= 25
        cdf[idx,"LON"] <- cdf[idx,"gamLON"] >= 20 & cdf[idx,"gamLON"] <= 25
      }
      
      for(v in  1:length(parms)){ ## loop variables
        tmp <- subset(parest, (variable == parms[v] & source == 'Estimated'  &  REG == gamREGS[iq]) | (variable == parms[v] & source == 'Actual' & REG == cdf[idx,'REG']))
        bounds <- tmp %>% filter(source == 'Estimated' & REG ==  cdf[idx,'gamREG'] ) %>% select(lwr,upr)
        ptf[v] <- tmp$value[tmp$source == 'Actual' &  tmp$REG ==    cdf[idx,'REG']] >= bounds[1] &
          tmp$value[tmp$source == 'Actual'&  tmp$REG ==    cdf[idx,'REG']] <= bounds[2]
      } ## end variables

      cdf[idx,parms[1]] <- ptf[1];  cdf[idx,parms[2]] <- ptf[2]; idx <- idx + 1
    } ## end no tvar
  # } else if(scen =="tempvar_R1R2")  {   ## compare EARLY to R1, LATE to R2
  }else{
        for(t in 1:length(c('early','late'))){
          tmp0 <- subset(parest, source == 'Estimated' & cREG ==  c('early','late')[t])
          for(j in 1:length(unique(tmp0$REG))){ ## loop rows in this period
            for(v in  1:length(parms)){ ## loop variables
              tmp <- rbind(subset(tmp0, REG == unique(tmp0$REG)[j] & variable == parms[v]), subset(parest, variable == parms[v] & source == 'Actual' & 
                                                                                                     REG == ifelse(t == 1  | scen == 'NoBreaks', 'R1', ifelse(scen == 'F0L1S_R3','R3','R2'))))
              bounds <- tmp %>% filter(source == 'Estimated') %>% select(lwr,upr)
              ptf[v] <- tmp$value[tmp$source == 'Actual'] >= bounds[1] &tmp$value[tmp$source == 'Actual'] <= bounds[2]
            } ## end parms
            cdf[idx,'scen'] <- scen
            cdf[idx,'boot'] <- b
            cdf[idx,'REG'] <- as.factor(if(!is.na(REGS[iq])){REGS[iq]}else{last(REGS)})
            
            cdf[idx,'gamREG'] <- paste0(tmp$REG[tmp$source=='Estimated'],"_",tmp$cREG[tmp$source=='Estimated'])
            cdf[idx,'gamLAT'] <- breaksdf$lat_breaks
            cdf[idx,'gamLON'] <- breaksdf$lon_breaks
            cdf[idx,'gamYR'] <- breaksdf$yr_breaks
            
            cdf[idx,"LAT"] <- compareNA(cdf[idx,"gamLAT"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
            cdf[idx,"LON"] <- compareNA(cdf[idx,"gamLON"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
            cdf[idx,"YEAR"] <- compareNA(cdf[idx,"gamYR"],scenarios$TEMPORAL[scenarios$DESC == cdf[idx,'scen']][1])
            
            if(scen == "F0LMW")  { # won't match "overlap"
              cdf[idx,"LAT"] <- cdf[idx,"gamLAT"] >= 20 & cdf[idx,"gamLAT"] <= 25
              cdf[idx,"LON"] <- cdf[idx,"gamLON"] >= 20 & cdf[idx,"gamLON"] <= 25
            }
            
            cdf[idx,parms[1]] <- ptf[1];  cdf[idx,parms[2]] <- ptf[2];idx <- idx+1
          } ## end j rows
        } ## end periods
    } ## end if R1R2
      

 
      
    # } ## end if not R1R2
    


  #   }
  #   if(scen =="NoBreaks")  {   ## compare EARLY to R1, LATE to R2
  #     for(j in 1:length(unique(parest$REG))){ ## loop REGs
  #       tmp0 <- subset(parest, source == 'Estimated' & REG ==  unique(parest$REG)[j] )
  #       cdf[idx,'scen'] <- scen
  #       ptf <- NULL
  #       cdf[idx,'boot'] <- b
  #       ## qualitative compare for regional
  #       cdf[idx,'gamREG'] <- tmp0$REG[1]
  #       cdf[idx,'REG'] <- 'R1'
  #       
  #       ## direct comp for break locations
  #       cdf[idx,'gamLAT'] <- breaksdf$lat_breaks
  #       cdf[idx,'gamLON'] <- breaksdf$lon_breaks
  #       cdf[idx,'gamYR'] <- breaksdf$yr_breaks
  #       
  #       cdf[idx,"LAT"] <- compareNA(cdf[idx,"gamLAT"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
  #       cdf[idx,"LON"] <- compareNA(cdf[idx,"gamLON"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
  #       
  #       for(v in  1:length(parms)){ ## loop variables
  #         tmp <- rbind(subset(tmp0, REG == unique(parest$REG)[j] & variable == parms[v]), subset(parest, variable == parms[v] & source == 'Actual'))
  #         bounds <- tmp %>% filter(source == 'Estimated') %>% select(lwr,upr)
  #         ptf[v] <- tmp$value[tmp$source == 'Actual'] >= bounds[1] &tmp$value[tmp$source == 'Actual'] <= bounds[2]
  #       } ## end parms
  #     } ## end rows in tmp0
  #   } ## end if nobreaks
  #   
  #       
  #         cdf[idx,'scen'] <- scen
  #         ptf <- NULL
  #         cdf[idx,'boot'] <- b
  #         ## qualitative compare for regional
  #         cdf[idx,'gamREG'] <- if(!is.na(gamREGS[iq])){gamREGS[iq]}else{last(gamREGS)}
  #         cdf[idx,'REG'] <- if(!is.na(REGS[iq])){REGS[iq]}else{last(REGS)}
  #         
  #         ## direct comp for break locations
  #         cdf[idx,'gamLAT'] <- breaksdf$lat_breaks
  #         cdf[idx,'gamLON'] <- breaksdf$lon_breaks
  #         cdf[idx,'gamYR'] <- breaksdf$yr_breaks
  #         
  #         cdf[idx,"LAT"] <- compareNA(cdf[idx,"gamLAT"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
  #         cdf[idx,"LON"] <- compareNA(cdf[idx,"gamLON"],scenarios$SPATIAL[scenarios$DESC == cdf[idx,'scen']][1])
  #         cdf[idx,"YEAR"] <- compareNA(cdf[idx,"gamYR"],scenarios$TEMPORAL[scenarios$DESC == cdf[idx,'scen']][1])
  #         

  #         } ## end parms
  #     
  #       } ## end rows in tmp0
  #     } ## end time periods
  #   } ## end if R1R2
  #   # if(is.na(cdf[idx,"LAT"])) stop(paste0("Throwing NAs in spatial comparison ", scen, " ", idx, "\n"))
  #   
  #   idx <- idx+1
  # } ## end gamREGS
  # idx <- idx+1
  # cdf0 <- cdf %>% filter(!is.na(scen)) 
  # return(list(cdf0,idx))
# }
