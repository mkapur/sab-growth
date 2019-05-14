## compareBreaks
## function to record for each GAM or STARboot whether the right breakpoints were detected AND whether the coverage includes the truth

  compareNA <- function(v1,v2) {
    same <- (v1 == v2) | (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
  }
  parms <- c('L1','L2','k','Linf')[1:2]
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
  
  parest$REG <- factor(parest$REG[order(unique(parest$REG))], levels=Rlevs,ordered = TRUE)
  parest$lwr <- parest$value - 1.96*parest$sd
  parest$upr <- parest$value + 1.96*parest$sd
  
  if(scen != 'tempvar_R1R2'){
    gamREGS <-  factor(unique(dat$gamREG)[order(unique(dat$gamREG))], levels=Rlevs)
  } else {
    # gamREGS <- factor(unique(dat$cREG),levels=Rlevs)
    gamREGS <- factor(unique(dat$cREG)[order(unique(dat$gamREG))], levels= Rlevs, ordered=TRUE)
  }

  REGS <-  factor(unique(dat$REG)[order(unique(dat$REG))], levels= Rlevs, ordered=TRUE) ## what was used to GENERATE
  
  if(is.na(breaksdf$yr_breaks)  ){
    for(iq in 1:length(gamREGS)){ ## loop regions, this leads to varying # samples (denoms)
      
      t.gamREG <- factor(if(!is.na(gamREGS[iq])){gamREGS[iq]}else{last(gamREGS)},levels = Rlevs, ordered = T)
      t.REG <- factor(if(!is.na(REGS[iq])){REGS[iq]}else{last(REGS)}, levels = Rlevs, ordered = T)
        
      ptf <- NULL
      cdf[idx,'scen'] <- scen
      cdf[idx,'boot'] <- b
      cdf[idx,'REG'] <- t.REG ## if anything other than R1, will compare to appropriate RX (2 or 3)
      # cdf[idx,'REG'] <- ifelse(factor(gsub("_.*", "\\1", gamREGS[iq]),levels = Rlevs) %in% REGS, REGS[REGS == factor(gsub("_.*", "\\1", gamREGS[iq]),levels = Rlevs)],
  
      cdf[idx,'gamREG'] <- t.gamREG
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
        # if(){
          # tmp <- subset(parest, (variable == parms[v] & source == 'Estimated'  &  REG == factor(gsub("_.*", "\\1", gamREGS[iq]),levels = Rlevs)) | 
          #                 (variable == parms[v] & source == 'Actual' & REG == factor(cdf[idx,'REG'], levels = Rlevs)))
          tmp <- subset(parest, (variable == parms[v] & source == 'Estimated'  &  REG == gsub("_.*", "\\1", t.gamREG)) | 
                          (variable == parms[v] & source == 'Actual' & REG ==  ifelse(scen == 'tempvar_R1R2','R1',t.REG)))

          bounds <- tmp %>% filter(source == 'Estimated') %>% select(lwr,upr)
        # } else{
          # tmp <- subset(parest, (variable == parms[v] & source == 'Estimated'  &  REG == gamREGS[iq]) | 
          #                 (variable == parms[v] & source == 'Actual' & REG == factor(cdf[idx,'REG'], levels = Rlevs)))
          # bounds <- tmp %>% filter(source == 'Estimated' & REG ==  factor(cdf[idx,'gamREG'], levels = Rlevs)) %>% select(lwr,upr)
          
        # }
        # ptf[v] <- tmp$value[tmp$source == 'Actual' &  tmp$REG ==   factor(cdf[idx,'REG'], levels = Rlevs)] >= bounds[1] &
        #   tmp$value[tmp$source == 'Actual' &  tmp$REG ==     factor(cdf[idx,'REG'], levels = Rlevs)] <= bounds[2]
        ptf[v] <- tmp$value[tmp$source == 'Actual' ] >= bounds[1] &
          tmp$value[tmp$source == 'Actual' ] <= bounds[2]
      } ## end variables
      cdf[idx,parms[1]] <- ptf[1];  cdf[idx,parms[2]] <- ptf[2]; idx <- idx + 1
    } ## end iq
  } else{
    parest$ccREG <- ifelse(parest$source == 'Estimated', paste0(parest$REG,"_",parest$cREG),NA)
    unCCREG <- unique(parest$ccREG[!is.na(parest$ccREG)])
    for(iq in 1:length(unCCREG)){ ## loop regions, this leads to varying # samples (denoms) per boot
      tmp0 <- subset(parest, source == 'Estimated' & ccREG == unique( parest$ccREG)[iq])
      for(v in  1:length(parms)){ ## loop variables
        ## generate IFELSE based on tree. OK to do REG[1] because second is just L2
        if(!(scen %in% c('tempvar_R1R2','F0L1S_R3'))){
          tmp1 <- subset(parest, variable == parms[v] & source == 'Actual' & REG == ifelse(tmp0$REG[1] == 'R1' |scen == 'NoBreaks','R1','R2'))
        } else if(scen == 'F0L1S_R3'){
          tmp1 <- subset(parest, variable == parms[v] & source == 'Actual' & REG == ifelse(tmp0$REG[1] == 'R1','R1','R3'))
        } else if(scen == 'tempvar_R1R2'){
          tmp1 <- subset(parest, variable == parms[v] & source == 'Actual' & REG == ifelse(tmp0$cREG[1] == 'early','R1','R2'))
        }
        tmp <- rbind(subset(tmp0, variable == parms[v]), tmp1)
        bounds <- tmp %>% filter(source == 'Estimated') %>% select(lwr,upr)
        ptf[v] <- tmp$value[tmp$source == 'Actual'] >= bounds[1] &tmp$value[tmp$source == 'Actual'] <= bounds[2]
      } ## end parms
      # for(t in 1:length(c('early','late'))){
        #   parest$ccREG <- with(parest, paste0(REG,"_",cREG))
        #   # tmp0 <- subset(parest, source == 'Estimated' & ccREG ==  c('early','late')[t])
        #   tmp0 <- subset(parest, source == 'Estimated')
        # # for(j in 1:length(unique(tmp0$REG))){ ## loop REGs in this period
        #   for(v in  1:length(parms)){ ## loop variables
        #     tmp <- rbind(subset(tmp0, REG == unique(tmp0$REG)[j] & variable == parms[v]), subset(parest, variable == parms[v] & source == 'Actual' &
        #                                                                                            REG == ifelse(t == 1  | scen == 'NoBreaks', 'R1', ifelse(scen == 'F0L1S_R3','R3','R2'))))
        #     bounds <- tmp %>% filter(source == 'Estimated') %>% select(lwr,upr)
        #     ptf[v] <- tmp$value[tmp$source == 'Actual'] >= bounds[1] &tmp$value[tmp$source == 'Actual'] <= bounds[2]
        #   } ## end parms
          cdf[idx,'scen'] <- scen
          cdf[idx,'boot'] <- b
          cdf[idx,'REG'] <- unique(tmp1$REG)
            # as.factor(if(!is.na(REGS[iq])){REGS[iq]}else{last(REGS)})

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
          } ## end overlap
          cdf[idx,parms[1]] <- ptf[1];  cdf[idx,parms[2]] <- ptf[2];idx <- idx+1
          
        # } ## end j rows
      # } ## end periods
    } ## end ccREG
  } ## end else
  
  
  
  
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
