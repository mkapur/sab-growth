## wrapper function for getBreaks and makeMod, given Age

source("./makeMod.R");source("./getBreaks.R")


# bootBreaks <- function(Age = 6, nboot = 100, testrows =  c(1:nrow(fLevs))){

  ldf <- data.frame()
  
  ## loop through boots, build GAMS and id Breakpoints ----
  for(l in 1:length(testrows)){
    
    
    ## get scenario name
    # scen0 <- paste0(fLevs[l,'DESC'])
    # scen1 <-  ifelse(is.na(fLevs[l, 5]),
    #                  paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
    #                  paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))
    scen0 <- testrows[l]; scen1 <- paste0(rep('MED',3),collapse = "_")
    scen <- paste(scen0,scen1,sep = "_")


    for(b in 1:nboot){ ## loop boots
      dat <- read.csv(paste0("C:/users/",compname,"/dropbox/uw/sab-growth/ipopsim/gendata/",
                             scen,"_",b,'.csv'))  %>% filter(Age == age)
      cat(scen,' boot ',b,' Age ',age,'\n')
      
      outdir0 <-  paste0("C:/users/",compname,"/dropbox/uw/sab-growth/gam/plots/", scen)
      if(!exists(outdir0)) dir.create(outdir0)
      outdir <- paste0(outdir0,"/boot_",b); if(!exists(outdir)) dir.create(outdir)
      
      png(paste0(outdir,"/rawData.png"), width = 10, height = 7, units = 'in', res = 420)
      par(mfrow = c(1,3))
      with(dat, hist(Length_cm, main = paste0('Age ',age,' Length_cm')))
      with(dat, plot(Length_cm ~ Latitude_dd, main = paste0('Age ',age,' Length_cm vs Lat')))
      with(dat, plot(Length_cm ~ Year, main = paste0('Age ',age,' Length_cm vs Year')))
      graphics.off()
      
      mod <- makeMod(scenario = scen,dat)
      bdf <-  getBreaks(mod,dat, scen)
      
      ## fill NAs in bdf for binding
      if (length(bdf[[1]]) == 0){ ## fill NA for empty
        bdf[[1]] <- NA
      }
      # } else if (length(bdf[[2]]) == 0){
      #   bdf[[2]] <- NA
      # }
      tldf <- cbind(as.numeric(bdf[[1]]), as.character(rep(scen, length(bdf[[1]]))), rep(b, length(bdf[[1]])))
      # tldf <- cbind(as.numeric(bdf[[2]]), as.character(rep(scen, length(bdf[[2]]))), rep(b, length(bdf[[2]])))
      # ydf <- rbind(ydf,tydf);
      ldf <- rbind(ldf, tldf)
      
    } ## end boots
    
  } ## end fLevs
  
  ## save
  names(ldf) <- c('lat_breaks','scen','boot')
  ldf$lat_breaks2 <- factor(ldf$lat_breaks, levels=c(paste(1:50),NA))
  write.csv(ldf, paste0(getwd(),"/summary_tables/ldf_raw_age_",age,".csv"),row.names = F)
  
  # return(ldf)
# }