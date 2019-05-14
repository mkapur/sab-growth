require(mgcv)
## wrapper function for getBreaks and makeMod, given Age
# source("./makeMod.R")
source("./functions/getBreaks.R")

ldf <- data.frame()
  
  ## loop through boots, build GAMS and ID Breakpoints ----
  for(l in 1:length(testrows)){
    
    scen0 <- testrows[l]
    scen <- scen0
    outdir0 <-  paste0("C:/users/",compname,"/dropbox/uw/sab-growth/gam/plots/", 
                       scen)
    if(!exists(outdir0)) dir.create(outdir0)
    
    for(b in 1:nboot){ ## loop boots
      dat <- read.csv(paste0("./IBM_output/datasets/",scen,"_",b,'.csv'))  %>% filter(Age == age)
      cat( testrows[l],' boot ',b,' Age ',age,'\n')

      mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd), data = dat) 
      save(mod, file = paste0("./GAM_output/mod_",scen,"_",b,".rds")) ## save on the fly
      bdf <-  getBreaks(mod,dat, scen)
      
      ## fill NAs in bdf for binding
      for(i in 1:length(bdf)){
        if (length(bdf[[i]]) == 0){ ## fill NA for empty
          bdf[[i]] <- NA
        }
      }
  
      tldf <- cbind(as.numeric(bdf[[1]]), as.numeric(bdf[[2]]), as.numeric(bdf[[3]]), as.character(rep(scen, length(bdf[[1]]))), rep(b, length(bdf[[1]])))
      
      ldf <- rbind(ldf, tldf)
 
    } ## end boots
    
  } ## end fLevsd
  ## save
  names(ldf) <- c('yr_breaks','lat_breaks','lon_breaks','scen','boot')
  ldf$lat_breaks2 <- factor(ldf$lat_breaks, levels=c(paste(1:50),NA))
  ldf$lon_breaks2 <- factor(ldf$lon_breaks, levels=c(paste(1:50),NA))
  write.csv(ldf, file = paste0("./GAM_output/ldf_raw_a",age,".csv"),row.names = F)
  cat('plotted & saved break tabulations in GAM_output \n')
