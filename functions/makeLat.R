
makeLat<-function(dat){
  for(i in 1:nrow(dat)){  
    if(is.na(sptl)){
      dat$Latitude_dd[i] <- runif(1, 0.0, 50.0); dat$REG[i] <- 'R1' ## uniform range
      dat$Longitude_dd[i] <- runif(1, 0.0, 50.0); dat$REG[i] <- 'R1' ## uniform range
    } else if(sptl == 25){ ## 25 is shorthand for single uniform break
      dat$Latitude_dd[i] <-  ifelse(dat$REG[i] == 'R1', runif(1,0.0,25.0), runif(1,25.00001,50.0))
      dat$Longitude_dd[i] <-  ifelse(dat$REG[i] == 'R1', runif(1,0.0,25.0), runif(1,25.00001,50.0))
    } else if(sptl == 48){ 
      dat$Latitude_dd[i] <- runif(1, 0.0, 50.0)#; dat$REG[i] <- 'R1' ## uniform range
      # dat$Latitude_dd[i] <-  ifelse(dat$REG[i] == 'R1', runif(1,0.0,49.0), runif(1,49.00001,50.0))
      dat$Longitude_dd[i] <-  ifelse(dat$REG[i] == 'R1', runif(1,0.0,47.9999), runif(1,48.0,50.0))
    } else if(sptl == "overlap"){ ## 20 is shorthand for overlapping zones
      dat$Latitude_dd[i] <-  ifelse(dat$REG[i] == 'R1', runif(1,0.0,25.0), runif(1,20.0,50.0))
      dat$Longitude_dd[i] <-  ifelse(dat$REG[i] == 'R1', runif(1,0.0,25.0), runif(1,20.0,50.0))
    } ## end else
  } ## end nrow
  return(dat)
}

