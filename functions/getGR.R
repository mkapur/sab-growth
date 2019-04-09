getGR <- function(tempdf,breaksdf){
  blat <- breaksdf$lat_breaks2
  blon <- breaksdf$lon_breaks2
  byr <- breaksdf$yr_breaks
  if(length(blat) == 1 & length(blon) == 1 & length(byr)==1){  ## assign regions for comparison with endpoints
    for(i in 1:nrow(tempdf)){ ## loop unique breaks
    tempdf[i,'Period'] <- ifelse(!is.na(byr),ifelse(tempdf[i,'Year'] < byr[1], 'early','late'),'pooled')
      
      tempdf$gamREG[i] <- ifelse(tempdf[i,"Latitude_dd"] >= blat &
                                   tempdf[i,"Longitude_dd"] >= blon, "R3",
                                 ifelse(tempdf[i,"Latitude_dd"] >= blat &
                                          tempdf[i,"Longitude_dd"] < blon, "R2",
                                        ifelse(tempdf[i,"Latitude_dd"] < blat &
                                                 tempdf[i,"Longitude_dd"] < blon, "R1","R4")))
      
      if (is.na(blat)  & is.na(blon) & tempdf[i,'Period'] == 'early'){tempdf$gamREG[i] <- 'R1'; next()}
      if (is.na(blat)  & is.na(blon) & tempdf[i,'Period'] == 'late'){tempdf$gamREG[i] <- 'R2'; next()}
      if (is.na(blat)  & is.na(blon)){tempdf$gamREG[i] <- 'R1'; next()}

      if (is.na(blat)) {
        if (tempdf[i, "Longitude_dd"] >= blon &  is.na(blat)) { tempdf$gamREG[i] <- 'R3' }
        if (tempdf[i, "Longitude_dd"] < blon &is.na(blat)) { tempdf$gamREG[i] <- 'R1'}
      } else if (is.na(blon)) {
        if (tempdf[i, "Latitude_dd"] >= blat &  is.na(blon)) { tempdf$gamREG[i] <- 'R3' }
        if (tempdf[i, "Latitude_dd"] < blat &  is.na(blon)) {tempdf$gamREG[i] <- 'R1' }
      }
      
    } ## end loop
  } ## end if length of all breaks == 1
  else{
    for(i in 1:nrow(tempdf)){ ## loop unique breaks -- we are working with two of each due to time breaks (SAB ONLY)
      if(tempdf[i,"Latitude_dd"] >= blat[2]){ ## northernmost set
        if(tempdf[i,"Longitude_dd"] <= blon[2]){
          tempdf$gamREG[i] <- "R5"
        } else if(tempdf[i,"Longitude_dd"] > blon[2] & tempdf[i,"Longitude_dd"] <= blon[1]){
          tempdf$gamREG[i] <- "R4"
        } else if(tempdf[i,"Longitude_dd"] > blon[2] & tempdf[i,"Longitude_dd"] > blon[1]){
          tempdf$gamREG[i] <- "R3"
        }
      } else if(tempdf[i,"Latitude_dd"] > blat[1] & tempdf[i,"Longitude_dd"] > blon[1]){
        tempdf$gamREG[i] <- "R2"
      } else{
        tempdf$gamREG[i] <- "R1"
      }
      tempdf[i,'Period'] <- ifelse(tempdf[i,'Year'] < byr[1], 'early','late')
      
    } ## end loop
    
    
    
  } ## end else
  
  
  return(tempdf)
}
## sanity check
# tempdf %>% group_by(gamREG) %>% dplyr::summarise(mean(Longitude_dd))
# 