getGR <- function(tempdf,breaksdf){
  blat <- breaksdf$lat_breaks2
  blon <- breaksdf$lon_breaks2
  byr <- breaksdf$yr_breaks
  for(i in 1:nrow(tempdf)){ ## loop unique breaks
    tempdf$gamREG[i] <- ifelse(tempdf[i,"Latitude_dd"] >= blat & 
                                 tempdf[i,"Longitude_dd"] >= blon, "R3", 
                               ifelse(tempdf[i,"Latitude_dd"] >= blat & 
                                        tempdf[i,"Longitude_dd"] < blon, "R2",
                                      ifelse(tempdf[i,"Latitude_dd"] < blat & 
                                               tempdf[i,"Longitude_dd"] < blon, "R1","R4")))
    if (is.na(blat)  & is.na(blon)){tempdf$gamREG[i] <- 'R1'; next()}
    if (is.na(blat)) {
      if (tempdf[i, "Longitude_dd"] >= blon &  is.na(blat)) { tempdf$gamREG[i] <- 'R3' }
      if (tempdf[i, "Longitude_dd"] < blon &is.na(blat)) { tempdf$gamREG[i] <- 'R1'}
    } else if (is.na(blon)) {
      if (tempdf[i, "Latitude_dd"] >= blat &  is.na(blon)) { tempdf$gamREG[i] <- 'R3' }
      if (tempdf[i, "Latitude_dd"] < blat &  is.na(blon)) {tempdf$gamREG[i] <- 'R1' }
    }
    tempdf[i,'Period'] <- ifelse(tempdf[i,'Year'] < byr, 'early','late')
    
  }
  return(tempdf)
}
## sanity check
# tempdf %>% group_by(gamREG) %>% dplyr::summarise(mean(Longitude_dd))
# 