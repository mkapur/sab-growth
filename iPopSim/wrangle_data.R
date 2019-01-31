build_simComp <- function(out_file, l=l, dat0 = NA){
  
  if(is.na(dat0)){
  ## identify boot directories for this F run
  boots  <- list.dirs(out_file, recursive = T)[grep('boot',list.dirs(out_file, recursive = T))]
  
  ## bind all rows
  dat <- list.files(boots, full.names = T)[grep('SAA',list.files(boots, full.names = T))] %>%
    lapply(read.table,sep=",",header=T) %>%
    reduce(bind_rows) 
  
  ## make summary dfs and plot 
  dat0 <- dat %>%
    sample_n(10000) 
  scenID <-
    ifelse(is.na(fLevs[l, 5]),
           paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
           paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))  ## extract interesting columns and save 
  dat0 %>% 
    mutate(Length_cm = fish_size/10) %>% 
    select(Year, Age, AgeE, Length_cm, FSIM, FMORT, REG) %>% 
    write.csv(., paste0(out_file,"/",scenID,"sim_Comp.csv"),row.names = F)

  }
  
  ## for age comps (true values only, aging error is not binnable)
  scenID <-
    ifelse(is.na(fLevs[l, 5]),
           paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
           paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))  
  
  dat1 <- dat0 %>%
    group_by(Year, Age) %>% 
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    filter(Year %in% c(50:99))
  
  dat1$Flev <- NA
  
  
  
  if(is.na(fLevs[l,5])){ ## if in 25 year blocks fill both
    dat1$Flev[dat1$Year %in% 0:49] <- "ZERO"
    dat1$Flev[dat1$Year %in% 50:75] <- as.character(fLevs[l,'FMORT_1'])
    dat1$Flev[dat1$Year %in% 76:101]<- as.character(fLevs[l,'FMORT_2'])
  } else {
    dat1$Flev[dat1$Year %in% 0:49] <- "ZERO"
    dat1$Flev[dat1$Year %in% 50:67] <- as.character(fLevs[l,'FMORT_1'])
    dat1$Flev[dat1$Year %in% 68:85]<- as.character(fLevs[l,'FMORT_2'])
    dat1$Flev[dat1$Year %in% 86:101]<- as.character(fLevs[l,'FMORT_3'])
  }
  
  
  ## length comps, use rounded vals for binning (only for illustration)
  dat2 <- dat0 %>%
    mutate(fish_size = round(fish_size)) %>%
    group_by(Year, fish_size) %>% 
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    filter(Year %in% c(50:99)) #%>%
    # mutate(block = ifelse(Year < 75, 'early','late'),
    #        Flev = tolower(ifelse(block == 'early', strsplit(as.character(dat$FSIM[1])," ")[[1]][1], strsplit(as.character(dat$FSIM[1])," ")[[1]][2])))
  dat2$Flev <- NA
  if(is.na(fLevs[l,5])){ ## if in 25 year blocks fill both
    dat2$Flev[dat2$Year %in% 0:49] <- "ZERO"
    dat2$Flev[dat2$Year %in% 50:75] <- as.character(fLevs[l,'FMORT_1'])
    dat2$Flev[dat2$Year %in% 76:101]<- as.character(fLevs[l,'FMORT_2'])
  } else {
    dat2$Flev[dat2$Year %in% 0:49] <- "ZERO"
    dat2$Flev[dat2$Year %in% 50:67] <- as.character(fLevs[l,'FMORT_1'])
    dat2$Flev[dat2$Year %in% 68:85]<- as.character(fLevs[l,'FMORT_2'])
    dat2$Flev[dat2$Year %in% 86:101]<- as.character(fLevs[l,'FMORT_3'])
  }
  
  
  datN <- dat1 %>% group_by(Year) %>% summarise(totN = sum(n))
  

  return(list(dat1,dat2,datN,scenID))
  
}


plotComps <- function(dat1,dat2,datN,scenID, saveloc = NA){
 

  ## plot ages
 p <- ggplot(dat1, aes(x = Age, y = freq, fill= Flev)) +
    theme_bw() +
    theme(panel.grid = element_blank())+
    geom_area(stat = 'identity', alpha = 0.5) +
    geom_line()+
    facet_wrap( ~ Year, ncol = 4) +
    scale_x_continuous(limits = c(0,max(dat1$Age))) +
    scale_y_continuous(breaks = seq(0,0.4,0.2),labels = seq(0,0.4,0.2)) +
    annotate("text", x = 0.9*max(dat1$Age), y = 0.35, label = paste0("n = ",datN$totN)) +
    labs(y = 'proportion', title = paste0("Age comps") )+
    scale_fill_brewer(palette = "Dark2", name = "F Level")
  
  if(is.na(saveloc)){
    ggsave(p, file = paste0(getwd(),"/plots/",paste0(scenID,"_ageComps.jpg")), dpi = 480, height = 10, width = 8, unit = 'in')
  } else{
    ggsave(p, file = paste0(out_file,"/",paste0(scenID, "_ageComps.jpg")), dpi = 480, height = 10, width = 8, unit = 'in')
  }
  
  ## plot lengths
 p <-  ggplot(dat2, aes(x = fish_size*10, y = freq, fill= Flev)) +
    theme_bw() +
    theme(panel.grid = element_blank())+
    geom_area(stat = 'identity', alpha = 0.5) +
    geom_line()+
    facet_wrap( ~ Year, ncol = 4) +
    # scale_x_continuous(limits = c(50,250)) +
    # scale_y_continuous(breaks = seq(0,0.25,0.1),labels = seq(0,0.25,0.1)) +
    annotate("text", x = max(dat2$fish_size)*9, y = 0.1, label = paste0("n = ",datN$totN)) +
    labs(y = 'proportion', x= "Length (mm)", title = paste0('Length comps'))+
    scale_fill_brewer(palette = "Dark2", name = "F Level")
  
  if(is.na(saveloc)){
    ggsave(p, file = paste0(getwd(),"/plots/",paste0(scenID,"_lenComps.jpg")), dpi = 480, height = 10, width = 8, unit = 'in')
  } else{
    ggsave(p, file = paste0(out_file,"/",paste0(scenID,"_lenComps.jpg")), dpi = 480, height = 10, width = 8, unit = 'in')
  } 
}
 
