build_simComp <- function(out_file, l=l, dat0 = NA){
  
  # if(is.na(dat0)){
  ## identify boot directories for this F run
  # boots  <- list.dirs(out_file, recursive = T)[grep('boot',list.dirs(out_file, recursive = T))]
  
  ## bind all rows
  # dat <- list.files(boots, full.names = T)[grep('SAA',list.files(boots, full.names = T))] %>%
  #   lapply(read.table,sep=",",header=T) %>%
  #   reduce(bind_rows)
  # 
  # ## make summary dfs and plot 
  # dat0 <- dat %>%
  #   sample_n(10000)
  
  
  ## just read boot 1 cause issues with stoke
  # dat0 <- read.table(paste0(out_file,"/boot_1/IBM_SAA_MATRIX.txt"),header = T,sep = ",")
  dat1 <- dat2 <- datN <- list()
  
  for(b in 1:nboot){
    dat0 <- read.table(paste0(out_file,"/boot_",b,"/IBM_SAA_MATRIX.txt"),sep=",",header=T)
    # scenID <-
    #   ifelse(is.na(fLevs[l, 5]),
    #          paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
    #          paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))  ## extract interesting columns and save 
    dat0 %>% 
      mutate(Length_cm = fish_size) %>% 
      select(Year, Age, AgeE, Length_cm, FMORT, REG) %>% 
      write.csv(., paste0(out_file,"/boot_",b,"/",b,"_sim_Comp.csv"),row.names = F)
    
    # }
    
    ## for age comps (true values only, aging error is not binnable)
    # scenID <-
    #   ifelse(is.na(fLevs[l, 5]),
    #          paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
    #          paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))  
    
    dat1[[b]] <- dat0 %>%
      group_by(Year, Age) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n)) %>% 
      filter(Year %in% c(50:99))
    
    dat1[[b]]$Flev <- 'No Fishing'
    
    # if(is.na(fLevs[l,5])){ ## if in 25 year blocks fill both
    #   dat1[[b]]$Flev[dat1[[b]]$Year %in% 0:49] <- "ZERO"
    #   dat1[[b]]$Flev[dat1[[b]]$Year %in% 50:75] <- as.character(fLevs[l,'FMORT_1'])
    #   dat1[[b]]$Flev[dat1[[b]]$Year %in% 76:101]<- as.character(fLevs[l,'FMORT_2'])
    # } else {
    #   dat1[[b]]$Flev[dat1[[b]]$Year %in% 0:49] <- "ZERO"
    #   dat1[[b]]$Flev[dat1[[b]]$Year %in% 50:67] <- as.character(fLevs[l,'FMORT_1'])
    #   dat1[[b]]$Flev[dat1[[b]]$Year %in% 68:85]<- as.character(fLevs[l,'FMORT_2'])
    #   dat1[[b]]$Flev[dat1[[b]]$Year %in% 86:101]<- as.character(fLevs[l,'FMORT_3'])
    # }
    
    
    ## length comps, use rounded vals for binning (only for illustration)
    dat2[[b]] <- dat0 %>%
      mutate(fish_size = round(fish_size)) %>%
      group_by(Year, fish_size) %>% 
      dplyr::summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% 
      filter(Year %in% c(50:99)) #%>%
    # mutate(block = ifelse(Year < 75, 'early','late'),
    #        Flev = tolower(ifelse(block == 'early', strsplit(as.character(dat$FSIM[1])," ")[[1]][1], strsplit(as.character(dat$FSIM[1])," ")[[1]][2])))
    dat2[[b]]$Flev <- 'No Fishing'
    # if(is.na(fLevs[l,5])){ ## if in 25 year blocks fill both
    #   dat2[[b]]$Flev[dat2[[b]]$Year %in% 0:49] <- "ZERO"
    #   dat2[[b]]$Flev[dat2[[b]]$Year %in% 50:75] <- as.character(fLevs[l,'FMORT_1'])
    #   dat2[[b]]$Flev[dat2[[b]]$Year %in% 76:101]<- as.character(fLevs[l,'FMORT_2'])
    # } else {
    #   dat2[[b]]$Flev[dat2[[b]]$Year %in% 0:49] <- "ZERO"
    #   dat2[[b]]$Flev[dat2[[b]]$Year %in% 50:67] <- as.character(fLevs[l,'FMORT_1'])
    #   dat2[[b]]$Flev[dat2[[b]]$Year %in% 68:85]<- as.character(fLevs[l,'FMORT_2'])
    #   dat2[[b]]$Flev[dat2[[b]]$Year %in% 86:101]<- as.character(fLevs[l,'FMORT_3'])
    # }
    # 
    
    datN[[b]] <- dat1[[b]] %>% group_by(Year) %>% summarise(totN = sum(n))
    
  } ## end boots
  
  return(list(dat1,dat2,datN))
  
}


plotComps <- function(dat1,dat2,datN, saveloc = NA, nboots = nboot){
  
  for(b in 1:nboot){
    
    ## plot ages
    p <- ggplot(dat1[[b]], aes(x = Age, y = freq, fill= Flev)) +
      theme_bw() +
      theme(panel.grid = element_blank())+
      geom_area(stat = 'identity', alpha = 0.5) +
      geom_line()+
      facet_wrap( ~ Year, ncol = 4) +
      # scale_x_continuous(limits = c(0,max(dat1[[b]]$Age))) +
      # scale_y_continuous(breaks = seq(0,0.4,0.2),labels = seq(0,0.4,0.2)) +
      # annotate("text", x = 9*max(dat1[[b]]$Age), y = 0.35, label = paste0("n = ",datN[[b]]$totN)) +
      labs(y = 'proportion', title = paste0("Age comps") )+
      scale_fill_brewer(palette = "Dark2", name = "F Level")
    
    # if(is.na(saveloc)){
      ggsave(p, file = paste0(out_file,"/boot_",b,"/",b,"_ageComps.jpg"), dpi = 480, height = 10, width = 8, unit = 'in')
    # } else{
      # ggsave(p, file = paste0(out_file,"/",paste0(b, "_ageComps.jpg")), dpi = 480, height = 10, width = 8, unit = 'in')
    # }
    
    ## plot lengths
    p <-  ggplot(dat2[[b]], aes(x = fish_size, y = freq, fill= Flev)) +
      theme_bw() +
      theme(panel.grid = element_blank())+
      geom_area(stat = 'identity', alpha = 0.5) +
      geom_line()+
      facet_wrap( ~ Year, ncol = 4) +
      scale_x_continuous(limits = c(50,250)) +
      # scale_y_continuous(breaks = seq(0,0.25,0.1),labels = seq(0,0.25,0.1)) +
      # annotate("text", x = max(dat2[[b]]$fish_size)*9, y = 0.1, label = paste0("n = ",datN[[b]]$totN)) +
      labs(y = 'proportion', x= "Length (cm)", title = paste0('Length comps'))+
      scale_fill_brewer(palette = "Dark2", name = "F Level")
    
    # if(is.na(saveloc)){
      ggsave(p, file = paste0(out_file,"/boot_",b,"/",b,"_lenComps.jpg"), dpi = 480, height = 10, width = 8, unit = 'in')
    # } else{
      # ggsave(p, file = paste0(out_file,"/",paste0(b, "_lenComps.jpg")), dpi = 480, height = 10, width = 8, unit = 'in')
    # }
    
  }## end nboot
  
} ## end function 


