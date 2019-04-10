build_simComp <- function(out_file, l=l, dat0 = NA){
  
  dat1 <- dat2 <- datN <- list()
  
  for(b in 1:nboot.temp){
    dat0 <- read.table(paste0(out_file,"/boot_",b,"/IBM_SAA_MATRIX.txt"),sep=",",header=T)
 
    dat0 %>% 
      mutate(Length_cm = fish_size) %>% 
      select(Year, Age, AgeE, Length_cm, FMORT, REG) %>% 
      write.csv(., paste0(out_file,"/boot_",b,"/",b,"_sim_Comp.csv"),row.names = F)
    
 
    dat1[[b]] <- dat0 %>%
      group_by(Year, Age) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::mutate(freq = n / sum(n)) %>% 
      filter(Year %in% c(50:99))
    
    dat1[[b]]$Flev <- 'No Fishing'
    

    ## length comps, use rounded vals for binning (only for illustration)
    dat2[[b]] <- dat0 %>%
      mutate(fish_size = round(fish_size)) %>%
      group_by(Year, fish_size) %>% 
      dplyr::summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>% 
      filter(Year %in% c(50:99)) #%>%
   
    dat2[[b]]$Flev <- 'No Fishing'

    
    datN[[b]] <- dat1[[b]] %>% group_by(Year) %>% summarise(totN = sum(n))
    
  } ## end boots
  
  return(list(dat1,dat2,datN))
  
}


plotComps <- function(dat1,dat2,datN, saveloc = NA, nboot = nboot.temp){
  
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
    
      ggsave(p, file = paste0(out_file,"/boot_",b,"/",b,"_ageComps.jpg"), dpi = 480, height = 10, width = 8, unit = 'in')

    
    ## plot lengths
    p <-  ggplot(dat2[[b]], aes(x = fish_size, y = freq, fill= Flev)) +
      theme_bw() +
      theme(panel.grid = element_blank())+
      geom_area(stat = 'identity', alpha = 0.5) +
      geom_line()+
      facet_wrap( ~ Year, ncol = 4) +
      scale_x_continuous(limits = c(50,250)) +
      labs(y = 'proportion', x= "Length (cm)", title = paste0('Length comps'))+
      scale_fill_brewer(palette = "Dark2", name = "F Level")
          ggsave(p, file = paste0(out_file,"/boot_",b,"/",b,"_lenComps.jpg"), dpi = 480, height = 10, width = 8, unit = 'in')

    
  }## end nboot
  
} ## end function 



