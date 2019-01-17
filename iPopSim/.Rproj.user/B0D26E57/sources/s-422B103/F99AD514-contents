
plotComps <- function(out_file){
  ## identify boot directories for this F run
  boots  <- list.dirs(out_file, recursive = T)[grep('boot',list.dirs(out_file, recursive = T))]
  
  ## bind all rows
  dat <- list.files(boots, full.names = T)[grep('SAA',list.files(boots, full.names = T))] %>%
    map(read.table,sep=",",header=T) %>%
    reduce(bind_rows) 
  
  dim(dat)
  ## make summary dfs and plot 
  dat0 <- dat %>%
    sample_n(10000) 
  
  
  ## for age comps (true values only, aging error is not binnable)
  dat1 <- dat0 %>%
    group_by(Year, Age) %>% 
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    filter(Year %in% c(50:99)) %>%
    mutate(block = ifelse(Year < 75, 'early','late'),
           Flev = tolower(ifelse(block == 'early', strsplit(as.character(dat$FSIM[1])," ")[[1]][1], strsplit(as.character(dat$FSIM[1])," ")[[1]][2])))
  
  ## length comps, use rounded vals for binning (only for illustration)
  dat2 <- dat0 %>%
    mutate(fish_size = round(fish_size)) %>%
    group_by(Year, fish_size) %>% 
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>% 
    filter(Year %in% c(50:99)) %>%
    mutate(block = ifelse(Year < 75, 'early','late'),
           Flev = tolower(ifelse(block == 'early', strsplit(as.character(dat$FSIM[1])," ")[[1]][1], strsplit(as.character(dat$FSIM[1])," ")[[1]][2])))
  
  
  
  datN <- dat1 %>% group_by(Year) %>% summarise(totN = sum(n))
  
  ## extract interesting columns and save 
  dat0 %>% 
    mutate(Length_cm = fish_size/10) %>% 
    select(Year, Age, AgeE, FSIM, Length_cm) %>% 
    write.csv(., paste0(out_file,"/",tolower(paste0(unique(dat$FSIM))),"sim_ageComp.csv"),row.names = F)
  
  ## plot ages
  ggplot(dat1, aes(x = Age, y = freq, fill= Flev)) +
    theme_bw() +
    theme(panel.grid = element_blank())+
    geom_area(stat = 'identity', alpha = 0.5) +
    geom_line()+
    facet_wrap( ~ Year, ncol = 4) +
    scale_x_continuous(limits = c(0,max(dat1$Age))) +
    scale_y_continuous(breaks = seq(0,0.4,0.2),labels = seq(0,0.4,0.2)) +
    annotate("text", x = 0.9*max(dat1$Age), y = 0.35, label = paste0("n = ",datN$totN)) +
    labs(y = 'proportion', title = paste0('Age comps, Last 25 Years') )+
    scale_fill_brewer(palette = "Dark2", name = "F Level")
  
  ggsave(last_plot(), file = paste0(getwd(),"/plots/",tolower(paste0(unique(dat$FSIM))),"_ageComps.jpg"), dpi = 480, height = 10, width = 8, unit = 'in')
  
  ## plot lengths
  ggplot(dat2, aes(x = fish_size, y = freq, fill= Flev)) +
    theme_bw() +
    theme(panel.grid = element_blank())+
    geom_area(stat = 'identity', alpha = 0.5) +
    geom_line()+
    facet_wrap( ~ Year, ncol = 4) +
    scale_x_continuous(limits = c(50,250)) +
    scale_y_continuous(breaks = seq(0,0.25,0.1),labels = seq(0,0.25,0.1)) +
    annotate("text", x = 0.9*250, y = 0.1, label = paste0("n = ",datN$totN)) +
    labs(y = 'proportion', x= "Length (mm)", title = paste0('Length comps, Last 25 Years') )+
    scale_fill_brewer(palette = "Dark2", name = "F Level")
  
  ggsave(last_plot(), file = paste0(getwd(),"/plots/",tolower(paste0(unique(dat$FSIM))),"_lenComps.jpg"), dpi = 480, height = 10, width = 8, unit = 'in')
}
