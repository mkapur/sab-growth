rm(list = ls())
setwd("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/gam")
source("./makeMod.R");source("./getBreaks.R")
fLevs <- read.csv('C:/Users/Maia Kapur/Dropbox/UW/sab-growth/iPopSim/inputs/scenarios.csv',na.strings = 'NA') ## manual file
nboot <- 25
require(mgcv);require(dplyr)

breaks_df <- ydf <-  ldf <- data.frame()

for(l in c(1:6,35,36)){
  
  ## get scenario name
  scen0 <- paste0(fLevs[l,'DESC'])
  scen1 <-  ifelse(is.na(fLevs[l, 5]),
                   paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
                   paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))
  scen <- paste(scen0,scen1,sep = "_")
  
  for(b in 1:nboot){ ## loop boots
    dat <- read.csv(paste0("C:/users/maia kapur/dropbox/uw/sab-growth/ipopsim/gendata/",
                           scen,"_",b,'.csv')) %>% filter(Age == 4)
    
    outdir0 <-  paste0("C:/users/maia kapur/dropbox/uw/sab-growth/gam/plots/", scen)
    if(!exists(outdir0)) dir.create(outdir0)
    outdir <- paste0(outdir0,"/boot_",b); if(!exists(outdir)) dir.create(outdir)
    
    png(paste0(outdir,"/rawData.png"), width = 10, height = 7, units = 'in', res = 420)
    par(mfrow = c(1,3))
    with(dat, hist(Length_cm, main = 'Age 4 Length_cm'))
    with(dat, plot(Length_cm ~ Latitude_dd, main = 'Age 4 Length_cm vs Lat'))
    with(dat, plot(Length_cm ~ Year, main = 'Age 4 Length_cm vs Year'))
    graphics.off()
    
    mod <- makeMod(scenario = scen,dat)
    bdf <-  getBreaks(mod,dat, scen)
    
    ## fill NAs in bdf for binding
    if (length(bdf[[1]]) == 0){ ## fill NA for empty
      bdf[[1]] <- NA
    } else if (length(bdf[[2]]) == 0){
      bdf[[2]] <- NA
    }
    tydf <- cbind(as.numeric(bdf[[1]]), as.character(rep(scen, length(bdf[[1]]))), rep(b, length(bdf[[1]])))
    tldf <- cbind(as.numeric(bdf[[2]]), as.character(rep(scen, length(bdf[[2]]))), rep(b, length(bdf[[2]])))
   ydf <- rbind(ydf,tydf); ldf <- rbind(ldf, tldf)
    
  } ## end boots

} ## end fLevs

names(ydf) <- c('year_breaks','scen','boot');   names(ldf) <- c('lat_breaks','scen','boot')

ldf$lat_breaks2 <- factor(ldf$lat_breaks, levels=c(paste(1:50),NA))
ldf %>% 
  group_by(scen,  lat_breaks2) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>%
  ggplot(., aes(x = lat_breaks2, y = freq)) +
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(limits = c(paste(1:50),NA),breaks = c(paste(seq(1,50,4)),NA)) +
  facet_wrap(~scen) + 
  labs(main = 'breaks identified', x = 'break location (latitude)', main = 'spatial breaks')

ldf %>% 
  group_by(scen,  lat_breaks) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n)) %>%
  write.csv(paste0(getwd(),"/summary_tables/ldf_prop.csv"),row.names = F)

ydf$year_breaks2 <- factor(ydf$year_breaks, levels=c(paste(1:50),NA))

ydf %>% 
  group_by(scen,  year_breaks2) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n)) %>% 
  ggplot(., aes(x = year_breaks2, y = freq)) +
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(limits = c(paste(1:50),NA),breaks = c(paste(seq(1,50,4)),NA)) +
  facet_wrap(~scen) + 
  labs(main = 'breaks identified', x = 'break location (year)', main = 'temporal breaks')


ydf %>% 
  group_by(scen,  year_breaks) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n)) %>%
  write.csv(paste0(getwd(),"/summary_tables/ydf_prop.csv"),row.names = F)
