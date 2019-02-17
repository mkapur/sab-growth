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
    # scen <- paste(scen0,scen1,sep = "_")
    scen <- scen0
    outdir0 <-  paste0("C:/users/",compname,"/dropbox/uw/sab-growth/gam/plots/", scen)
    if(!exists(outdir0)) dir.create(outdir0)

    for(b in 1:nboot){ ## loop boots
      dat <- read.csv(paste0("C:/users/",compname,"/dropbox/uw/sab-growth/ipopsim/gendata/",
                             scen,"_",b,'.csv'))  %>% filter(Age == age)
      cat(scen,' boot ',b,' Age ',age,'\n')
      
  
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
  write.csv(ldf, paste0(getwd(),"/summary_tables/ldf_raw_a",age,".csv"),row.names = F)
  
  # return(ldf)
# }
  # ldf <- read.csv("summary_tables/ldf_raw_a6.csv")
  
  ## tabulate restuls ----
  
  ## assign levels so it plots in order
  # corlevs <- data.frame(scen = c("F0L1S_25_MED_MED_MED","F0L1S_30_MED_MED_MED",
  #                    "F0L1S_49_MED_MED_MED","F0LMW_MED_MED_MED","NoBreaks_MED_MED_MED"),
  #                  scen2 = c("Break at 25?","Break at 30?",
  #                    "Break at 49?",  "Overlap 20? - 25?", "No Breaks"))
  # 
  
  levels(ldf$scen) <-  c( "Break at 25 deg.", "Break at 49 deg.",
                         "Low Contrast at 25 deg.", "Overlap 20-25 deg.","No Breaks")
  
  # levels(ldf$scen) <- paste('Break At ',sub('.*\\_', "\\1", levels(ldf$scen)))
  ldf$lat_breaks2 <- factor(ldf$lat_breaks, levels=c(paste(1:50),NA))
  

  
  # trueb <- c(factor(NA),25,30,NA,49) ## true breaks, FILL OVERLAP BELOW
  trueb <- c(25,49,25,NA,NA)
  ntrue <- ldf %>% group_by(scen) %>% summarise(ct = n()) %>% data.frame() %>% mutate(trueb)
  ntrue <- rbind(ntrue, data.frame(scen = rep(levels(ldf$scen)[4],6),ct = rep(NA,6), trueb = 20:25))
  
  write.csv(ntrue, paste0(getwd(),"/summary_tables/ntrue_a6.csv"),row.names = F)
  
   ldf %>%
    group_by(scen,  lat_breaks2) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    ggplot(., aes(x = lat_breaks2, y = freq)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme(panel.grid = element_blank())+
     # scale_y_continuous(limits = c(0,max(ldfprop$freq))) +
    scale_x_discrete(limits = c(paste(1:55),NA),breaks = c(paste(seq(0,50,5)),NA)) +
    geom_vline(data = ntrue, aes(xintercept =trueb), col = 'red', linetype = 'dashed') +
    facet_wrap(~scen,ncol = 1) +
    labs(main = 'breaks identified', y = 'frequency', x = 'break location (deg. latitude)', main = 'spatial breaks')

  ggsave(last_plot(), file = "./plots/ldf_a6.png")
  
  ## get max values
  ldfprop <- read.csv("summary_tables/ldf_raw_a6.csv") %>%
    group_by(scen,  lat_breaks) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    group_by(scen) %>%
    filter(freq == max(freq)) %>% data.frame()
  
  write.csv(ldfprop, paste0(getwd(),"/summary_tables/ldf_prop_a6.csv"),row.names = F)
  
  cat('plotted & saved break tabulations \n')