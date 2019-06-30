## PlotMaster
## Figures for Manuscript/Chapter 1

require(dplyr)
require(ggplot2)
require(gridExtra)
require(RColorBrewer)
library(ggpubr)
require(DescTools)
require(mapdata)
require(mgcv)




## GAM dataviz  ----

## manual file
age <- 6

## how many NA years (should be most!!)
# ldf %>% group_by(scen) %>% dplyr::summarise(naY = sum(is.na(yr_breaks)), n = n()) %>% mutate(naY/n)

# "g) 20% Higher k, Break at 25 deg.""c) Scenario 3 20% Higher k, Break at 25 deg.",
## panel plot of one gendata for each scenario ---
source("./functions/Deriv.R")
scenarios <- read.csv('./input_data/scenarios.csv',na.strings = 'NA') %>% filter(DESC != 'F0L1S_R3')
scens <-  unique(scenarios$DESC)
scens.title <- scens.title0 <- scens
levels(scens.title) <-  c("g) Scenario 2 Break at 25 deg.", "g) Scenario 4 Break at 48 deg.", 
                          "g) Scenario 3 Overlap 20-25 deg.","g) Scenario 1 No Breaks",
                          "g) Scenario 5 Temporal Break Year 50")
levels(scens.title0) <-  c("b) Scenario 2 Break at 25 deg.", "d) Scenario 4 Break at 48 deg.", 
                           "c) Scenario 3 Overlap 20-25 deg.","a) Scenario 1 No Breaks",
                          "e) Scenario 5 Temporal Break Year 50")

# cdf_gam %>% filter(LAT == TRUE & LON == TRUE & YEAR == TRUE & scen == 'F0L1S_49')
bootpicks <-c(5,33,3,83,12,2)
plist0 <- list(); idx0 <- 1
for(i in 1:length(scens)){
  plist  <- list(); idx <- 1
  ## plot raw data
  scen <- scens[i]
  tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_",bootpicks[i],".csv")) %>% filter(Age == 6)

  breaksdf <- read.csv( paste0("./GAM_output/ldf_raw_a6.csv")) %>% filter(scen == scens[i] & boot == bootpicks[i])
  load(file = paste0("./GAM_output/mod_",scen,"_",bootpicks[i],".rds")) ## mod for boot 4
  load(file = paste0("./GAM_output/m2d_",scen,"_",bootpicks[i],".rds")) ## m2.d for boot 4 (contains derivative info)
  
  Terms <- c("Year","Latitude_dd","Longitude_dd")
  
  if(i < 5){ ## plot tempvar separately
  plist0[[idx0]]  <- ggplot(tempdf, aes(y = Latitude_dd, x = Longitude_dd)) +
    theme_classic() +
    # theme(legend.position = 'none') +
    scale_y_continuous(limits = c(-10,60), breaks = seq(0,50,10)) +
    geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
    scale_fill_gradient2(low = 'dodgerblue2', mid = "dodgerblue3",
                         high ="dodgerblue4", midpoint = 200, guide = 'legend') +
    scale_size_continuous(range = c(10, 15), guide = 'legend') +
    labs(x = 'Longitude',y = 'Latitude',
         fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
         title = paste0(scens.title0[i])) #+
    # geom_label(x = 10, y = 50, label = paste0('n = ',nrow(tempdf)), 
    #            size = 4, fill = 'white')
  idx0 <- idx0+1
  }
  
  if(i == 5){ ## plot tempvar separately
    plist0[[5]]  <- ggplot(tempdf, aes(x = Year, y = Length_cm)) +
      theme_classic() +
      # theme(legend.position = 'right') +
      scale_y_continuous(limits = c(125,450)) +
      geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
      scale_fill_gradient2(low = 'dodgerblue2', mid = "dodgerblue3",
                           high ="dodgerblue4", midpoint = 200, guide = 'legend') +
      scale_size_continuous(range = c(10, 15), guide = 'legend') +
      labs(x = 'Year', y = 'Length of Age Six Fish (cm)',
           fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
           title = paste0(scens.title0[i])) #+
  }
    ggsave(plot = last_plot(),  file = paste0("./figures/rawdat_",scen,".png"),
  width = 5, height = 5, units = 'in', dpi = 480)


  ## map showing new breaks
  if(i < 5){
    plist[[idx]] <- ggplot(tempdf, aes(y = Latitude_dd, x = Longitude_dd)) +
      theme_classic() +
      theme(legend.position = 'right') +
       scale_y_continuous(limits = c(0,50), breaks = seq(0,50,10)) +
      geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
      scale_fill_gradient2(low = 'dodgerblue2', mid = "dodgerblue3",
                           high ="dodgerblue4", midpoint = 200, guide = 'legend') +
      scale_size_continuous(range = c(10, 15), guide = 'legend') +
      labs(x = 'Longitude',y = 'Latitude',
           fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
           title = paste0(scens.title[i], ' Data and GAM-detected breaks ')) +
      geom_hline(yintercept = breaksdf$lat_breaks, lwd = 1.1, linetype = 'dashed', col = 'red') +
      geom_vline(xintercept = breaksdf$lon_breaks, lwd = 1.1, linetype = 'dashed', col = 'red')
    idx <- idx + 1
    }
    if(i == 5){ ## plot tempvar separately
    
      plist[[idx]] <- ggplot(tempdf, aes(x = Year, y = Length_cm)) +
        theme_classic() +
        theme(legend.position = 'right') +
        scale_y_continuous(limits = c(125,450)) +
        geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
        scale_fill_gradient2(low = 'dodgerblue2', mid = "dodgerblue3",
                             high ="dodgerblue4", midpoint = 200, guide = 'legend') +
        scale_size_continuous(range = c(10, 15), guide = 'legend') +
        labs(x = 'Year', y = 'Length of Age Six Fish (cm)',
             fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
             title = paste0(scens.title0[i])) +
        # geom_label(x = 10, y = 400, label = paste0('n = ',nrow(tempdf)), 
        #            size = 4, fill = 'white') +
        geom_vline(xintercept = breaksdf$yr_breaks, lwd = 1.1, linetype = 'dashed', col = 'red')
      idx <- idx + 1
      }
  # ## plot GAM smooth and deriv -- will have to be from single, sorta-representative boot
  for(t in 1:length(Terms)){
    pd <- plot(mod,   select = t, scheme  =2, lwd  =2, main = paste0(Terms[t],' Smoother'), cex.axis = 2, ylim = c(-10,ifelse(t != 3,10,500)))
    temp0 <- pd[t] ## get this smoother
    temp <- data.frame(cbind(temp0[[1]]$x,temp0[[1]]$fit, temp0[[1]]$se)); names(temp) = c('x','fit','se')
    if(t == 1){
      xlims <- c(0,100)
    } else{ xlims <- c(0,50)}
    ## plot smooth
    plist[[idx]] <-  ggplot(temp, aes(x = x, y = fit)) +
      theme_minimal() +
      theme(panel.grid = element_blank())+
      geom_line(lwd = 1.1) +
      scale_x_continuous(limits = xlims) +
      
      geom_line(aes(y= fit-se), linetype = 'dashed') +
      geom_line(aes(y= fit+se), linetype = 'dashed') +
      geom_rug(sides = 'b') +
      labs(x = Terms[t], y = "smoother", title = paste0(letters[idx-1],") ",'Smoother for ',Terms[t] ))
    idx <- idx+1
    CI <- confint(m2.d, term = Terms[t])
    m2.dtemp <- data.frame(cbind(m2.d$eval[,Terms[t]], m2.d[[Terms[t]]]$deriv,
                                 CI[[1]]$upper, CI[[1]]$lower)); names(m2.dtemp) = c('x','deriv','upper','lower')
    ## plot deriv
    # if(i < 5){
    plist[[idx]] <-  ggplot(m2.dtemp,aes(x = x, y = deriv))    +
      theme_minimal() +
      theme(panel.grid = element_blank())+
      geom_line(lwd = 1.1) +
      scale_x_continuous(limits = xlims) +
            geom_hline(yintercept = 0, col = 'grey22') +
      geom_line(aes(y= upper), linetype = 'dashed') +
      geom_line(aes(y= lower), linetype = 'dashed') +
      geom_vline(xintercept = breaksdf[,c('yr_breaks','lat_breaks','lon_breaks')[t]], 
                 col = 'red', lwd = 1.1, linetype = 'dashed') +
      labs(x = Terms[t], y = "f'(x)", title = paste0(letters[idx-1],") ",'First Derivative for ', Terms[t]))
    idx <- idx+1

  } ## end Terms


  lay <- rbind(c(2,2,1,1,1),
               c(3,3,1,1,1),
               c(4,4,1,1,1),
               c(5,5,1,1,1),
               c(6,6,1,1,1),
               c(7,7,1,1,1))
  grid.arrange(grobs = plist[1:7], layout_matrix = lay) %>%
  ggsave(plot = .,  file = paste0("./figures/GAM_analysis_",scen,".png"), width = 11, height = 8, units = 'in', dpi = 480)

}

ggarrange(plotlist = plist0, ncol=2, nrow=3, common.legend = TRUE, legend="bottom") %>%
  ggsave(plot = .,  file = paste0("./figures/rawdata_compile.png"), width = 11, height = 8, units = 'in', dpi = 480)


# a6df <- list.files("./IBM_output/datasets/", full.names = T) %>%
#   lapply(read.csv) %>%  
#   bind_rows()  %>%
#     filter(Age == 6)

# dim(a6df)/(100*length(unique(scenarios$DESC))) ## 100 datasets times 5 simulations -- getting average per ds


## GAM propagg----
cdfprop <- read.csv(paste0('./gam_output/cdf_prop_',Sys.Date()-2,'.csv')) 
levels(cdfprop$scen) <- c("Break at 25 deg.", "Break at 48 deg.", "Overlap 20-25 deg.","No Breaks","Temporal Break Year 50")

# cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
levels(cdfprop$variable) <- c('L1','L2' )
cdfprop2 <- cdfprop%>% filter(!(variable %in% c('Both L1 and L2')) & scen != "20% Higher k, Break at 25 deg.")

# cdfprop$variable <- factor(cdfprop$variable, levels=c('L1','L2','Both L1 and L2'))
plist1 <- list()
plist1[[1]] <- ggplot(cdfprop2, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = c(0.6,0.9)) +
  scale_fill_viridis_d()+
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Coverage Probability', fill = 'Scenario', 
       title = 'a) Coverage Probability for Endpoints of Growth Curve') +
  geom_bar(stat = 'identity',width=0.6, position = position_dodge(width=0.7)) +
  facet_wrap(~variable)
# ggsave(plot = last_plot(),  file = paste0("./figures/cdfprop.png"), width = 9, height = 6, units = 'in', dpi = 480)

cdfaccu <- read.csv(paste0('./gam_output/cdf_accu_',Sys.Date()-2,'.csv'))
levels(cdfaccu$scen) <-c("Break at 25 deg.", "Break at 48 deg.", "Overlap 20-25 deg.","No Breaks","Temporal Break Year 50")
# cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
levels(cdfaccu$variable) <- c('Latitude', 'Longitude' ,'Year')
# cdfaccu$variable <- factor(cdfaccu$variable, levels=c('L1','L2','Both L1 and L2'))

cdfaccu2 <- cdfaccu %>% filter(!(variable %in% c('Lat, Long and Year','Both Latitude and Longitude'))& scen != "20% Higher k, Break at 25 deg.")
# cdfaccu$scen  <- factor(cdfaccu$scen , levels = cdfaccu$scen [order(cdfprop$prop  )])
plist1[[2]] <- ggplot(cdfaccu2, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = 'none') +
  scale_fill_viridis_d()+
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Proportion Detected Accurate Breaks', fill = 'Scenario', title = 'b) Proportion Detected Accurate Breaks') +
  geom_bar(stat = 'identity',width=0.5, position = position_dodge(width=0.5)) +
  facet_wrap(~variable, nrow = 1)

ggarrange(plotlist = plist1, ncol=1, nrow=2, common.legend = TRUE, legend="bottom") %>%
  ggsave(plot = .,  file = paste0("./figures/GAM_cdfprob_",Sys.Date()-2,".png"), width = 11, height = 8, units = 'in', dpi = 480)


## STARS propagg ----
cdfprop <- read.csv(paste0('./stars_output/STARS_cdf_prop_2019-06-10.csv'))
levels(cdfprop$scen) <- c("Break at 25 deg.", "Break at 48 deg.", "Overlap 20-25 deg.","No Breaks","Temporal Break Year 50")
# cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
levels(cdfprop$variable) <- c('L1','L2','Both L1 and L2')
cdfprop2 <- cdfprop %>% filter(!(variable %in% c('Both L1 and L2'))& scen != "20% Higher k, Break at 25 deg.")

# cdfprop$variable <- factor(cdfprop$variable, levels=c('L1','L2','Both L1 and L2'))
plist1 <- list()
plist1[[1]] <- ggplot(cdfprop2, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = c(0.6,0.9)) +
  scale_fill_grey() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Coverage Probability', fill = 'Scenario', 
       title = 'a) Coverage Probability for Endpoints of Growth Curve') +
  geom_bar(stat = 'identity',width=0.6, position = position_dodge(width=0.7)) +
  facet_wrap(~variable)
# ggsave(plot = last_plot(),  file = paste0("./figures/cdfprop.png"), width = 9, height = 6, units = 'in', dpi = 480)

cdfaccu <- read.csv(paste0('./STARS_output/STARS_cdf_accu_2019-06-10.csv'))
levels(cdfaccu$scen) <- c("Break at 25 deg.", "Break at 48 deg.", "Overlap 20-25 deg.","No Breaks","Temporal Break Year 50")
# cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
levels(cdfaccu$variable) <- c('Latitude', 'Longitude' ,'Year')
# cdfaccu$variable <- factor(cdfaccu$variable, levels=c('L1','L2','Both L1 and L2'))

cdfaccu2 <- cdfaccu %>% filter(!(variable %in% c('Lat, Long and Year','Both Latitude and Longitude'))& scen != "20% Higher k, Break at 25 deg.")
# cdfaccu$scen  <- factor(cdfaccu$scen , levels = cdfaccu$scen [order(cdfprop$prop  )])
plist1[[2]] <- ggplot(cdfaccu2, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = 'none') +
  # scale_fill_viridis_d()+
  scale_fill_grey() +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Proportion Detected Accurate Breaks', fill = 'Scenario', title = 'b) Proportion Detected Accurate Breaks') +
  geom_bar(stat = 'identity',width=0.5, position = position_dodge(width=0.5)) +
  facet_wrap(~variable, ncol = 3)

ggarrange(plotlist = plist1, ncol=1, nrow=2, common.legend = TRUE, legend="bottom") %>%
  ggsave(plot = .,  file = paste0("./figures/STARS_cdfprob_",Sys.Date()-2,".png"), width = 11, height = 8, units = 'in', dpi = 480)


rbind(cdfprop, cdfaccu) %>% 
  group_by(scen) %>%
  
  summarise(mean(prop)) 



## GAM Histogram of detected breaks----
cdf_gam <- read.csv("GAM_output/cdf_2019-06-10.csv") %>% filter(scen != 'F0L1S_R3')
ntrue <- read.csv("input_data/ntrue_a6.csv", na.strings = 'NA') %>% filter(scen != 'F0L1S_R3') %>% mutate(value = 1)
levels(cdf_gam$scen) <- levels(ntrue$scen) <- c("Scenario 2", "Scenario 4", 
                                                "Scenario 3","Scenario 1",
                                                "Scenario 5")

plist <- list(); idx = 1
for(l in 1:length(unique(cdf_gam$scen))){
  scentemp <- unique(cdf_gam$scen)[l]
  # scen <- paste0(ldfprop$scen[ldfprop$scentemp == scentemp]) ## filenames
  cat(l,idx,"\n")
  tmp0 <- cdf_gam %>%  filter(scen == scentemp) %>% select(scen, gamLAT, gamLON, gamYR) %>% melt(id = c('scen'))
  levels(tmp0$variable) <- c('Latitude','Longitude','Year')
  tmp0$value <- factor(tmp0$value, levels=c(paste(1:100),NA))
  
  tmp <- tmp0 %>%
    group_by(variable, value) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n), scen = scentemp) 
  
  ## plot lat long first, then year
  plist[[idx]]  <- ggplot(subset(tmp, variable != 'Year'), aes(x = value, y = freq)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          strip.text = element_text(size= 9))+
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = c(paste(1:51),NA),
                     breaks = c(paste(seq(0,51,10)),NA)) +
    geom_rect(aes(xmin = ifelse(variable != 'Year', 20,-10),
                  xmax = ifelse(variable != 'Year', 25,-10), ymin = 0, ymax = Inf),
              fill=ifelse(l == 3,"red",NA), alpha = 1E-2) +
    facet_wrap(~scen+variable,ncol = 2,drop=TRUE) +
    labs( y = 'frequency', x = 'detected break location', title = '') +
    geom_bar(data = subset(ntrue, scen == scentemp & variable != 'Year'),
               aes(x = factor(trueb), y = value), col = 'red', stat = 'identity', alpha = 0.5) 
  idx <- idx +1
  
  plist[[idx]]  <- ggplot(subset(tmp, variable == 'Year'), aes(x = value, y = freq)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          strip.text = element_text(size= 9))+
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = c(paste(1:101),NA),
                     breaks = c(paste(seq(0,101,10)),NA)) +
    geom_rect(aes(xmin = ifelse(variable != 'Year', 20,-10),
                  xmax = ifelse(variable != 'Year', 25,-10), ymin = 0, ymax = Inf), 
              fill=ifelse(l == 3,"red",NA), alpha = 1E-2) +
    facet_wrap(~scen+variable,ncol = 3) +
    labs( y = 'frequency', x = 'detected break location', title = '') +
    geom_bar(data =  subset(ntrue, scen == scentemp & variable == 'Year'),
             aes(x = factor(trueb), y = value), col = 'red', stat = 'identity', alpha = 0.5) 
  idx <- idx +1
}

grid.arrange(grobs = plist, ncol = 4) %>%
  ggsave(plot = .,   file = paste0("./figures/GAM_hist_breaks.png"), width = 10, height = 12, units = 'in', dpi = 480)


## STARS Histogram of detected breaks----
cdf_stars <- read.csv("STARS_output/STARS_cdf_2019-06-10.csv") %>% filter(scen != 'F0L1S_R3')
ntrue <- read.csv("input_data/ntrue_a6.csv", na.strings = 'NA') %>% filter(scen != 'F0L1S_R3') %>% mutate(value = 1)
levels(cdf_stars$scen) <- levels(ntrue$scen) <- c("Scenario 2", "Scenario 4", 
                                                "Scenario 3","Scenario 1",
                                                "Scenario 5")
plist <- list(); idx = 1
for(l in 1:length(unique(cdf_stars$scen))){
  scentemp <- unique(cdf_stars$scen)[l]
  # scen <- paste0(ldfprop$scen[ldfprop$scentemp == scentemp]) ## filenames
  cat(l,idx,"\n")
  tmp0 <- cdf_stars %>%  filter(scen == scentemp) %>% select(scen, gamLAT, gamLON, gamYR) %>% melt(id = c('scen'))
  levels(tmp0$variable) <- c('Latitude','Longitude','Year')
  tmp0$value <- factor(tmp0$value, levels=c(paste(1:100),NA))
  
  tmp <- tmp0 %>%
    group_by(variable, value) %>%
    dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n), scen = scentemp) 
  
  ## plot lat long first, then year
  plist[[idx]]  <- ggplot(subset(tmp, variable != 'Year'), aes(x = value, y = freq)) +
    geom_bar(stat = 'identity', col = 'skyblue') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          strip.text = element_text(size= 9))+
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = c(paste(1:51),NA),
                     breaks = c(paste(seq(0,51,10)),NA)) +
    geom_rect(aes(xmin = ifelse(variable != 'Year', 20,-10),
                  xmax = ifelse(variable != 'Year', 25,-10), ymin = 0, ymax = Inf),
              fill=ifelse(l == 3,"red",NA), alpha = 1E-2) +
    facet_wrap(~scen+variable,ncol = 2,drop=TRUE) +
    labs( y = 'frequency', x = 'detected break location', title = '') +
    geom_bar(data = subset(ntrue, scen == scentemp & variable != 'Year'),
             aes(x = factor(trueb), y = value), col = 'red', stat = 'identity', alpha = 0.5) 
  idx <- idx +1
  
  plist[[idx]]  <- ggplot(subset(tmp, variable == 'Year'), aes(x = value, y = freq)) +
    geom_bar(stat = 'identity', col = 'skyblue') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          strip.text = element_text(size= 9))+
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = c(paste(1:101),NA),
                     breaks = c(paste(seq(0,101,10)),NA)) +
    geom_rect(aes(xmin = ifelse(variable != 'Year', 20,-10),
                  xmax = ifelse(variable != 'Year', 25,-10), ymin = 0, ymax = Inf), 
              fill=ifelse(l == 3,"red",NA), alpha = 1E-2) +
    facet_wrap(~scen+variable,ncol = 3) +
    labs( y = 'frequency', x = 'detected break location', title = '') +
    geom_bar(data =  subset(ntrue, scen == scentemp & variable == 'Year'),
             aes(x = factor(trueb), y = value), col = 'red', stat = 'identity', alpha = 0.5) 
  idx <- idx +1
}

grid.arrange(grobs = plist, ncol = 4) %>% 
  ggsave(plot = .,   file = paste0("./figures/STARS_hist_breaks.png"), width = 10, height = 12, units = 'in', dpi = 480)

breaksdf <- read.csv(paste0("./stars_output/STARS_breaksdf_2019-04-29.csv"))

## Trajectories of individual fish from IBM ----
# L <- list.files("./IBM_output/", full.names = T, recursive = T)[grep("_55/IBM_SAA", list.files("./IBM_output/", full.names = T,recursive = T))]  
# O = lapply(L, function(x) {
#   DF <- read.csv(x, header = T, sep = ",")
#   DF$ID <- paste0(basename(dirname(dirname(dirname(x)))))
#   return(DF)})
# df1 <- do.call(rbind, O)
# df1$REG[df1$REG == 'R0'] <- 'R1'
# df2 <- df1 %>% group_by(ID,Age,REG) %>% summarise(meanL = mean(fish_size))
# df2 <- subset(df1,  cohort = 44) %>% 
#   filter(!is.na(REG) & REG != 'R3'  & ID == 'F0L1S_25')
#   # group_by(ID, Age) %>% 
#   # summarise(meanL = mean(fish_size), minL = min(fish_size), maxL = max(fish_size))
# df3 <- read.csv("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/IBM_output/datasets/F0L1S_25_44.csv") %>% filter( Length_cm < 275)
# ggplot(df3, aes(x = Age, y = Length_cm, color = REG)) +
#   theme_classic()+
#   theme(legend.position = c(0.9,0.1), 
#         legend.text = element_text(size = 14), 
#         legend.title = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14)) +
#   scale_color_viridis_d(guide = "legend") +
#   geom_point(size = 3, alpha = 0.2) +
#   # geom_line() +
#   labs(x = 'Age (years)', y = 'Length (cm)', color = "Growth Regime") 
# ggsave(plot = last_plot(),  file = paste0("./figures/ibm_growth.png"), width = 8, height = 6, units = 'in', dpi = 480)
a=

## Make Appendix Table A2 ----

## read stuff
cdf_gam <- read.csv(paste0("GAM_output/cdf_",Sys.Date()-2,".csv"))
cdfaccu_gam <- read.csv(paste0("GAM_output/cdf_accu_",Sys.Date()-2,".csv"))
cdfprop_gam <- read.csv(paste0("GAM_output/cdf_prop_",Sys.Date()-2,".csv"))

cdfaccu_threeq <- read.csv(paste0("GAM_output/cdf_accu_",Sys.Date()-2,"threequ",".csv"))
cdfaccu_half <- read.csv(paste0("GAM_output/cdf_accu_",Sys.Date()-2,"half",".csv"))

cdfprop_threeq <-read.csv(paste0("GAM_output/cdf_prop_",Sys.Date()-2,"threequ",".csv"))
cdfprop_half <-read.csv(paste0("GAM_output/cdf_prop_",Sys.Date()-2,"half",".csv"))

TA2 <- data.frame()
for(i in 1:length(unique(cdfaccu_gam$scen))){
  TA2[i,"scen"] <- unique(cdfaccu_gam$scen)[i]
  
  ## OG values [first two cols]
  tmp_accu <- subset(cdfaccu_gam, scen == unique(cdfaccu_gam$scen)[i])
  tmp_prop <- subset(cdfprop_gam, scen == unique(cdfprop_gam$scen)[i])
  TA2[i,"L1L2_orig"] <- paste(tmp_prop$prop, collapse = ", ")
  TA2[i,"latlonyr_orig"] <- paste(tmp_accu$prop, collapse = ", ")
  
  ## Halved Sens
  tmp_accu <- subset(cdfaccu_half, scen == unique(cdfaccu_half$scen)[i])
  tmp_prop <- subset(cdfprop_half, scen == unique(cdfprop_half$scen)[i])
  TA2[i,"L1L2_half"] <- paste(tmp_prop$prop, collapse = ", ")
  TA2[i,"latlonyr_half"] <- paste(tmp_accu$prop, collapse = ", ")
  
  ## 3/4 Sens
  tmp_accu <- subset(cdfaccu_threeq, scen == unique(cdfaccu_threeq$scen)[i])
  tmp_prop <- subset(cdfprop_threeq, scen == unique(cdfprop_threeq$scen)[i])
  TA2[i,"L1L2_3q"] <- paste(tmp_prop$prop, collapse = ", ")
  TA2[i,"latlonyr_3q"] <- paste(tmp_accu$prop, collapse = ", ")
}
TA2$ORD <- c(2,4,3,1,5)
TA2[order(TA2$ORD),] %>% select(-ORD) %>% write.csv(paste0("figures/table_a2_",Sys.Date()-2,".csv"), row.names = F)
# cdf_stars <- read.csv(paste0("STARS_output/STARS_cdf_",Sys.Date()-2,".csv"))
# cdfaccu_stars <- read.csv(paste0("STARS_output/STARS_cdf_accu_",Sys.Date()-2,".csv"))
# cdfprop_stars <- read.csv(paste0("STARS_output/STARS_cdf_prop_",Sys.Date()-2,".csv"))



## GAM/STARS summarise MISSED proportions ----
# cdf_gam <- read.csv(paste0("GAM_output/cdf_2019-06-10.csv"))
# cdfaccu_gam <- read.csv(paste0("GAM_output/cdf_accu_2019-05-28.csv"))
# cdfprop_gam <- read.csv(paste0("GAM_output/cdf_prop_2019-05-28.csv"))
# 
# cdf_stars <- read.csv(paste0("STARS_output/STARS_cdf_",Sys.Date()-2,".csv"))
cdfaccu_stars <- read.csv(paste0("STARS_output/STARS_cdf_accu_2019-06-10.csv"))
cdfprop_stars <- read.csv(paste0("STARS_output/STARS_cdf_prop_2019-06-10.csv"))
# 
# cdf_gam %>% 
#   filter(L1 == FALSE) %>%
#   group_by(scen) %>% 
#   summarise( meanL1 = mean(abs(L1_MISS)))
# 
# cdf_gam %>% 
#   filter(L2 == FALSE) %>%
#   group_by(scen) %>% 
#   summarise( meanL2 = mean(abs(L2_MISS)))
# 
# cdf_stars %>% 
#   filter(L1 == FALSE) %>%
#   group_by(scen) %>% 
#   summarise( meanL1 = mean(abs(L1_MISS)))
# 
# cdf_stars %>% 
#   filter(L2 == FALSE) %>%
#   group_by(scen) %>% 
#   summarise( meanL2 = mean(abs(L2_MISS)))
# 
# rbind(cdf_gam,cdf_stars) %>% group_by(method,scen, YEAR) %>% dplyr::summarise(n=n())
# cdf_gam %>% filter(YEAR == FALSE) %>% group_by(scen, gamYR) %>% dplyr::summarise(n=n()) %>% 
#   ggplot(., aes(x = gamYR, y = n)) + geom_histogram(stat='identity')
# 
# 
# cdf_gam %>% 
#   filter(scen == 'F0L1S_25' & LAT == FALSE) %>% 
#   group_by(gamLAT) %>% summarise(n=n())

## many false edges were nearby
# cdf_gam %>%
#   filter(scen == 'F0L1S_49' & (LON == FALSE | LAT == FALSE)) %>%
#   group_by(gamLAT,gamLON) %>% dplyr::summarise(n=n())
# cdf_gam %>% 
#   filter(scen == 'F0L1S_49' & LAT == TRUE) %>% 
#   group_by(gamLAT) %>% dplyr::summarise(n=n())


## "no discernable pattern in spurious breaks"
# cdf %>% filter(scen %in% c('NoBreaks','tempvar_R1R2') & (LAT == FALSE  | LON == FALSE)) %>% group_by(scen, gamLAT,gamLON) %>% dplyr::summarise(n = n())
# ## what proportion of wrong ones are only one off?
# cdf %>% filter(!(scen %in% c('NoBreaks','tempvar_R1R2')) & (LAT == FALSE  | LON == FALSE)) %>% 
#   group_by(scen) %>% dplyr::summarise(n = n())
# # 
cdf_gam %>% filter(scen %in% c('F0L1S_25','F0L1S_R3') & (LAT == FALSE  | LON == FALSE)) %>%
  group_by(scen) %>% dplyr::summarise(n = n(),
                                      inbounds = sum(gamLAT %in% c(24,26) |gamLON %in% c(24,26)  )) %>% mutate(inbounds/n)

# cdf %>% filter(YEAR == FALSE) %>% group_by(scen, gamYR) %>% dplyr::summarise(n = n())
# cdf %>% filter(scen == 'F0LMW'& LAT == FALSE) %>% group_by(gamLAT) %>% dplyr::summarise(n = n())
# cdf %>% group_by(YEAR) %>% dplyr::summarise(n = n()) 
# 
# 
# ## 16+156 /280 is 61% correct for YEAR when relaxed (49 or 51)
## 116 [50] + 15+156 / 280+116 [total] = 73% between 49 and 51
cdf_gam %>% filter(scen == 'tempvar_R1R2') %>% group_by(YEAR) %>% dplyr::summarise(n = n())
cdf_gam %>% filter(scen == 'tempvar_R1R2' & YEAR == FALSE & gamYR %in% c(49,51)) %>% group_by(gamYR) %>% dplyr::summarise(n = n())
# 
# cdf_gam %>% filter(scen == 'tempvar_R1R2' & YEAR == FALSE ) %>% group_by(gamYR) %>% dplyr::summarise(n = n())
# 
# 
# 
# ## percent accuracy by cat
cdfaccu_stars %>% filter(scen != 'F0L1S_R3' & scen != 'F0L1S_49') %>% group_by(variable) %>%  dplyr::summarise(avg = mean(prop))
cdfaccu_gam %>% filter(scen != 'F0L1S_R3' & scen != 'F0L1S_49') %>% group_by(variable) %>%  dplyr::summarise(avg = mean(prop))


cdfprop_stars %>% filter(scen != 'F0L1S_R3' & scen != 'F0L1S_49') %>% group_by(variable) %>% dplyr::summarise(avg = mean(prop))
cdfprop_gam %>%  filter(scen != 'F0L1S_R3' & scen != 'F0L1S_49') %>%group_by(variable) %>%  dplyr::summarise(avg = mean(prop))

## SAB STUFF ## ----

## average # of age 6 fish per dataset----
# a6df <- list.files("./IBM_output/datasets", full.names = T) %>%
#   lapply(read.csv) %>%
#   bind_rows() %>%
#   filter(Age == 6)
# dim(a6df)/(100*length(unique(scenarios$DESC))) ## 100 datasets times 5 simulations -- getting average per ds

## average # of age X SAB per sex in fulldat
# load(paste0("./input_data/gam_data_sab_0415.rda")) ## full_data -- made using gam_dataprep NOT 15k subsample
# full_data %>% filter(Age %in% c(4,6,10,30)) %>% 
#   group_by(Age, Sex) %>% summarise(n=n()) %>% write.csv(.,paste0("./output_data/n_sab_sex_age.csv"),row.names=F)


## Figure 8 SAB map for a6 and a10 breaks, likely redone by Elliot ----
# usa <- map_data("world")
#   dat <- all_data %>% filter(Age == 30 & Sex == 'M')
#     ggplot() + 
#     geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
#     coord_quickmap() +
#     scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
#     scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
#     theme_minimal() +
#     theme(panel.grid.major = element_blank(),
#           axis.title =element_blank(), 
#           legend.position = c(0.1,0.2),
#           legend.key.size = unit(0.75, "cm"),
#           legend.text = element_text(size = 12)) +
#     geom_hline(yintercept = c(36,50),lwd = 1.1, linetype = 'dashed', col = 'red') +
#     geom_vline(xintercept =  -130, lwd = 1.1, linetype = 'dashed', col = 'red') +
#       ## ecosystem break
#     geom_vline(xintercept =  -145, lwd = 1.1, linetype = 'dashed', col = 'blue') +
#     geom_rect(fill = 'white', aes(xmin = -180, xmax = -130.1, ymin = 30, ymax = 49.9)) +
#     geom_point(data = dat, aes(x = Longitude_dd, y = Latitude_dd, size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
#     scale_fill_viridis_c(guide = "legend") +
#     labs(fill = paste0("Length (cm)"),
#          size = paste0("Length (cm)"),
#          title = "GAM- and Ecosystem-based Regions with Age 30 Fish") +
#       geom_label(aes(x = c(-155,-138,-120,-120,-120), 
#                     y = c(65,65,65,40,32),
#                     label = c(paste('Region ',5:1)) ),size = 6, col = 'black',fill = 'white',show.legend   =FALSE)
#     ggsave(plot = last_plot(),  file = paste0("./figures/sab_zones.png"), width = 8, height = 6, units = 'in', dpi = 480)
#     
# df1<-all_data %>% filter(Age %in% c(4,6,30)) %>% 
#   group_by(REG, Age) %>% summarise(n = n()) 
# 
# df1 %>%
#   group_by(Age) %>% summarise(mn = mean(n))

## predicts and parest for sab ----
## parest: see if actually different. PHASE 1 = pre-merge, PHASE 2- temporal merges
# for(phase in c("phase1","phase2")){
#   parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_2019-04-15_",phase,'.csv')) %>% 
#     filter(variable ==  "Linf") %>% 
#     mutate(source = 'Estimated') %>%
#     mutate(REG2 = gsub("_.*", "\\1", REG),
#            REG3 = as.numeric(substr(REG,2,2)),
#            REG4 = sub('_([^_]*)$', '',REG),
#            Sex = gsub(".*_", "\\1", REG),
#            lwr = value - 1.96 * sd,
#            upr = value + 1.96 * sd,
#            matchcol = 'black')
#   # parest$value <- exp(parest$value - (parest$sd^2)/2) ## bias correction
#   parest <- parest[order(parest$REG2,parest$Sex),]
#   
#   # parest <- rbind(parest, read.csv("./input_data/true_sab_vals.csv")) 
#   # levels(parest$REG) <- c("ALL","R1","AK","R2","BC","R1","WC")
#   # parest$REG <- factor(parest$REG ,levels=c("ALL","R1","WC","R2","BC","R3","AK",'R4'))
#   
#   ## check if CIs overlap
#   if(phase == 'phase1'){
#     for(i in seq(1,nrow(parest),2)) { ## will go by sex
#       parest$match[i:(i+1)] <- ifelse(c(parest$lwr[i], parest$upr[i]) %overlaps%  c(parest$lwr[i+1], parest$upr[i+1]),
#                                       paste0('OVERLAP ', parest$REG2[i+1]), 'NO OVERLAP')
#     }
#     
#   } else if(phase == 'phase2') {
#     for(i in 1:nrow(parest)){ ## we want to compare w adjacent regions
#       if(parest$REG3 < 5){
#         tmp <- subset(parest, Sex == parest$Sex[i] & REG3 == parest$REG3[i]+1)
#         parest$match[i] <-ifelse(any(c(tmp$lwr, tmp$upr) %overlaps% c(parest$lwr[i], parest$upr[i])),paste0('OVERLAP ', tmp$REG2[1]),'NO OVERLAP')
#       } ## end R < 5
#     } ## end parest rows
#   } ## end else phase2
#   
#   parest$matchcol <- parest$match != 'NO OVERLAP'
#   parest$Sex <- factor(parest$Sex)
#   levels(parest$Sex) <- c('Females','Males')
#   ggplot(parest, aes(x = REG4, y = value, col = matchcol))+
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           legend.position = 'right',
#           legend.background = element_blank(),
#           axis.text = element_text(size = 10,angle = 45),
#           axis.title = element_text(size = 10),
#           legend.text = element_text(size = 10),
#           strip.text = element_text(size=14))+
#     scale_y_continuous(limits = c(0,100)) +
#     scale_color_manual(values = c('black','red')) +
#     geom_point() +
#     geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
#     labs(x = 'Spatiotemporal x Sex Stratum', y = "", 
#          col = ifelse(phase == 'phase1', 'CI overlap within Region + Sex','CI Overlap Adjacent Region x Sex'),
#          title = paste0(phase," Linf Estimates")) +
#     facet_wrap(~Sex )
# 
# ggsave(plot = last_plot(),  file = paste0("./figures/sab_parest_",Sys.Date()-2,"_",phase,".png"), width = 10, height = 8, units = 'in', dpi = 480)
# write.csv(parest, file = paste0("./GAM_output/overlap_",Sys.Date()-2,"_",phase,".csv"),row.names=F)


## Figure 9 SAB Fits ----
# ypreds <- read.csv(paste0("./GAM_output/SAB_predicts_2019-04-15_",phase,".csv"))
# 
# ypreds$gamREG <- paste0('Region ',ypreds$gamREG)
# levels(ypreds$Sex) <- c('Females','Males')
# levels(ypreds$Period) <- c('pre-2010','2010-Present','All Years')
# for(i in 1:nrow(ypreds)){
#   ypreds$Period[i] <-  ifelse(length(grep('pool', ypreds$cREG[i])) == 1, 
#                               'All Years', paste(ypreds$Period[i]))
# }
# 
# ggplot(ypreds, aes(x = Age, y = Predicted, col = REG, linetype = Period )) +
#   theme_classic() +
#   theme(panel.grid = element_blank(),
#         legend.position = 'right',
#         legend.background = element_blank(),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         legend.text = element_text(size = 10),
#         strip.text = element_text(size=10))+
#   scale_alpha(guide = 'none') +
#   scale_y_continuous(limits = c(0,110)) +
#   scale_color_brewer(palette =  'Accent')+
#   geom_point(alpha = 0.5, aes(y = Length_cm)) +
#   geom_line(lwd = 1.1, col = 'black')+
#   labs(y = 'Length (cm)', x= 'Age (years)', col = "Actual Data Source") +
#   scale_linetype_manual(values=c("solid", "dashed", "dotted"))+
#   facet_wrap(~gamREG + Sex , ncol = 4)
# 
# ggsave(plot = last_plot(),  
#        file = paste0("./figures/sab_fits_",Sys.Date()-2,"_",phase,".png"), 
#        width = 10, height = 12, units = 'in', dpi = 520)
# cat(phase," done \n")
# } ## end phase



# all_data %>% filter(Age %in% c(4,6,30)) %>% group_by(Age,Sex,REG) %>% summarise(n = n())
