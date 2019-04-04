## PlotMaster
## Figures for Manuscript/Chapter 1

require(dplyr)
require(ggplot2)
require(gridExtra)
require(RColorBrewer)

scenarios <- read.csv('./input_data/scenarios.csv',na.strings = 'NA') ## manual file
age <- 6
## average ## of age X fish per dataset

#   a6df <- list.files("./IBM_output/datasets", full.names = T) %>%
#   lapply(read.csv) %>%
#   bind_rows() %>%
#   filter(Age == 6)
# dim(a6df)/(100*length(unique(scenarios$DESC))) ## 100 datasets times 5 simulations -- getting average per ds

## how many NA years (should be most!!)
ldf %>% group_by(scen) %>% summarise(naY = sum(is.na(yr_breaks)), n = n()) %>% mutate(naY/n)


## panel plot of one gendata for each scenario ---
scens <-  unique(scenarios$DESC)
scens.title <- scens.title0 <- scens
levels(scens.title) <-  c("e) Break at 25 deg.", "e) Break at 49 deg.",
                           "e) Low Contrast at 25 deg.", "e) Overlap 20-25 deg.","e) No Breaks")
levels(scens.title0) <-  c("b) Break at 25 deg.", "e) Break at 49 deg.",
                          "c) Low Contrast at 25 deg.", "d) Overlap 20-25 deg.","a) No Breaks")

plist0 <- list(); idx0 <- 1
for(i in 1:length(scens)){
  plist  <- list(); idx <- 1
  ## plot raw data
  # num <- floor(runif(1,1,100))
  scen <- scens[i]
  tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_4.csv")) %>% filter(Age == 6)

  breaksdf <- read.csv( paste0("./GAM_output/ldf_raw_a",age,".csv")) %>% filter(scen == scens[i] & boot == 4)
  load(file = paste0("./GAM_output/mod_",scen,"_4.rds")) ## mod for boot 4
  load(file = paste0("./GAM_output/m2d_",scen,"_4.rds")) ## m2.d for boot 4 (contains derivative info)
  
  Terms <- c("Year","Latitude_dd","Longitude_dd")
  
  ## raw data save separately
  plist0[[idx0]]  <- ggplot(tempdf, aes(x = Latitude_dd, y = Longitude_dd)) +
    theme_classic() +
    theme(legend.position = ifelse(idx0 == 5,'right','none')) +
    scale_y_continuous(limits = c(0,50)) +
    geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
    scale_fill_viridis_c(guide = "legend") +
    scale_size_continuous(range = c(1, 15)) +
    labs(x = 'Longitude',y = 'Latitude',
         fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
         title = paste0(scens.title0[i], ' Raw Data')) +
    geom_label(x = 2, y = 50, label = paste0('n = ',nrow(tempdf)), 
               size = 4, fill = 'white')
  idx0 <- idx0+1

  #   ggsave(plot = last_plot(),  file = paste0("./figures/rawdat_",scen,".png"), width = 5, height = 5, units = 'in', dpi = 480)


  ## map showing new breaks
  plist[[idx]] <- ggplot(tempdf, aes(x = Latitude_dd, y = Longitude_dd)) +
    theme_classic() +
    theme(legend.position = 'right') +
    scale_y_continuous(limits = c(0,50)) +
    geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
    scale_fill_viridis_c(guide = "legend") +
    scale_size_continuous(range = c(1, 15)) +
    labs(x = 'Longitude',y = 'Latitude',
         fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
         title = paste0(scens.title[i], ' Data and GAM-detected breaks ')) +
    geom_hline(yintercept = breaksdf$lat_breaks, lwd = 1.1, linetype = 'dashed', col = 'red') +
    geom_vline(xintercept = breaksdf$lon_breaks, lwd = 1.1, linetype = 'dashed', col = 'red')
  idx <- idx + 1
  # ## plot GAM smooth and deriv -- will have to be from single, sorta-representative boot
  for(t in 2:length(Terms)){
    pd <- plot(mod,   select = t, scheme  =2, lwd  =2, main = paste0(Terms[t],' Smoother'), cex.axis = 2, ylim = c(-10,ifelse(t != 3,10,500)))
    temp0 <- pd[t] ## get this smoother
    temp <- data.frame(cbind(temp0[[1]]$x,temp0[[1]]$fit, temp0[[1]]$se)); names(temp) = c('x','fit','se')

    ## plot smooth
    plist[[idx]] <-  ggplot(temp, aes(x = x, y = fit)) +
      theme_minimal() +
      theme(panel.grid = element_blank())+
      geom_line(lwd = 1.1) +
      geom_line(aes(y= fit-se), linetype = 'dashed') +
      geom_line(aes(y= fit+se), linetype = 'dashed') +
      geom_rug(sides = 'b') +
      labs(x = Terms[t], y = "smoother", title = paste0(letters[idx-1],") ",'Smoother for ',Terms[t] ))
    idx <- idx+1
    CI <- confint(m2.d, term = Terms[t])
    m2.dtemp <- data.frame(cbind(m2.d$eval[,Terms[t]], m2.d[[Terms[t]]]$deriv, CI[[1]]$upper, CI[[1]]$lower)); names(m2.dtemp) = c('x','deriv','upper','lower')
    ## plot deriv
    plist[[idx]] <-  ggplot(m2.dtemp,aes(x = x, y = deriv))    +
      theme_minimal() +
      theme(panel.grid = element_blank())+
      geom_line(lwd = 1.1) +
      geom_hline(yintercept = 0, col = 'grey22') +
      geom_line(aes(y= upper), linetype = 'dashed') +
      geom_line(aes(y= lower), linetype = 'dashed') +
      geom_vline(xintercept = breaksdf[,c(NA,'lat_breaks','lon_breaks')[t]], col = 'red', lwd = 1.1, linetype = 'dashed') +
      labs(x = Terms[t], y = "f'(x)", title = paste0(letters[idx-1],") ",'First Derivative for ', Terms[t]))
    idx <- idx+1
  }


  lay <- rbind(c(2,2,1,1,1),
               c(3,3,1,1,1),
               c(4,4,1,1,1),
               c(5,5,1,1,1)) #,
               # c(6,6,1,1,1),
               # c(7,7,1,1,1))
  grid.arrange(grobs = plist, layout_matrix = lay) %>%
  ggsave(plot = .,  file = paste0("./figures/analysis_p1_",scen,".png"), width = 11, height = 8, units = 'in', dpi = 480)


}

lay <- rbind(c(1,2,3),
             c(1,2,3),
             c(4,5,5),
             c(4,5,5))
grid.arrange(grobs = plist0, layout_matrix = lay) %>%
  ggsave(plot = .,  file = paste0("./figures/rawdata_compile.png"), width = 11, height = 8, units = 'in', dpi = 480)
  


a6df <- list.files("./IBM_output/datasets/", full.names = T) %>%
  lapply(read.csv) %>%  
  bind_rows()  %>%
    filter(Age == 6)

dim(a6df)/(100*length(unique(scenarios$DESC))) ## 100 datasets times 5 simulations -- getting average per ds

## Plot proportion agreegments ----
cdfprop <- read.csv(paste0('./gam_output/cdf_prop_',Sys.Date(),'.csv'))
levels(cdfprop$scen) <-c("Break at 25 deg.", "Break at 49 deg.",
                                                  "Low Contrast at 25 deg.", 
                                  "Overlap 20-25 deg.","No Breaks")
cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
plist1 <- list()
plist1[[1]] <- ggplot(cdfprop, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 8),
        legend.position = c(0.9,0.75)) +
  scale_fill_grey()+
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Coverage Probability', fill = 'Spatial Scenario', 
       title = 'a) Coverage Probability for Endpoints of Growth Curve') +
  geom_bar(stat = 'identity',width=0.6, position = position_dodge(width=0.7)) +
  facet_wrap(~variable)
# ggsave(plot = last_plot(),  file = paste0("./figures/cdfprop.png"), width = 9, height = 6, units = 'in', dpi = 480)



cdfaccu <- read.csv(paste0('./gam_output/cdf_accu_',Sys.Date(),'.csv'))
levels(cdfaccu$scen) <-c("Break at 25 deg.", "Break at 49 deg.",
                         "Low Contrast at 25 deg.",
                         "Overlap 20-25 deg.","No Breaks")
levels(cdfaccu$variable) <- c('Both Latitude and Longitude','Latitude','Longitude' )

# cdfaccu$scen  <- factor(cdfaccu$scen , levels = cdfaccu$scen [order(cdfprop$prop  )])
plist1[[2]] <- ggplot(cdfaccu, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text( size = 8),
        legend.position = 'none') +
  scale_fill_grey()+
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Proportion Detected Accurate Spatial Breaks', fill = 'Spatial Scenario', title = 'b) Proportion Detected Accurate Breaks') +
  geom_bar(stat = 'identity',width=0.5, position = position_dodge(width=0.5)) +
  facet_wrap(~variable)

# ggsave(plot = last_plot(),  file = paste0("./figures/cdfaccu.png"), width = 9, height = 6, units = 'in', dpi = 480)

lay <- rbind(c(1,1,1,1), 
             c(2,2,2,2))
grid.arrange(grobs = plist1, layout_matrix = lay) %>%
  ggsave(plot = .,  file = paste0("./figures/cdfprob_",Sys.Date(),".png"),
         width = 15, height = 10, units = 'in', dpi = 480)


rbind(cdfprop, cdfaccu) %>% group_by(scen) %>%
  summarise(mean(prop)) 
## Trajectories of individual fish from IBM ----
L <- list.files("./IBM_output/", full.names = T, recursive = T)[grep("_55/IBM_SAA", list.files("./IBM_output/", full.names = T,recursive = T))]  
O = lapply(L, function(x) {
  DF <- read.csv(x, header = T, sep = ",")
  DF$ID <- paste0(basename(dirname(dirname(dirname(x)))))
  return(DF)})
df1 <- do.call(rbind, O)
df1$REG[df1$REG == 'R0'] <- 'R1'
df2 <- df1 %>% group_by(ID,Age,REG) %>% summarise(meanL = mean(fish_size))
df2 <- subset(df1,  cohort = 44) %>% filter(!is.na(REG))
  # group_by(ID, Age) %>% 
  # summarise(meanL = mean(fish_size), minL = min(fish_size), maxL = max(fish_size))

ggplot(df2, aes(x = Age, y = fish_size, color = REG, group = cohort)) +
  theme_classic()+
  theme(legend.position = c(0.9,0.1), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  scale_color_viridis_d(guide = "legend") +
  geom_point(size = 3, alpha = 0.2) +
  labs(x = 'Age (years)', y = 'Length (cm)', color = "Growth Regime") 
ggsave(plot = last_plot(),  file = paste0("./figures/ibm_growth.png"), width = 8, height = 6, units = 'in', dpi = 480)
