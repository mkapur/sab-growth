## PlotMaster
## Figures for Manuscript/Chapter 1

require(dplyr)
require(ggplot2)
require(gridExtra)
require(RColorBrewer)
library(ggpubr)

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
levels(scens.title) <-  c("g) Break at 25 deg.", "g) Break at 49 deg.",
                           "g) Low Contrast at 25 deg.", "g) Overlap 20-25 deg.","g) No Breaks","f) Temporal Break Year 50")
levels(scens.title0) <-  c("b) Break at 25 deg.", "e) Break at 49 deg.",
                          "c) Low Contrast at 25 deg.", "d) Overlap 20-25 deg.","a) No Breaks","f) Temporal Break Year 50")

plist0 <- list(); idx0 <- 1
for(i in 1:length(scens)){
  plist  <- list(); idx <- 1
  ## plot raw data
  scen <- scens[i]
  tempdf <- read.csv(paste0("./IBM_output/datasets/",scen,"_4.csv")) %>% filter(Age == 6)

  breaksdf <- read.csv( paste0("./GAM_output/ldf_raw_a",age,".csv")) %>% filter(scen == scens[i] & boot == 4)
  load(file = paste0("./GAM_output/mod_",scen,"_4.rds")) ## mod for boot 4
  load(file = paste0("./GAM_output/m2d_",scen,"_4.rds")) ## m2.d for boot 4 (contains derivative info)
  
  Terms <- c("Year","Latitude_dd","Longitude_dd")
  
  ## raw data save separately ----
  if(i < 6){ ## plot tempvar separately
    
  plist0[[idx0]]  <- ggplot(tempdf, aes(x = Latitude_dd, y = Longitude_dd)) +
    theme_classic() +
    theme(legend.position = 'none') +
    scale_y_continuous(limits = c(-10,60), breaks = seq(0,50,10)) +
    geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
    scale_fill_viridis_c(guide = "legend") +
    scale_size_continuous(range = c(1, 15)) +
    labs(x = 'Longitude',y = 'Latitude',
         fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
         title = paste0(scens.title0[i], ' Raw Data')) +
    geom_label(x = 10, y = 50, label = paste0('n = ',nrow(tempdf)), 
               size = 4, fill = 'white')
  idx0 <- idx0+1
  }
  
  if(i == 6){ ## plot tempvar separately
    plist0[[idx0]]  <- ggplot(tempdf, aes(x = Year, y = Length_cm)) +
      theme_classic() +
      theme(legend.position = 'right') +
      scale_y_continuous(limits = c(125,450)) +
      geom_point(aes(size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
      scale_fill_viridis_c(guide = "legend") +
      scale_size_continuous(range = c(1, 15)) +
      labs(x = 'Year', y = 'Length of Age Six Fish (cm)',
           fill = "Length of Age-6 Fish (cm)",size  = "Length of Age-6 Fish (cm)",
           title = paste0(scens.title0[i], ' Raw Data')) +
      geom_label(x = 10, y = 400, label = paste0('n = ',nrow(tempdf)), 
                 size = 4, fill = 'white')
    idx0 <- idx0+1
  }
    ggsave(plot = last_plot(),  file = paste0("./figures/rawdat_",scen,".png"),
  width = 5, height = 5, units = 'in', dpi = 480)


  ## map showing new breaks
  plist[[idx]] <- ggplot(tempdf, aes(x = Latitude_dd, y = Longitude_dd)) +
    theme_classic() +
    theme(legend.position = 'right') +
    scale_y_continuous(limits = c(0,50), breaks = seq(0,50,10)) +
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
  for(t in 1:length(Terms)){
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
    m2.dtemp <- data.frame(cbind(m2.d$eval[,Terms[t]], m2.d[[Terms[t]]]$deriv,
                                 CI[[1]]$upper, CI[[1]]$lower)); names(m2.dtemp) = c('x','deriv','upper','lower')
    ## plot deriv
    plist[[idx]] <-  ggplot(m2.dtemp,aes(x = x, y = deriv))    +
      theme_minimal() +
      theme(panel.grid = element_blank())+
      geom_line(lwd = 1.1) +
      geom_hline(yintercept = 0, col = 'grey22') +
      geom_line(aes(y= upper), linetype = 'dashed') +
      geom_line(aes(y= lower), linetype = 'dashed') +
      geom_vline(xintercept = breaksdf[,c('yr_breaks','lat_breaks','lon_breaks')[t]], 
                 col = 'red', lwd = 1.1, linetype = 'dashed') +
      labs(x = Terms[t], y = "f'(x)", title = paste0(letters[idx-1],") ",'First Derivative for ', Terms[t]))
    idx <- idx+1
  }


  lay <- rbind(c(2,2,1,1,1),
               c(3,3,1,1,1),
               c(4,4,1,1,1),
               c(5,5,1,1,1),
               c(6,6,1,1,1),
               c(7,7,1,1,1))
  grid.arrange(grobs = plist, layout_matrix = lay) %>%
  ggsave(plot = .,  file = paste0("./figures/analysis_",scen,".png"), width = 11, height = 8, units = 'in', dpi = 480)

}

ggarrange(plotlist = plist0, ncol=2, nrow=3, common.legend = TRUE, legend="bottom") %>%
  ggsave(plot = .,  file = paste0("./figures/rawdata_compile.png"), width = 11, height = 8, units = 'in', dpi = 480)


a6df <- list.files("./IBM_output/datasets/", full.names = T) %>%
  lapply(read.csv) %>%  
  bind_rows()  %>%
    filter(Age == 6)

dim(a6df)/(100*length(unique(scenarios$DESC))) ## 100 datasets times 5 simulations -- getting average per ds

## Plot proportion agreegments ----
cdfprop <- read.csv(paste0('./gam_output/cdf_prop_',Sys.Date(),'.csv'))
levels(cdfprop$scen) <- c("Break at 25 deg.", "Break at 49 deg.",
                                                  "Low Contrast at 25 deg.", 
                                  "Overlap 20-25 deg.","No Breaks",
                          "Temporal Break at Year 50")
# cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
levels(cdfprop$variable) <- c('Both L1 and L2','L1','L2' )
cdfprop2 <- cdfprop# %>% filter(!(variable %in% c('Both L1 and L2')))

# cdfprop$variable <- factor(cdfprop$variable, levels=c('L1','L2','Both L1 and L2'))
plist1 <- list()
plist1[[1]] <- ggplot(cdfprop2, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = c(0.9,0.75)) +
  scale_fill_viridis_d()+
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Coverage Probability', fill = 'Scenario', 
       title = 'a) Coverage Probability for Endpoints of Growth Curve') +
  geom_bar(stat = 'identity',width=0.6, position = position_dodge(width=0.7)) +
  facet_wrap(~variable)
# ggsave(plot = last_plot(),  file = paste0("./figures/cdfprop.png"), width = 9, height = 6, units = 'in', dpi = 480)

cdfaccu <- read.csv(paste0('./gam_output/cdf_accu_',Sys.Date(),'.csv'))
levels(cdfaccu$scen) <- c("Break at 25 deg.", "Break at 49 deg.",
                          "Low Contrast at 25 deg.", 
                          "Overlap 20-25 deg.","No Breaks","Temporal Break at Year 50")
# cdfprop$scen  <- factor(cdfprop$scen , levels = cdfprop$scen [order(cdfprop$prop )])
levels(cdfaccu$variable) <- c('Lat, Long and Year','Both Latitude and Longitude','Latitude', 'Longitude' ,'Year')
# cdfaccu$variable <- factor(cdfaccu$variable, levels=c('L1','L2','Both L1 and L2'))

cdfaccu2 <- cdfaccu# %>% filter(!(variable %in% c('Lat, Long and Year','Both Latitude and Longitude')))
# cdfaccu$scen  <- factor(cdfaccu$scen , levels = cdfaccu$scen [order(cdfprop$prop  )])
plist1[[2]] <- ggplot(cdfaccu2, aes(x = scen, y = prop, fill = scen)) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = 'none') +
  scale_fill_viridis_d()+
  scale_y_continuous(limits = c(0,1)) +
  labs(x = '',y = 'Proportion Detected Accurate Spatial Breaks', fill = 'Scenario', title = 'b) Proportion Detected Accurate Breaks') +
  geom_bar(stat = 'identity',width=0.5, position = position_dodge(width=0.5)) +
  facet_wrap(~variable, nrow = 2)

ggarrange(plotlist = plist1, ncol=1, nrow=2, common.legend = TRUE, legend="bottom") %>%
  ggsave(plot = .,  file = paste0("./figures/cdfprob_",Sys.Date(),".png"), width = 11, height = 8, units = 'in', dpi = 480)


rbind(cdfprop, cdfaccu) %>% 
  group_by(scen) %>%
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


## combo SAB map for a6 and a10 breaks
  dat <- all_data %>% filter(Age == 30 & Sex == 'M')
    ggplot() + 
    geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
    coord_quickmap() +
    scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
    scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          axis.title =element_blank(), 
          legend.position = c(0.1,0.2),
          legend.key.size = unit(0.75, "cm"),
          legend.text = element_text(size = 12)) +
    geom_hline(yintercept = c(36,50),lwd = 1.1, linetype = 'dashed', col = 'red') +
    geom_vline(xintercept =  -130, lwd = 1.1, linetype = 'dashed', col = 'red') +
      ## ecosystem break
    geom_vline(xintercept =  -145, lwd = 1.1, linetype = 'dashed', col = 'blue') +
    geom_rect(fill = 'white', aes(xmin = -180, xmax = -130.1, ymin = 30, ymax = 49.9)) +
    geom_point(data = dat, aes(x = Longitude_dd, y = Latitude_dd, size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
    scale_fill_viridis_c(guide = "legend") +
    labs(fill = paste0("Length (cm)"),
         size = paste0("Length (cm)"),
         title = "GAM- and Ecosystem-based Regions with Age 30 Fish") +
      geom_label(aes(x = c(-155,-138,-120,-120,-120), 
                    y = c(65,65,65,40,32),
                    label = c(paste('Region ',5:1)) ),size = 6, col = 'black',fill = 'white',show.legend   =FALSE)
    ggsave(plot = last_plot(),  file = paste0("./figures/sab_zones.png"), width = 8, height = 6, units = 'in', dpi = 480)
    
df1<-all_data %>% filter(Age %in% c(4,6,30)) %>% 
  group_by(REG, Age) %>% summarise(n = n()) 

df1 %>%
  group_by(Age) %>% summarise(mn = mean(n))

## predicts and parest for sab ----
## parest: see if actually different. PHASE 1 = pre-merge, PHASE 2- temporal merges
for(phase in c("phase1","phase2")){
  parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),"_",phase,'.csv')) %>% 
    filter(variable ==  "Linf") %>% 
    mutate(source = 'Estimated') %>%
    mutate(REG2 = gsub("_.*", "\\1", REG),
           REG3 = as.numeric(substr(REG,2,2)),
           Sex = gsub(".*_", "\\1", REG),
           lwr = value - 1.96 * sd,
           upr = value + 1.96 * sd,
           matchcol = 'black')
  # parest$value <- exp(parest$value - (parest$sd^2)/2) ## bias correction
  parest <- parest[order(parest$REG2,parest$Sex),]
  
  # parest <- rbind(parest, read.csv("./input_data/true_sab_vals.csv")) 
  # levels(parest$REG) <- c("ALL","R1","AK","R2","BC","R1","WC")
  # parest$REG <- factor(parest$REG ,levels=c("ALL","R1","WC","R2","BC","R3","AK",'R4'))
  
  ## check if CIs overlap
  if(phase == 'phase1'){
    for(i in seq(1,nrow(parest),2)) { ## will go by sex
      parest$match[i:(i+1)] <- ifelse(c(parest$lwr[i], parest$upr[i]) %overlaps%  c(parest$lwr[i+1], parest$upr[i+1]),
                                      paste0('OVERLAP ', parest$REG2[i+1]), 'NO OVERLAP')
    }
    
  } else if(phase == 'phase2') {
    for(i in 1:nrow(parest)){ ## we want to compare w adjacent regions
      if(parest$REG3 < 5){
        tmp <- subset(parest, Sex == parest$Sex[i] & REG3 == parest$REG3[i]+1)
        parest$match[i] <-ifelse(any(c(tmp$lwr, tmp$upr) %overlaps% c(parest$lwr[i], parest$upr[i])),paste0('OVERLAP ', tmp$REG2[1]),'NO OVERLAP')
      } ## end R < 5
    } ## end parest rows
  } ## end else phase2
  
  parest$matchcol <- parest$match != 'NO OVERLAP'
  ggplot(parest, aes(x = REG, y = value, col = matchcol))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = c(0.9,0.1),
          legend.background = element_blank(),
          axis.text = element_text(size = 10,angle = 45),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,100)) +
    scale_color_manual(values = c('black','red')) +
    geom_point() +
    geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    labs(x = 'Spatiotemporal x Sex Stratum', y = "", 
         col = ifelse(phase == 'phase1', 'CI overlap within Region + Sex','CI Overlap Adjacent Region x Sex'),
  title = paste0(phase," Linf Estimates")) +
  facet_wrap(~variable)

ggsave(plot = last_plot(),  file = paste0("./figures/sab_parest_",phase,".png"), width = 10, height = 8, units = 'in', dpi = 480)
write.csv(parest, file = paste0("./GAM_output/overlap_",Sys.Date(),"_",phase,".csv"),row.names=F)

## fits
ypreds <- read.csv(paste0("./GAM_output/SAB_predicts_",Sys.Date(),"_",phase,".csv"))

ypreds$gamREG <- paste0('Region ',ypreds$gamREG)
levels(ypreds$Sex) <- c('Females','Males')
levels(ypreds$Period) <- c('pre-2010','2010-Present','All Years')
for(i in 1:nrow(ypreds)){
  ypreds$Period[i] <-  ifelse(length(grep('pool', ypreds$cREG[i])) == 1, 
                              'All Years', paste(ypreds$Period[i]))
}

ggplot(ypreds, aes(x = Age, y = Predicted, col = REG )) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size=10))+
  scale_alpha(guide = 'none') +
  scale_y_continuous(limits = c(0,110)) +
  scale_color_brewer(palette =  'Accent')+
  geom_point(alpha = 0.5, aes(y = Length_cm)) +
  geom_line(lwd = 1.1, col = 'black')+
  labs(y = 'Length (cm)', col = "Actual Data Source") +
  facet_wrap(~gamREG + Sex + Period, ncol = 4)
ggsave(plot = last_plot(),  
       file = paste0("./figures/sab_fits_",phase,".png"), 
       width = 10, height = 12, units = 'in', dpi = 520)
} ## end phase
## "no discernable pattern in spurious breaks"
cdf %>% filter(scen %in% c('NoBreaks','tempvar_R1R2') & (LAT == FALSE  | LON == FALSE)) %>% group_by(scen, gamLAT,gamLON) %>% dplyr::summarise(n = n())
## what proportion of wrong ones are only one off?
cdf %>% filter(!(scen %in% c('NoBreaks','tempvar_R1R2')) & (LAT == FALSE  | LON == FALSE)) %>% 
  group_by(scen) %>% dplyr::summarise(n = n())

cdf %>% filter(scen %in% c('F0L1S_25','F0L1S_R3') & (LAT == FALSE  | LON == FALSE)) %>% 
  group_by(scen) %>% dplyr::summarise(n = n(), inbounds = sum(gamLAT %in% c(24,26)   |gamLON %in% c(24,26)  )) %>% mutate(inbounds/n) 

cdf %>% filter(YEAR == FALSE) %>% group_by(scen, gamYR) %>% dplyr::summarise(n = n())
cdf %>% filter(scen == 'F0LMW'& LAT == FALSE) %>% group_by(gamLAT) %>% dplyr::summarise(n = n())
cdf %>% group_by(YEAR) %>% dplyr::summarise(n = n()) 


   