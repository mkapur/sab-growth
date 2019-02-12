## Code to test GAM method on datasets generated via iPopSim (IBM)
## This will fit a gam for each replicate, identify break points, fit new life history params
## and plot comparisons
## Maia Sosa Kapur kapurm@uw.edu WiSp 2019

rm(list = ls())

require(mgcv);require(dplyr);require(ggplot2); require(TMB)
library(gridExtra); library(grid); library(lattice)

compname <- c("Maia Kapur","mkapur")[1]
# setwd(paste0("C:/Users/",compname,"/Dropbox/UW/sab-growth/gam"))

## Build mods, get breakpoints

fLevs <- read.csv(paste0("C:/Users/",compname,"/Dropbox/UW/sab-growth/iPopSim/inputs/scenarios.csv"),na.strings = 'NA') ## manual file
age <- 6; nboot <- 100; testrows <- unique(fLevs$DESC) ##
source("bootBreaks.R") ## about 10 mins

## read in CSV generated above
# ldf <- read.csv("summary_tables/ldf_raw_age_6.csv")

## assign levels so it plots in order
levels(ldf$scen) <-  c("No Breaks", "Break at 25°", "Break at 30°", 
                       "Overlap 20° - 25°", "Break at 49°")

trueb <- c(factor(NA),25,30,NA,49) ## true breaks, FILL OVERLAP BELOW
ntrue <- ldf %>% group_by(scen) %>% summarise(ct = n()) %>% data.frame() %>% mutate(trueb)
ntrue <- rbind(ntrue, data.frame(scen = rep("Overlap 20° - 25°",6),ct = rep(NA,6), trueb = 20:25))


 ldf %>% 
  group_by(scen,  lat_breaks2) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  ggplot(., aes(x = lat_breaks2, y = freq)) +
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(limits = c(paste(1:55),NA),breaks = c(paste(seq(0,50,5)),NA)) +
  geom_vline(data = ntrue, aes(xintercept =trueb), col = 'red', linetype = 'dashed') +
  facet_wrap(~scen,ncol = 1) + 
  labs(main = 'breaks identified', y = 'frequency', x = 'break location (° latitude)', main = 'spatial breaks')

ggsave(last_plot(), file = "./plots/ldf_a6.png")

## get max values
ldfprop <- read.csv("summary_tables/ldf_raw_age_6.csv") %>%
  group_by(scen,  lat_breaks) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>% 
  group_by(scen) %>%
  filter(freq == max(freq)) %>% data.frame()

write.csv(ldfprop, paste0(getwd(),"/summary_tables/ldf_prop_a6.csv"),row.names = F)

## Now Aggregate each dataset at IDd breaks and assign R1/R2 ----
# gs[[l]] <- breakhist ## storage for plots

# p <- qplot(1,1)
# p2 <- xyplot(1~1)
# r <- rectGrob(gp=gpar(fill="grey90"))
# t <- breakhist
# grid.arrange(t, p)
# 
# grid.arrange(grobs=gs[1:3], ncol=2, widths = 1:2, 
#              heights=unit(c(1,10), c("in", "mm")))

for(l in 1:nrow(ldfprop)){
  
  ## get scenario name
  # scen0 <- paste0(fLevs[l,'DESC'])
  # scen1 <-  ifelse(is.na(fLevs[l, 5]),
  #                  paste(fLevs[l, 3], fLevs[l, 4], sep = "_"),
  #                  paste(fLevs[l, 3], fLevs[l, 4], fLevs[l, 5],  sep = "_"))
  # scen <- paste(scen0,scen1,sep = "_")
  scen <- paste0(ldfprop[l,'scen'])
  
  # for(b in 1:nboot){ ## loop boots
    dat <- read.csv(paste0("C:/users/",compname,"/dropbox/uw/sab-growth/ipopsim/gendata/",
                           scen,"_",b,'.csv'))
    
    
    ## now re-aggregate the data ----
    ## generate DES matrix of vectors and a KEY for later comparison
    DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat)) ## H1 is all zeros
    ## get region
    # temp0 <- with(dat, ifelse(Latitude_dd >= 49 & Year < 2005 ,0,
    #                                ifelse(Latitude_dd >= 49 & Year >= 2005,1,
    #                                       ifelse(Latitude_dd < 49 & Year < 2005,2,3))))
    # ## get region and sex combos as factors
    # DES <- as.numeric(factor(paste0(temp0,all_data[,'Sex'])))-1
    DES <- ifelse(!is.na(dat$REG), as.numeric(dat$REG),1)-1 ## overwrite NA for nobreaks
    KEY <- paste(scen,DES,sep = "_")
    keybase <- unique(ifelse(!is.na(dat$REG), as.character(dat$REG),"R0")) ## text regions
    # KEY <- paste0(c("North_Early","North_Late","South_Early","South_Late")[(temp0+1)],"_",all_data[,'Sex'])
    # keybase <- apply(expand.grid(c("North_Early","North_Late","South_Early","South_Late"), c('F','M')), 1,
    #                  paste, collapse="_")
    
    # save(DES, file = paste0(getwd(),"/data/DES_gam.rda"))
    # save(KEY, file = paste0(getwd(),"/data/KEY_gam.rda"))
    # load( paste0(getwd(),"/data/DES_gam.rda")) ## DES
    # load( paste0(getwd(),"/data/KEY_gam.rda")) ## KEY
    
    ## run TMB parest----
    compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
    dyn.load(dynlib("sptlvb"))
    
    dat0 <- rep0 <- NULL ## later storage
    nStrata <- length(unique(DES))
    
    data <-
      list(
        Length_cm = dat[,"Length_cm"],
        Age = dat[,"Age"],
        DES = as.vector(DES),
        nStrata = nStrata
      )
    
    parameters <-
      list(
        log_Linf = rep(log(200), nStrata),
        log_k = rep(log(0.1), nStrata),
        t0 = rep(0.1, nStrata),
        log_Sigma = 0.1
      )
    
    # Now estimate everything
    map <- NULL
    model <- MakeADFun(data, parameters,  DLL="sptlvb",silent=T,map=map)
    fit <- nlminb(
      model$par,
      model$fn,
      model$gr,
      control = list(
        rel.tol = 1e-12,
        eval.max = 100000,
        iter.max = 10000
      )
    )
    # for (k in 1:3)  fit <- nlminb(model$env$last.par.best, model$fn, model$gr) ## start at last-best call, for stability
    
    best <- model$env$last.par.best
    
    rep <- sdreport(model)
    dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) ## each 6 cols is new sim
    rep0 <- bind_rows(rep0,
                      bind_cols(
                        data.frame(names(rep$value)),
                        data.frame(rep$value),
                        data.frame(rep$sd),
                        # data.frame(c(rep(unique(
                        #   KEY
                        # ), 3), rep("ALL", 1)))
                        data.frame(c(rep(keybase
                                         , 3), rep("ALL", 1)))
                        
                      ))
    
    ## reformat outputs ----
    names(rep0) <- c('variable', 'value','sd', 'REG')
    # rep0 <- rep0 %>% mutate(
    #   st =sub('_([^_]*)$', '',  ID),
    #   Sex = sub(".*_ *(._?)", "\\1", ID)
    # ) %>% select(-ID)
    
    write.csv(rep0, file = paste0(getwd(),"/results/",scen,"_parEst_gam_",Sys.Date(),'.csv'),row.names = F)
    
    ypreds0 <- cbind(dat0,dat) %>% data.frame()  
    names(ypreds0)[1] <- c('Predicted')
    ypreds0$REG <- as.factor(sub('_([^_]*)$', '', KEY))
    
    write.csv(ypreds0,  paste0(getwd(),"/results/",scen,"_predicts",Sys.Date(),".csv"),row.names = F)
    ypreds <- read.csv( paste0(getwd(),"/results/",scen,"_predicts",Sys.Date(),".csv"))
    ## plotting ----
    
    ## fits
    ggplot(ypreds, aes(x = Age, y = Predicted )) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            legend.position = c(0.9,0.1),
            axis.text = element_text(size = 18),
            legend.text = element_text(size = 14),
            strip.text = element_text(size=14))+
      scale_y_continuous(limits = c(0,100)) +
      # scale_x_continuous(limits = c(0,50)) +
      scale_color_manual(values = c("#d8b365","#5ab4ac"))+
      scale_alpha(guide = 'none') +
      geom_point(alpha = 0.2, aes(y = Length_cm)) +
      # geom_line(lwd = 1.1) +
      labs(y = 'Length (cm)', col = "") +
      facet_wrap(~REG)
    # ggsave(file = paste0(getwd(),"/plots/fits_gam",Sys.Date(),".png"), 
    #        plot = last_plot(), height = 8, width = 10, unit = 'in', dpi = 520)
    
    
    ## parest by sex and region
    parest <- read.csv(paste0(getwd(),"/results/",scen,"_parEst_gam_",Sys.Date(),'.csv')) %>%
      filter(variable != "Sigma" & variable != 't0') %>% mutate(source = 'Estimated')
    parest <- rbind(parest, read.csv("true_vals.csv"))
    ## exponentiate logk
    parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
    parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))
    
    # levels(ldf$scen) <-  c("No Breaks", "Break at 25°", "Break at 30°", 
    #                        "Overlap 20° - 25°", "Break at 49°")
    # 
    ldf %>% 
      group_by(scen,  lat_breaks2) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n)) %>%
      filter(scen == scen) %>%
      ggplot(., aes(x = lat_breaks2, y = freq)) +
      geom_bar(stat = 'identity') + 
      theme_bw() + 
      theme(panel.grid = element_blank())+
      scale_x_discrete(limits = c(paste(1:55),NA),breaks = c(paste(seq(0,50,5)),NA)) +
      geom_vline(data = ntrue, aes(xintercept =trueb), col = 'red', linetype = 'dashed') +
      # facet_wrap(~scen,ncol = 1) + 
      labs( y = 'frequency', x = 'break location (° latitude)', title = scen)
    
    
    plist[[l]] <- ggplot(parest, aes(x = REG, y = value, col = source))+
      theme_bw() +
      theme(panel.grid = element_blank(), 
            # legend.position = c(0.25,0.9),
            legend.text = element_text(size = 14),
            strip.text = element_text(size=14),
            axis.text = element_text(angle = 45))+
      scale_color_manual(values = c("red","black"))+
      geom_point() +
      geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
      labs(x = 'Spatial Stratum', y = "", col = "") +
      facet_wrap(~variable, scales = "free_y") 
    
}
    
    ggsave(file = paste0(getwd(),"/plots/",scen,"_parplot_gam",Sys.Date(),".png"), 
           plot = last_plot(), height = 6, width = 8, unit = 'in', dpi = 520)
    
    

## dep----
# ydf$year_breaks2 <- factor(ydf$year_breaks, levels=c(paste(1:50),NA))
# 
# ydf %>% 
#   group_by(scen,  year_breaks2) %>% 
#   summarise(n = n()) %>% 
#   mutate(freq = n/sum(n)) %>% 
#   ggplot(., aes(x = year_breaks2, y = freq)) +
#   geom_bar(stat = 'identity') + 
#   theme_bw() + 
#   theme(panel.grid = element_blank())+
#   scale_x_discrete(limits = c(paste(1:50),NA),breaks = c(paste(seq(1,50,4)),NA)) +
#   facet_wrap(~scen) + 
#   labs(main = 'breaks identified', x = 'break location (year)', main = 'temporal breaks')
# ggsave(last_plot(), file = "./plots/ydf_a6.png")


# ydf %>% 
#   group_by(scen,  year_breaks) %>% 
#   summarise(n = n()) %>% 
#   mutate(freq = n/sum(n)) %>%
#   write.csv(paste0(getwd(),"/summary_tables/ydf_prop.csv"),row.names = F)
