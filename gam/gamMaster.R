## Code to test GAM method on datasets generated via iPopSim (IBM)
## This will fit a gam for each replicate, identify break points, fit new life history params
## and plot comparisons
## Maia Sosa Kapur kapurm@uw.edu WiSp 2019

rm(list = ls())

require(mgcv);require(dplyr);require(ggplot2); require(TMB)
library(gridExtra); library(grid); library(lattice)

compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))

compname <- c("Maia Kapur","mkapur")[1]

## Build mods, get breakpoints
scenarios <- read.csv(paste0("C:/Users/",compname,"/Dropbox/UW/sab-growth/iPopSim/inputs/scenarios.csv"),na.strings = 'NA') ## manual file
age <- 6; nboot <- 100; testrows <- unique(scenarios$DESC) ##
# source("bootBreaks.R") ## about 10 mins -- don't need to do this >1x

## read in CSV generated above
ldf <- read.csv("summary_tables/ldf_raw_a6.csv")
ldfprop <- read.csv("summary_tables/ldf_prop_a6.csv")
ntrue <- read.csv("summary_tables/ntrue_a6.csv")


## Now Aggregate each dataset at IDd breaks and assign R1/R2 ----
for(l in 1:nrow(ldfprop)){ 
  for(b in 1:nboot){

  scen <- paste0(ldfprop[l,'scen'])
  
  # for(b in 1:nboot){ ## loop boots
  dat <- read.csv(paste0("C:/users/",compname,"/dropbox/uw/sab-growth/ipopsim/gendata/",
                         scen,"_",b,'.csv'))
  
  
  getGR <- function(ldf,dat){
    
    temp <- subset(ldf, scen == ldfprop[l,'scen'] & boot == b)
    for(i in 1:nrow(temp)){ ## loop unique breaks
    sptl <- temp[i,'lat_breaks2'] ## diagnosed break
    for(i in 1:nrow(dat)){  
      if(is.na(sptl)){
        dat[i,"gamREG"] <- "R1"
      } else if (dat[i,"Latitude_dd"] < sptl) {dat[i,"gamREG"] <- "R1"}
      else if(dat[i,"Latitude_dd"] >= sptl & "R3" %in% levels(dat$REG)){
        dat[i,"gamREG"] <- "R3"
      }    else if(dat[i,"Latitude_dd"] >= sptl & "R2" %in% levels(dat$REG)){
        dat[i,"gamREG"] <- "R2"
      }
    }
    ## sanity check
    dat %>% group_by(gamREG) %>% summarise(mean(Latitude_dd))
    
    
  }
  # dat$Length_cm <- dat$Length_cm*10 ## currently data are in in wrong units
  if(scen == 'NoBreaks') dat$REG <- as.factor('R1')
  
  ## now re-aggregate the data ----
  ## generate DES matrix of vectors and a KEY for later comparison
  DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat)) 
  DES <- ifelse(!is.na(dat$REG), as.numeric(dat$REG),1)-1 ## overwrite NA for nobreaks
  KEY <- paste(scen,DES,sep = "_")
  keybase <- unique(ifelse(!is.na(dat$REG), as.character(dat$REG),"R1")) ## text regions
  keybase[keybase == 'R0'] <- "R1"


  
  ## run TMB parest----

  
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

                      data.frame(c(rep(keybase
                                       , 3), rep("ALL", 1)))
                      
                    ))
  
  ## reformat outputs ----
  names(rep0) <- c('variable', 'value','sd', 'REG')
  write.csv(rep0, file = paste0(getwd(),"/results/",scen,"_parEst_gam_",Sys.Date(),'.csv'),row.names = F)
  
  ypreds0 <- cbind(dat0,dat) %>% data.frame()  
  names(ypreds0)[1] <- c('Predicted')

  write.csv(ypreds0,  paste0(getwd(),"/results/",scen,"_predicts",Sys.Date(),".csv"),row.names = F)
  # ypreds <- read.csv( paste0(getwd(),"/results/",scen,"_predicts",Sys.Date(),".csv"))
  
  cat("Fit TMB model ",scen," & saved outputs \n")
  ## plotting ----
  
  ## fits
  # ggplot(ypreds, aes(x = Age, y = Predicted )) +
  #   theme_bw() +
  #   theme(panel.grid = element_blank(),
  #         legend.position = c(0.9,0.1),
  #         axis.text = element_text(size = 18),
  #         axis.title = element_text(size = 18),
  #         legend.text = element_text(size = 14),
  #         strip.text = element_text(size=14))+
  #   # scale_y_continuous(limits = c(0,100)) +
  #   scale_x_continuous(limits = c(0,15)) +
  #   scale_color_manual(values = c("#d8b365","#5ab4ac"))+
  #   scale_alpha(guide = 'none') +
  #   geom_point(alpha = 0.2, aes(y = Length_cm)) +
  #   geom_line(lwd = 1.1) +
  #   labs(y = 'Length (cm)', col = "") +
  #   facet_wrap(~REG)
  # ggsave(file = paste0(getwd(),"/plots/",scen,"_fits_gam",Sys.Date(),".png"),
  #        plot = last_plot(), height = 8, width = 10, unit = 'in', dpi = 520)
  
  
  ## parest by sex and region
  # parest <- read.csv(paste0(getwd(),"/results/",scen,"_parEst_gam_",Sys.Date(),'.csv')) %>%
  #   filter(variable != "Sigma" & variable != 't0') %>% mutate(source = 'Estimated')
  # parest <- rbind(parest, read.csv("true_vals.csv"))
  # ## exponentiate logk
  # parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
  # parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))
  # 
  # ggplot(parest, aes(x = REG, y = value, col = source))+
  #   theme_bw() +
  #   theme(panel.grid = element_blank(), 
  #         # legend.position = c(0.25,0.9),
  #         legend.text = element_text(size = 14),
  #         strip.text = element_text(size=14),
  #         axis.text = element_text(angle = 45))+
  #   scale_color_manual(values = c("red","black"))+
  #   geom_point() +
  #   geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
  #   labs(x = 'Spatial Stratum', y = "", col = "") +
  #   facet_wrap(~variable, scales = "free_y") 
  # ggsave(file = paste0(getwd(),"/plots/",scen,"_parplot_gam",Sys.Date(),".png"), 
  #        plot = last_plot(), height = 6, width = 8, unit = 'in', dpi = 520)
  # 
}



## loop plotting ----
plist <- list(); idx <- 1
ldf$lat_breaks2 <- factor(ldf$lat_breaks, levels=c(paste(1:50),NA))

# ldfprop$scentemp <-  c( "Break at 25?", "Break at 30?",  "Break at 49?", "Break at 49?",
#                        "Overlap 20? - 25?","No Breaks")
levels(ldf$scen) <-  c( "Break at 25 deg.", "Break at 49 deg.",
                        "Low Contrast at 25 deg.", "Overlap 20-25 deg.","No Breaks")

ldfprop$scentemp  <-  c( "Break at 25 deg.", "Break at 49 deg.","Break at 49 deg.",
                                              "Low Contrast at 25 deg.", "Overlap 20-25 deg.","No Breaks")
## should match
for(l in 1:length(unique(ldf$scen))){
  
  
  scentemp <- unique(ldf$scen)[l]
  scen <- paste0(ldfprop$scen[ldfprop$scentemp == scentemp]) ## filenames
  
  ## plot hist
  plist[[idx]]  <-  ldf %>%
    filter(scen == scentemp) %>%
    group_by(lat_breaks2) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    ggplot(., aes(x = lat_breaks2, y = freq)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = c(0.9,0.1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = c(paste(1:55),NA),breaks = c(paste(seq(0,50,5)),NA)) +
    facet_wrap(~scen,ncol = 1) +
    geom_vline(data = subset(ntrue, scen == scentemp), aes(xintercept =trueb), col = 'red', linetype = 'dashed') +
    labs(main = 'breaks identified', y = 'frequency', x = 'break location (deg. latitude)', main = 'spatial breaks')
  
  idx <- idx +1
  
  ## plot estimates
  parest <- read.csv(paste0(getwd(),"/results/",scen[1],"_parEst_gam_",Sys.Date(),'.csv')) %>%
    filter(variable != "Sigma" & variable != 't0') %>% mutate(source = 'Estimated')
  
  ## bind only regions of use
  parest <- rbind(parest, read.csv("true_vals.csv") %>% filter(REG %in% parest$REG))
  
  ## exponentiate logk
  parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
  parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))
  levels(parest$REG) <- 
    c('Single Region','Regime 1',
      ifelse(levels(parest$REG)[3] == 'R2', 'Regime 2','Regime 3'),
      ifelse(levels(parest$REG)[4] == 'R3', 'Regime 3','Regime 2'))
      
  
  plist[[idx]] <- ggplot(parest, aes(x = REG, y = value, col = source))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position =  ifelse(l == 1, 'left','none'),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size=14))+
    scale_color_manual(values = c("red","black"))+
    geom_point() +
    geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    labs(x = 'Spatial Stratum', y = "", col = "") +
    facet_wrap(~variable, scales = "free_y") 
  idx <- idx + 1
  ## plot fits
  ypreds <- read.csv( paste0(getwd(),"/results/",scen[1],"_predicts",Sys.Date(),".csv"))
  levels(ypreds$REG) <-  c('Regime 1',ifelse(levels(ypreds$REG)[2] == 'R2', 'Regime 2', 'Regime 3'))
  ## fits
  # mpred <- ypreds %>% group_by(Age,REG) %>% summarise(pred = mean(Predicted))
  plist[[idx]] <- ggplot(ypreds, aes(x = Age, y = Predicted )) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = c(0.9,0.1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          strip.text = element_text(size=14))+
    scale_y_continuous(limits = c(0,500)) +
    scale_x_continuous(limits = c(0,15)) +
    scale_color_manual(values = c("#d8b365","#5ab4ac"))+
    scale_alpha(guide = 'none') +
    geom_point(alpha = 0.2, aes(y = Length_cm), col = 'red') +
    geom_line(lwd = 1.1)+
    # geom_line(data = mpred, aes(y = pred),lwd = 1.1) +
    labs(y = 'Length (cm)', col = "") +
    facet_wrap(~REG)
  idx <- idx + 1
  
}
grid.arrange(grobs = plist, ncol = 3) %>%
ggsave(plot = ., 
       file = paste0(getwd(),"/plots/fullfigure.png"), width = 17, height = 11, units = 'in', dpi = 480)


ldf %>% filter(scen == "Overlap 20-25 deg.") %>% mutate(overl = lat_breaks %in% c(20:25)) %>% group_by(overl) %>% summarise(n = n()) 
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
