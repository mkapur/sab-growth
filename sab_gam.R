## Updated GAM method for SAB
## kapurm spring 2019
## source from gist
rm(list = ls())
compname <- c("Maia Kapur","mkapur")[1]

source(paste0(getwd(),"/functions/Deriv.R"))
compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))
usa <- map_data("world") 
source("./functions/getGR.R")


compname <- c("Maia Kapur","mkapur")[1]
require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2);library(TMB)
require(maps);require(mapdata); require(gridExtra)
## load data
load(paste0("./input_data/gam_data_sab_0315.rda")) ## all_data -- made using gam_dataprep and 15k subsample
all_data$Longitude_dd[all_data$Longitude_dd > 0] <- all_data$Longitude_dd[all_data$Longitude_dd > 0]*-1

for(AA in c(3,10)){
  for(SS in c("F","M")){
    dat <- all_data %>% filter( Age == AA  & Sex == SS) 
    mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd),data = dat)
    png(paste0("./figures/sab_gam_diagnostics_",AA,SS,".png"), width = 7, height = 5, units = 'in', res = 420)
    layout(matrix(1:4, ncol = 2))
    gam.check(mod)
    graphics.off()
    
    ## get & eval derivatives ----
    pdat <- dat
    pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
    p2 <- predict(mod, newdata = pdat) ## raw predicts
    pdat <- transform(pdat, predLen = p2, se2_lat = pTerm$se.fit[,2],se2_lon = pTerm$se.fit[,3], se2_yr = pTerm$se.fit[,1])

    df.res <- df.residual(mod)
    crit.t <- qt(0.025, df.res, lower.tail = FALSE)
    ## variances are additive just FYI
    pdat <- transform(pdat,
                      upper = predLen + (crit.t * (se2_lat+se2_lon+se2_yr)),
                      lower = predLen - (crit.t * (se2_lat+se2_lon+se2_yr)))
    
    breaksdf <- list()
    Terms <- c("Year","Latitude_dd","Longitude_dd")
    for(t in 1:length(Terms)){
      Term <- Terms[t]
      
      
      newD <- data.frame( Year = seq(1981,2017,length = 100),
                          Longitude_dd = seq(-186,-116, length = 100),
                          Latitude_dd = seq(0,64,length = 100)) ## whatever isn't modeled will just get ignored
      
      
      m2.d <- Deriv(mod, newdata = newD)
      m2.dci <- confint(m2.d, term = Term)
      
      crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
      crit.eval.se = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$se.deriv) ## use tails
      
      
      ## identify where mean crosses zero or falls out of bounds (NAs where it does)
      m2.dsig.zeros <- signifMK(x = m2.d$eval[[Term]],
                                d = m2.d[[Term]]$deriv,
                                upper = m2.dci[[Term]]$upper,
                                lower = m2.dci[[Term]]$lower,
                                crit.eval = crit.eval,
                                eval = 0)
      
      pix <- !is.na(m2.dsig.zeros)
      
      vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
      breaksdf[[t]] <- sort(c(unique(round(vals)))) ## get rounded unique
      ## fill NAs in bdf for binding
      for(i in 1:length(breaksdf)){
        if (length(breaksdf[[i]]) == 0){ ## fill NA for empty
          breaksdf[[i]] <- NA
        }
      }
      
      
    }

    plist <- list()
    idx <- 2
    ## extract plotting stuff
    pd <- plot(mod,   select = t, scheme  =2, lwd  =2, main = paste0(Terms[t],' Smoother'), cex.axis = 2, ylim = c(-10,ifelse(t != 3,10,500)))
    for(t in 1:length(Terms)){
      temp0 <- pd[t] ## get this smoother
      temp <- data.frame(cbind(temp0[[1]]$x,temp0[[1]]$fit, temp0[[1]]$se)); names(temp) = c('x','fit','se')
      plist[[idx]] <-  ggplot(temp, aes(x = x, y = fit)) +
        theme_minimal() +
        theme(panel.grid = element_blank())+
        geom_line(lwd = 1.1) +
        geom_line(aes(y= fit-se), linetype = 'dashed') +
        geom_line(aes(y= fit+se), linetype = 'dashed') +
        geom_rug(sides = 'b') +
        labs(x = Terms[t], y = "smoother", title = paste0(letters[idx-1],") ",'Smoother for ', temp0[[1]]$xlab)) 
      idx <- idx+1
      CI <- confint(m2.d, term = Terms[t])
      m2.dtemp <- data.frame(cbind(m2.d$eval[,Terms[t]], m2.d[[Terms[t]]]$deriv, CI[[1]]$upper, CI[[1]]$lower)); names(m2.dtemp) = c('x','deriv','upper','lower')
      plist[[idx]] <-  ggplot(m2.dtemp,aes(x = x, y = deriv))    +
        theme_minimal() +
        theme(panel.grid = element_blank())+
        geom_line(lwd = 1.1) +
        geom_hline(yintercept = 0, col = 'grey22') +
        geom_line(aes(y= upper), linetype = 'dashed') +
        geom_line(aes(y= lower), linetype = 'dashed') +
        geom_vline(xintercept = breaksdf[[t]], lwd = 1.1, linetype = 'dashed', col = 'red') +
        labs(x = Terms[t], y = "f'(x)", title =  paste0(letters[idx-1],") ",'First Derivative for ', temp0[[1]]$xlab)) 
      idx <- idx+1
    }
    
    
    plist[[1]] <- ggplot() + 
      geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
      coord_quickmap() +
      scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
      scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            axis.title =element_blank(), 
            legend.position = c(0.2,0.15)) +
      geom_hline(yintercept = breaksdf[[2]],lwd = 1.1, linetype = 'dashed', col = 'red')+
      geom_vline(xintercept = breaksdf[[3]], lwd = 1.1, linetype = 'dashed', col = 'red') +
      geom_point(data = dat, aes(x = Longitude_dd, y = Latitude_dd, size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
      scale_fill_viridis_c(guide = "legend") +
      labs(fill = paste0("Length of Age-",AA," ",SS," Fish (cm)"),
           size = paste0("Length of Age-",AA," ",SS," Fish (cm)")) +
      ggtitle(paste0('g) GAM-Estimated Regions, with raw data Age-',AA,", ",SS)) +
      geom_label(aes(x = c(-155,-125,-125), y = c(60,45,60), 
                     label = paste0("GAM Estimated Region ",c(3,1,2))),
                 fill = "white",col = 'black')
    
    lay <- rbind(c(2,2,3,3,NA),
                 c(4,4,5,5,NA),
                 c(6,6,1,1,1),
                 c(7,7,1,1,1))
    grid.arrange(grobs = plist, layout_matrix = lay)  %>% 
      ggsave(plot = .,  file = paste0("./figures/sab_smooth_map_",AA,SS,".png"), width = 12, height = 15, units = 'in', dpi = 480)
    
    
    
    ## now re-aggregate and do the fits on all data, based on model ----
    DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat)) 

    # for (i in 1:nrow(all_data)) {
    #   all_data[i,'gamREG'] <-
    #     ifelse(all_data[i, "Latitude_dd"] <=  breaksdf[[2]] & all_data[i, "Longitude_dd"] > breaksdf[[3]], 1,
    #            ifelse(all_data[i, "Latitude_dd"] > breaksdf[[2]] &  all_data[i, "Longitude_dd"] < breaksdf[[3]], 3, 
    #                   ifelse(all_data[i, "Latitude_dd"] <= breaksdf[[2]] &  all_data[i, "Longitude_dd"] < breaksdf[[3]], 4, 2)))
    #   
    # }
    breaksdf <- data.frame(lat_breaks2 = breaksdf[[2]], lon_breaks2 = breaksdf[[1]])
    dat <- getGR(tempdf = all_data, breaksdf)
    ## sanity check
    dat %>% group_by(gamREG) %>% summarise(mnlat = mean(Latitude_dd), mnlon = mean(Longitude_dd))
    
    DES <- ifelse(!is.na(dat$gamREG), as.numeric(dat$gamREG),1)-1 ## this is now numeric index, R3=slot 3 (idx 2)
    KEY <- paste("sab",DES,sep = "_")
    keybase <- paste0("R",1:length(unique(DES))) ## text regions
    
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
        log_Linf = rep(log(70), nStrata),
        log_k = rep(log(0.5), nStrata),
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
    dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
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
    
    rep0$value[rep0$variable %in% c('log_k','log_Linf')] <- exp(rep0$value[rep0$variable %in% c('log_k','log_Linf')])
    rep0$sd[rep0$variable == 'log_Linf'] <- exp(rep0$sd[rep0$variable == 'log_Linf'])
    rep0$variable <- factor(rep0$variable, levels = c("k","Linf","Sigma","t0","log_k","log_Linf")) ## enable new levels
    rep0$variable[rep0$variable == 'log_k'] <- 'k'
    rep0$variable[rep0$variable == 'log_Linf'] <- 'Linf'
    
    
    write.csv(rep0, file = paste0("./GAM_output/SAB_parEst_gam_",AA,SS,"_",Sys.Date(),'.csv'),row.names = F)
    
    ypreds0 <- cbind(dat0,dat) %>% data.frame()  
    names(ypreds0)[1] <- c('Predicted')
    
    write.csv(ypreds0,  paste0("./GAM_output/SAB_predicts",AA,SS,"_",Sys.Date(),".csv"),row.names = F)
    
    cat("Fit TMB model  & saved outputs \n")
    ## plotting ----
    
    
    
    ## plot estimates
    # parest <- read.csv(paste0("./GAM_output/SAB_parEst_gam_",AA,SS,"_",Sys.Date(),".csv")) %>%
    #   filter(variable != "Sigma" & variable != 't0') %>% mutate(source = 'Estimated')
    # 
    # ## bind only regions of use
    # parest <- rbind(parest, read.csv("./input_data/true_sab_vals.csv")) 
    # 
    # ## exponentiate logk
    # parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
    # parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))
    # # levels(parest$REG) <- c("ALL","R1","AK","R2","BC","R1","WC")
    # parest$REG <- factor(parest$REG ,levels=c("ALL","R1","WC","R2","BC","R3","AK",'R4'))
    # cat("Read in Parameter Estimates \n")
    # 
    # plist <- list()
    # plist[[1]] <- ggplot(parest, aes(x = REG, y = value, col = source))+
    #   theme_bw() +
    #   theme(panel.grid = element_blank(),
    #         legend.position = c(0.9,0.9),
    #         legend.background = element_blank(),
    #         axis.text = element_text(size = 10),
    #         axis.title = element_text(size = 10),
    #         legend.text = element_text(size = 10),
    #         strip.text = element_text(size=14))+
    #   scale_color_manual(values = c("red","black"))+
    #   geom_point() +
    #   geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    #   labs(x = 'Spatial Stratum', y = "", col = "") +
    #   facet_wrap(~variable, scales = "free_y") +
    #   ggtitle("a)")
    # 
    # ## plot fits
    # ypreds <- read.csv(  paste0("./GAM_output/SAB_predicts",AA,SS,"_",Sys.Date(),".csv"))
    # ypreds$gamREG <- paste0('GAM-defined Region ',ypreds$gamREG," ", ypreds$Sex)
    # ypreds$Sex <- SS
    # 
    # ## fits
    # plist[[2]] <- ggplot(ypreds, aes(x = Age, y = Predicted, col = REG )) +
    #   theme_bw() +
    #   theme(panel.grid = element_blank(),
    #         legend.position = c(0.9,0.15),
    #         legend.background = element_blank(),
    #         axis.text = element_text(size = 14),
    #         axis.title = element_text(size = 14),
    #         legend.text = element_text(size = 14),
    #         strip.text = element_text(size=14))+
    #   scale_alpha(guide = 'none') +
    #   geom_point(alpha = 0.2, aes(y = Length_cm)) +
    #   geom_line(lwd = 1.1, col = 'black')+
    #   labs(y = 'Length (cm)', col = "Actual Data Source", title = paste0("b) Predicted Length-at-Age ",AA,", ",SS)) +
    #   facet_wrap(~gamREG,ncol = 3)
    # 
    # 
    # lay <- rbind(c(1,1,1,1),
    #              c(1,1,1,1),
    #              c(2,2,2,2),
    #              c(2,2,2,2))
    # 
    # grid.arrange(grobs = plist, layout_matrix = lay) %>%  
    #   ggsave(plot = .,  file = paste0("./figures/sab_fits_a",AA,SS,".png"), width = 10, height = 12, units = 'in', dpi = 480)
    # cat("Plotted overall fits \n")
    
  }
}
