## Updated GAM method for SAB
## kapurm spring 2019
## source from gist
rm(list = ls())
require(mgcv);require(dplyr);require(ggplot2); require(TMB); library(reshape2)
library(gridExtra); library(grid); library(lattice)
compname <- c("Maia Kapur","mkapur")[1]

source(paste0(getwd(),"/functions/Deriv.R")); source("./functions/getGR.R")
# compile("sptlvbSel.cpp") ## will throw 'incomplete line' msg, ignore
# dyn.load(dynlib("sptlvbSel"))
compile("it10.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("it10"))

usa <- map_data("world") 


## load data
# load(paste0("./input_data/gam_data_sab_0315.rda")) ## all_data -- made using gam_dataprep and 15k subsample
# all_data$Longitude_dd[all_data$Longitude_dd > 0] <- all_data$Longitude_dd[all_data$Longitude_dd > 0]*-1

load(paste0("./input_data/gam_data_sab_0415.rda")) ## full_data -- made using gam_dataprep NOT 15k subsample
full_data$Longitude_dd[full_data$Longitude_dd > 0] <- full_data$Longitude_dd[full_data$Longitude_dd > 0]*-1

## Fit GAMS to ID breakpoints [don't need to repeat if data hasn't changed] ----
# for(AA in c(4,6,30)){
#   for(SS in c("F","M")){
#     dat <- all_data %>% filter( Age == AA  & Sex == SS)
#     # dat <- all_data %>% filter( Age == 6  | Age == 10  & Sex == SS)
#     
#     mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd),data = dat)
#     png(paste0("./figures/sab_gam_diagnostics_",AA,SS,".png"), width = 7, height = 5, units = 'in', res = 420)
#     layout(matrix(1:4, ncol = 2))
#     gam.check(mod)
#     graphics.off()
#     
#     ## get & eval derivatives ----
#     pdat <- dat
#     pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
#     p2 <- predict(mod, newdata = pdat) ## raw predicts
#     pdat <- transform(pdat, predLen = p2, se2_lat = pTerm$se.fit[,2],se2_lon = pTerm$se.fit[,3], se2_yr = pTerm$se.fit[,1])
#     
#     df.res <- df.residual(mod)
#     crit.t <- qt(0.025, df.res, lower.tail = FALSE)
#     ## variances are additive just FYI
#     pdat <- transform(pdat,
#                       upper = predLen + (crit.t * (se2_lat+se2_lon+se2_yr)),
#                       lower = predLen - (crit.t * (se2_lat+se2_lon+se2_yr)))
#     
#     breaksdf <- list()
#     Terms <- c("Year","Latitude_dd","Longitude_dd")
#     for(t in 1:length(Terms)){
#       Term <- Terms[t]
#       
#       
#       newD <- data.frame( Year = seq(1981,2017,length = 100),
#                           Longitude_dd = seq(-186,-116, length = 100),
#                           Latitude_dd = seq(0,64,length = 100)) ## whatever isn't modeled will just get ignored
#       
#       
#       m2.d <- Deriv(mod, newdata = newD)
#       m2.dci <- confint(m2.d, term = Term)
#       
#       crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
#       crit.eval.se = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$se.deriv) ## use tails
#       
#       
#       ## identify where mean crosses zero or falls out of bounds (NAs where it does)
#       m2.dsig.zeros <- signifMK(x = m2.d$eval[[Term]],
#                                 d = m2.d[[Term]]$deriv,
#                                 upper = m2.dci[[Term]]$upper,
#                                 lower = m2.dci[[Term]]$lower,
#                                 crit.eval = crit.eval,
#                                 eval = 0)
#       
#       pix <- !is.na(m2.dsig.zeros)
#       
#       vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
#       breaksdf[[t]] <- sort(c(unique(round(vals)))) ## get rounded unique
#       ## fill NAs in bdf for binding
#       for(i in 1:length(breaksdf)){
#         if (length(breaksdf[[i]]) == 0){ ## fill NA for empty
#           breaksdf[[i]] <- NA
#         }
#       }## end breaksdf
#       
#       
#     } ## end terms
#     
#     plist <- list()
#     idx <- 2
#     ## extract plotting stuff
#     pd <- plot(mod,   select = t, scheme  =2, lwd  =2, main = paste0(Terms[t],' Smoother'), cex.axis = 2, ylim = c(-10,ifelse(t != 3,10,500)))
#     for(t in 1:length(Terms)){
#       
#       ## plot smooth
#       temp0 <- pd[t] ## get this smoother
#       temp <- data.frame(cbind(temp0[[1]]$x,temp0[[1]]$fit, temp0[[1]]$se)); names(temp) = c('x','fit','se')
#       plist[[idx]] <-  ggplot(temp, aes(x = x, y = fit)) +
#         theme_minimal() +
#         theme(panel.grid = element_blank())+
#         geom_vline(xintercept = ifelse(is.na(breaksdf[[t]]),-999,breaksdf[[t]]), lwd = 1.1, linetype = 'dashed', col = 'red') +
#         scale_x_continuous(expand = c(0,0), limits = c(ifelse(t == 2, 30, min(temp$x)),max(temp$x))) +
#         geom_line(lwd = 1.1) +
#         geom_line(aes(y= fit-se), linetype = 'dashed') +
#         geom_line(aes(y= fit+se), linetype = 'dashed') +
#         geom_rug(sides = 'b') +
#         labs(x = Terms[t], y = "smoother", title = paste0(letters[idx-1],") ",'Smoother for ', temp0[[1]]$xlab))
#       
#       ## plot deriv
#       idx <- idx+1
#       CI <- confint(m2.d, term = Terms[t])
#       m2.dtemp <- data.frame(cbind(m2.d$eval[,Terms[t]], m2.d[[Terms[t]]]$deriv, CI[[1]]$upper, CI[[1]]$lower)); names(m2.dtemp) = c('x','deriv','upper','lower')
#       plist[[idx]] <-  ggplot(m2.dtemp,aes(x = x, y = deriv))    +
#         theme_minimal() +
#         theme(panel.grid = element_blank())+
#         geom_line(lwd = 1.1) +
#         geom_hline(yintercept = 0, col = 'grey22') +
#         geom_line(aes(y= upper), linetype = 'dashed') +
#         geom_line(aes(y= lower), linetype = 'dashed') +
#         geom_vline(xintercept = ifelse(is.na(breaksdf[[t]]),-999,breaksdf[[t]]), lwd = 1.1, linetype = 'dashed', col = 'red') +
#         scale_x_continuous(expand = c(0,0), limits = c(ifelse(t == 2, 30, min(m2.dtemp$x)),max(m2.dtemp$x))) +
#         labs(x = Terms[t], y = "f'(x)", title =  paste0(letters[idx-1],") ",'First Derivative for ', temp0[[1]]$xlab))
#       idx <- idx+1
#     }
#     
#     ## map with breaks
#     plist[[1]] <- ggplot() +
#       geom_polygon(data = usa, aes(x = long, y = lat, group = group)) +
#       coord_quickmap() +
#       scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
#       scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
#       theme_minimal() +
#       theme(panel.grid.major = element_blank(),
#             axis.title =element_blank(),
#             legend.position = c(0.2,0.15)) +
#       geom_hline(yintercept =ifelse(is.na(breaksdf[[2]]),-999,breaksdf[[2]]),lwd = 1.1, linetype = 'dashed', col = 'red') +
#       
#       
#       ## had to force this to plot off-map if NA cause throwing 'discrete' error
#       geom_vline(xintercept =  ifelse(is.na(breaksdf[[3]]),-999,breaksdf[[3]]), lwd = 1.1, linetype = 'dashed', col = 'red') +
#       # geom_vline(xintercept =  ifelse(SS == 'F',50,-999), lwd = 1.1, linetype = 'dashed', col = 'red') +
#       ## Owen request to see alternative breaks
#       geom_hline(yintercept =ifelse(SS == 'F' & AA == 6,50,36),lwd = 1, linetype = 'dotted', col = 'grey66') +
#       
#       # geom_point(data = dat, aes(x = Longitude_dd, y = Latitude_dd, size = Length_cm, fill = Length_cm), shape = 21, alpha = 0.7) +
#       scale_fill_viridis_c(guide = "legend") +
#       labs(fill = paste0("Length of Age-",AA," ",SS," Fish (cm)"),
#            size = paste0("Length of Age-",AA," ",SS," Fish (cm)")) +
#       ggtitle(paste0('g) GAM-Estimated Regions')) #+
#     
#     # ggtitle(paste0('g) GAM-Estimated Regions, with raw data Age-',AA,", ",SS)) #+
#     
#     
#     lay <- rbind(c(2,2,3,3,1,1),
#                  c(4,4,5,5,1,1),
#                  c(6,6,7,7,1,1))
#     grid.arrange(grobs = plist, layout_matrix = lay)  %>%
#       ggsave(plot = .,  file = paste0("./figures/sab_smooth_map_",AA,SS,".png"), width = 12, height = 15, units = 'in', dpi = 480)
#   } ## end SS
# } ## end AA



## now re-aggregate and do the fits on all data, based on model ----
# breaksdf <- data.frame(yr_breaks = breaksdf[[1]], lat_breaks2 = breaksdf[[2]], lon_breaks2 = breaksdf[[1]])
for(phase in c("phase1","phase2")[2]){

  breaksdf <- data.frame(yr_breaks = 2010, lat_breaks2 = c(36,50), lon_breaks2 =c(-130,-145)) ## auto
  dat <- getGR(tempdf = full_data, breaksdf) ## assign breakpoints
  dat$selType <- ifelse(dat$REG == 'BC', 2, 1) ## 2 does truncation, 1 leaves as-is
  
  ## add in selex for BC data, using DFO values
  selectivity <- function(Length){
    selex <- 1/(1+exp( 52.976 - Length)) ## email from Sam Johnson/
    return(selex)
  }
  dat$Sel <- ifelse(dat$REG == 'BC', selectivity(dat$Length_cm), 1)
  # dat$Sel <- ifelse(dat$Sel == 0, 1E-5, dat$Sel) ## to avoid div by zero
  
  DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat)) 
  cat(phase," Built Dat \n")
  
  ## sanity check
  # dat %>% group_by(gamREG) %>% summarise(mnlat = mean(Latitude_dd), mnlon = mean(Longitude_dd))
  # dat %>% group_by(gamREG,REG) %>% summarise(n = n())
  
  dat$cREG <- paste0(dat$gamREG,"_",dat$Period,"_",dat$Sex)
  
  ## manual override based on Linf ambiguity ----
  ## phase 2 
  ## males r1/r2
  if(phase == 'phase2'){
    dat$cREG[dat$cREG == "R1_early_M"  | dat$cREG == "R1_late_M" ] <- "R1_pool_M"
    dat$cREG[dat$cREG == "R2_early_M"  | dat$cREG == "R2_late_M" ] <- "R2_pool_M"
    ## both sexes r3/R4//r5
    dat$cREG[dat$cREG == "R3_early_M"  | dat$cREG == "R3_late_M" ] <- "R3_pool_M"
    # dat$cREG[dat$cREG == "R3_early_F"  | dat$cREG == "R3_late_F" ] <- "R3_pool_F"
    
    dat$cREG[dat$cREG == "R4_early_M"  | dat$cREG == "R4_late_M" ] <- "R4_pool_M"
    # dat$cREG[dat$cREG == "R4_early_F"  | dat$cREG == "R4_late_F" ] <- "R4_pool_F"
    
    dat$cREG[dat$cREG == "R5_early_M"  | dat$cREG == "R5_late_M" ] <- "R5_pool_M"
    dat$cREG[dat$cREG == "R5_early_F"  | dat$cREG == "R5_late_F" ] <- "R5_pool_F"
    ## females r4
    length(unique(dat$cREG)) ## 12
  }
  
  ## MK START HERE
  # load("./sabdat_Oct2019_formatted.Rda") 
  ## loads as "dat" saved up to this point since getGR is slow.
  # dat <- sample_n(dat, nrow(dat)*0.25)  %>% filter(selType == 2) ## testing denom
  # DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat))
  DES <- ifelse(!is.na(dat$cREG), as.numeric(as.factor(dat$cREG)),1)-1 ## this is now numeric index, R3=slot 3 (idx 2)
  # KEY <- paste("sab",DES,sep = "_")
  temp <- data.frame(DES = as.numeric(as.character(DES)), cREG = dat$cREG); temp <- unique(temp)
  keybase <- paste(as.factor(temp[order(temp$DES),'cREG']))

  dat0 <- rep0 <- NULL ## later storage
  nStrata <- length(unique(DES))
  

  
  ## this will assign a unique DES depending on period X sex X region -- whatever is in DES
  data <-
    list(
      Length_cm = dat[,"Length_cm"],
      Age = dat[,"Age"],
      DES = as.vector(DES), ## keep this for master iterations
      selType = dat[,'selType'],
      Sel = dat[,'Sel'],
      nStrata = nStrata,
      a2 = 30 ## wtf is this supposed to be
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
  model <- MakeADFun(data, parameters,  DLL="it10",silent=T,map=map)
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
  model$report()$denominator ## if we only ran seltype 2 points, this should NOT be 1.0
  best <- model$env$last.par.best
  rep <- sdreport(model)
  
  dat0 <- rep0 <- NULL ##  storage
  dat0 <- c(dat0, model$report()$ypreds %>% data.frame()) 
  rep0 <- bind_rows(rep0,
                    bind_cols(
                      data.frame(names(rep$value)),
                      data.frame(rep$value),
                      data.frame(rep$sd),
                      data.frame(c(rep(keybase, 5), rep("ALL", 1)))))
  
  ## reformat outputs ----
  names(rep0) <- c('variable', 'value','sd', 'REG')
  write.csv(rep0, file = paste0("./GAM_output/SAB_parEst_gam_",Sys.Date(),'_',phase,'.csv'),row.names = F)
  
  ypreds0 <- cbind(dat0,dat) %>% data.frame()  
  names(ypreds0)[1] <- c('Predicted')
  write.csv(ypreds0,  paste0("./GAM_output/SAB_predicts_",Sys.Date(),'_',phase,".csv"),row.names = F)
  
  cat(phase," Fit TMB model  & saved outputs \n")
}

# rep0 %>% filter(REG == 'R4_early_F')
