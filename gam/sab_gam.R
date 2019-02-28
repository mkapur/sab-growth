## Updated GAM method for SAB
## kapurm spring 2019
## source from gist
rm(list = ls())
source(paste0(getwd(),"/Deriv.R"))
compile("sptlvb.cpp") ## will throw 'incomplete line' msg, ignore
dyn.load(dynlib("sptlvb"))

compname <- c("Maia Kapur","mkapur")[1]
require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2);library(TMB)
require(maps);require(mapdata)
## load data
load(paste0("C:/Users/",compname,"/Dropbox/UW/sab-growth//data/gam_data_sab_227.rda")) ## all_data -- made using gam_dataprep
all_data$Longitude_dd[all_data$Longitude_dd > 0] <- all_data$Longitude_dd[all_data$Longitude_dd > 0]*-1
dat <- all_data %>% filter(Age == 4 & Sex == 'F') 
# with(dat, plot(Length_cm ~ Latitude_dd))

mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd), data = dat)
# mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) , data = dat)

## get & eval derivatives ----
llsmooth <- mod$smooth[2][[1]]
# pdat <- sample_n(dat,100)
pdat <- dat
pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
p2 <- predict(mod, newdata = pdat) ## raw predicts
pdat <- transform(pdat, predLen = p2, se2_lat = pTerm$se.fit[,2],se2_lon = pTerm$se.fit[,3], se2_yr = pTerm$se.fit[,1])
# pdat <- transform(pdat, predLen = p2, se2_lat = pTerm$se.fit[,2], se2_yr = pTerm$se.fit[,1])

df.res <- df.residual(mod)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
## variances are additive just FYI
pdat <- transform(pdat,
                  upper = predLen + (crit.t * (se2_lat+se2_lon+se2_yr)),
                  lower = predLen - (crit.t * (se2_lat+se2_lon+se2_yr)))
# pdat <- transform(pdat,
#                   upper = predLen + (crit.t * (se2_lat+se2_yr)),
#                   lower = predLen - (crit.t * (se2_lat+se2_yr)))



newD <- data.frame(Year = seq(1981,2017,length = 100),
                   Longitude_dd = seq(-186,-116, length = 100),
                   Latitude_dd = seq(0,64,length = 100)) ## whatever isn't modeled will just get ignored

breaksdf <- list()
Terms <- c("Year","Latitude_dd","Longitude_dd")
for(t in 1:length(Terms)){
  Term <- Terms[t]
  
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
   labs(x = Terms[t], y = "smoother", title = paste0('Smoother for ', temp0[[1]]$xlab)) 
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
   geom_vline(xintercept = breaksdf[[t]], col = 'red') +
   labs(x = Terms[t], y = "f'(x)", title = paste0('First Derivative for ', temp0[[1]]$xlab)) 
 idx <- idx+1
 # plist[[t+2]] <- plot.Deriv(m2.d, term = Terms[t],cex.axis = 2,  main = paste0('derivative of ', Terms[t])); abline(v = breaksdf[[t]], col = 'red')
}

 ## ! forcing breaksdf based on...intuition
 breaksdf[[2]] <- mean(breaksdf[[2]])
 breaksdf[[3]] <- mean(breaksdf[[3]][3:4])
 usa <- map_data("world") # we already did this, but we can do it again
 plist[[1]] <- ggplot() + 
   geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
   coord_quickmap() +
   scale_x_continuous(expand = c(0,0), limits = c(-180,-110), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
   scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
   theme_minimal() +
   theme(panel.grid.major = element_blank(), axis.title =element_blank()) +
   geom_hline(yintercept = breaksdf[[2]], col = 'red')+
   geom_vline(xintercept = breaksdf[[3]], col = 'red') +
   ggtitle('Midpoint of Lat-Long Breaks') +
   geom_label(aes(x = c(-150, -150,-120), y = c(60,45,45), 
                  label = paste0("GAM Estimated Region ",c(2,3,1))),
              fill = "white",col = 'black')

 
lay <- rbind(c(2,3,1,1),
             c(4,5,1,1),
             c(6,7,1,1))
grid.arrange(grobs = plist, layout_matrix = lay)  %>% 
 ggsave(plot = .,  file = paste0(getwd(),"/plots/sab_smooth_map_a62.png"), width = 15, height = 7, units = 'in', dpi = 480)



## now re-aggregate and do the fits ON ALL DATA ----
DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(dat)) 
## get region
# for(i in 2:length(breaksdf)){ ## loop spatial smooths -- just do simple for now


for (i in 1:nrow(all_data)) {
  all_data[i,'gamREG'] <-
    ifelse(all_data[i, "Latitude_dd"] <=  breaksdf[[2]] & all_data[i, "Longitude_dd"] > breaksdf[[3]], 1,
           ifelse(all_data[i, "Latitude_dd"] > breaksdf[[2]] &  all_data[i, "Longitude_dd"] < breaksdf[[3]], 2, 
                  ifelse(all_data[i, "Latitude_dd"] <= breaksdf[[2]] &  all_data[i, "Longitude_dd"] < breaksdf[[3]], 3, NA)))
                         
}


DES <- ifelse(!is.na(all_data$gamREG), as.numeric(all_data$gamREG),1)-1 ## overwrite NA for nobreaks
KEY <- paste("sab",DES,sep = "_")
keybase <- paste0("R",unique(all_data$gamREG)) ## text regions

dat0 <- rep0 <- NULL ## later storage
nStrata <- length(unique(DES))

data <-
  list(
    Length_cm = all_data[,"Length_cm"],
    Age = all_data[,"Age"],
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
write.csv(rep0, file = paste0(getwd(),"/results/SAB_parEst_1gam_",Sys.Date(),'.csv'),row.names = F)

ypreds0 <- cbind(dat0,all_data) %>% data.frame()  
names(ypreds0)[1] <- c('Predicted')

# write.csv(ypreds0,  paste0(getwd(),"/results/SAB_predicts",Sys.Date(),".csv"),row.names = F)
ypreds <- read.csv( paste0(getwd(),"/results/SAB_predicts",Sys.Date(),".csv"))

cat("Fit TMB model  & saved outputs \n")
## plotting ----



## plot estimates
parest <- read.csv(paste0(getwd(),"/results/SAB_parEst_1gam_",Sys.Date(),'.csv')) %>%
  filter(variable != "Sigma" & variable != 't0') %>% mutate(source = 'Estimated')

## bind only regions of use
parest <- rbind(parest, read.csv("true_sab_vals.csv")) 

## exponentiate logk
parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))
# levels(parest$REG) <- c("ALL","R1","AK","R2","BC","R1","WC")
parest$REG <- factor(parest$REG ,levels=c("ALL","R1","WC","R2","BC","R3","AK"))

plist <- list()
plist[[1]] <- ggplot(parest, aes(x = REG, y = value, col = source))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9,0.9),
        legend.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size=14))+
  scale_color_manual(values = c("red","black"))+
  geom_point() +
  geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
  labs(x = 'Spatial Stratum', y = "", col = "") +
  facet_wrap(~variable, scales = "free_y") 

## plot fits
ypreds <- read.csv( paste0(getwd(),"/results/SAB_predicts",Sys.Date(),".csv"))
ypreds$gamREG <- paste0('GAM-defined Region ',ypreds$gamREG," ", ypreds$Sex)
# ypreds$Sex <- paste0('GAM-defined Region ',ypreds$gamREG)

## fits
plist[[2]] <- ggplot(ypreds[ypreds$Sex == 'F',], aes(x = Age, y = Predicted, col = REG )) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.9),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14))+
  scale_alpha(guide = 'none') +
  geom_point(alpha = 0.2, aes(y = Length_cm)) +
  geom_line(lwd = 1.1, col = 'black')+
  labs(y = 'Length (cm)', col = "Actual Data Source") +
  facet_wrap(~gamREG,ncol = 3)


lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,2,2),
             c(2,2,2,2))
             
grid.arrange(grobs = plist, layout_matrix = lay) %>%  
   ggsave(plot = .,  file = paste0(getwd(),"/plots/sab_fits_a4.png"), width = 10, height = 12, units = 'in', dpi = 480)
  
