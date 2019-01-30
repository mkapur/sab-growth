## SCRIPT FOR 559 SAB GROWTH PROJECT
## kapurm@uw.edu Fa 19 SAFS

# https://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/

require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2);library(TMB)
# library(ggsidekick)

## data prep and explore ----
# setwd("C:/Users/mkapur/Dropbox/UW/sab-growth") ## fimbria
setwd("C:/Users/Maia Kapur/Dropbox/UW/sab-growth") ## sebastes
load( paste0(getwd(),"/data/gam_data.rda")) ## all_data

all_data$Year <- as.numeric(as.character(all_data$Year))
all_data$Sex <- as.factor(all_data$Sex)
all_data$Longitude_dd <- with(all_data, ifelse(Longitude_dd> 0, -1*Longitude_dd,Longitude_dd))
## change point gam in par est

## first fit with year only and check ACF
mod <- gam(Length_cm ~ s(Year, bs = "cc"), data = all_data)
# mod <- lm(Length_cm ~ Year, data = all_data)
# acf(resid(mod),  main = "ACF")
# 
# mod <- gam(Length_cm ~ Age + Sex + 
#              s(Year, bs = "cc") + s(Latitude_dd), 
#            data = all_data)
# summary(mod)

## plotting model
# layout(matrix(1:2, ncol = 2))
# plot(mod$model, scale = 0)
# layout(1)
# layout(matrix(1:2, ncol = 2))
# pacf(resid(mod), lag.max = 36, main = "pACF")
# layout(matrix(1:2, ncol = 2))

## now try with some AR structures and check AIC
# mod1 <- gam(Length_cm ~ Age + Sex +
#               s(Latitude_dd),
#             correlation = corAR1(form = ~ 1|Year, p = 1),
#             data = all_data)
# mod2 <- gam(Length_cm ~ Age + Sex + s(Year, bs = "cc") +
#               s(Latitude_dd),
#             correlation = corAR1(form = ~ 1|Year, p = 2),
#             data = all_data)
# mod3 <- gam(Length_cm ~ Age + Sex + s(Year, bs = "cc") +
#               s(Latitude_dd),
#             correlation = corAR1(form = ~ 1|Year, p = 3),
#            data = all_data)


MuMIn::model.sel(mod,mod1,mod2,mod3) ## no improvement

# save(mod, file = paste0(getwd(),'/spatial_gam.rda'))
load(paste0(getwd(),'/spatial_gam.rda'))

## check that new model reduced error structure
# layout(matrix(1:2, ncol = 2))
# res <- resid(mod, type = "deviance")
# acf(res, main = "ACF - AR(2) errors")
# pacf(res, main = "pACF- AR(2) errors")

png( file = paste0(getwd(),"/plots/gam_check.png"), height = 6, width = 8, units = 'in', res = 500)
layout(matrix(1:4, ncol = 2))
gam.check(mod)
dev.off()

## get & eval derivatives ----
llsmooth <- mod$smooth[2][[1]]
llsmooth$knots
## calc first derivatives
pdat <- sample_n(all_data,1000)
pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
p2 <- predict(mod, newdata = pdat) ## raw predicts
pdat <- transform(pdat, predLen = p2, se2_spt = pTerm$se.fit[,4], se2_yr = pTerm$se.fit[,3])
df.res <- df.residual(mod)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
## variances are additive just FYI
pdat <- transform(pdat,
                  upper = predLen + (crit.t * (se2_spt+se2_yr)),
                  lower = predLen - (crit.t * (se2_spt+se2_yr)))

## source from gist
source(paste0(getwd(),"/Deriv.R"))

## here I iterate sex as a sanity check but it is actually not necessary -- results will be the same since sex is a covariate
# breaksdf <- data.frame(Sex = NA, YrBreaks = NA, LatBreaks = NA)
breaksdf <- list()
idx <- 1
for(s in 1:2){
  # breaksdf[s,'Sex'] = c("M","F")[s]
  for(t in 1:2){
    
    newD <- data.frame(Age = mean(seq(0,94,length = 100)), ## use mean age?
                       Sex = c("M","F")[s], 
                       Year = seq(1981,2017,length = 100),
                       Longitude_dd = seq(-186,-116, length = 100),
                       Latitude_dd = seq(0,64,length = 100))
    
    dr <- Deriv(mod,newdata = newD)
    confint.Deriv(dr)
    
    Term <- c("Year","Latitude_dd")[t]
    
    m2.d <- Deriv(mod, newdata = newD)
    m2.dci <- confint(m2.d, term = Term)
    
    # crit.eval = mean(m2.d[[Term]]$deriv) ## use mean
    crit.eval = quantile(probs = c(0.025, 0.975), x=  m2.d[[Term]]$deriv) ## use tails
    ## identify where CI does NOT include some crit value (return NA where it does)
    m2.dsig <- signifD(m2.d$eval[[Term]], 
                       d = m2.d[[Term]]$deriv,
                       m2.dci[[Term]]$upper, 
                       m2.dci[[Term]]$lower, eval = crit.eval)
    pix <- c(which(!is.na(m2.dsig$incr)), which(!is.na(m2.dsig$decr)))
    vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
    breaksdf[[idx]] <- sort(c(unique(round(vals)))) ## get rounded unique
    # cat(idx)
    idx <- idx+1
    
  }
}

png( file = paste0(getwd(),"/plots/gam_smooths_update.png"), height = 6, width = 8, units = 'in', res = 500)
layout(matrix(1:4, ncol = 2))
plot(mod,  select  =1,  scheme  
     =2,  lwd  =2, main = 'Year Smoother', cex.axis = 2)
plot(mod,  select  =2,  scheme  =2, xlim = c(35,60),
     lwd  =2, main = 'Latitude Smoother', cex.axis = 2)
plot.Deriv(m2.d, term = 'Year', cex.axis = 2, main = 'derivative of year')
abline(v = breaksdf[[1]], col = 'red')
plot.Deriv(m2.d, term = 'Latitude_dd',
           cex.axis = 2, xlim = c(35,60), main = 'derivative of latitude')
abline(v = breaksdf[[2]], col = 'red')

dev.off()

## now re-aggregate the data ----
## generate DES matrix of vectors and a KEY for later comparison
DES <- KEY <-  matrix(NA, ncol = 1, nrow = nrow(all_data)) ## H1 is all zeros
## get region
# temp0 <- with(all_data, ifelse(Latitude_dd >= 49 & Year < 2005 ,0,
#                                ifelse(Latitude_dd >= 49 & Year >= 2005,1,
#                                       ifelse(Latitude_dd < 49 & Year < 2005,2,3))))
# ## get region and sex combos as factors
# DES <- as.numeric(factor(paste0(temp0,all_data[,'Sex'])))-1
# KEY <- paste0(c("North_Early","North_Late","South_Early","South_Late")[(temp0+1)],"_",all_data[,'Sex'])
keybase <- apply(expand.grid(c("North_Early","North_Late","South_Early","South_Late"), c('F','M')), 1,
      paste, collapse="_")

# save(DES, file = paste0(getwd(),"/data/DES_gam.rda"))
# save(KEY, file = paste0(getwd(),"/data/KEY_gam.rda"))
load( paste0(getwd(),"/data/DES_gam.rda")) ## DES
load( paste0(getwd(),"/data/KEY_gam.rda")) ## KEY

## run TMB parest----
compile("sptlvb.cpp")
dyn.load(dynlib("sptlvb"))

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
    t0 = rep(0, nStrata),
    log_Sigma = 0
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
names(rep0) <- c('variable', 'value','sd', 'ID')
rep0 <- rep0 %>% mutate(
  st =sub('_([^_]*)$', '',  ID),
  Sex = sub(".*_ *(._?)", "\\1", ID)
) %>% select(-ID)

write.csv(rep0, file = paste0(paste0(getwd(),"/results/parEst_gam_",Sys.Date(),'.csv')),row.names = F)

ypreds0 <- cbind(dat0,all_data) %>% data.frame()  
names(ypreds0)[1] <- c('Predicted')
ypreds0$REG <- as.factor(sub('_([^_]*)$', '', KEY))

write.csv(ypreds0,  paste0(getwd(),"/results/predicts",Sys.Date(),".csv"))
ypreds <- read.csv( paste0(getwd(),"/results/predicts",Sys.Date(),".csv"))
## plotting ----

## fits

ggplot(ypreds, aes(x = Age, y = Predicted, col = Sex )) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        legend.position = c(0.9,0.1),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14))+
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  scale_color_manual(values = c("#d8b365","#5ab4ac"))+
  scale_alpha(guide = 'none') +
  geom_point(alpha = 0.2, aes(y = Length_cm)) +
  geom_line(lwd = 1.1) +
  labs(y = 'Length (cm)', col = "") +
  facet_wrap(~REG)
ggsave(file = paste0(getwd(),"/plots/fits_gam",Sys.Date(),".png"), 
       plot = last_plot(), height = 8, width = 10, unit = 'in', dpi = 520)


## parest by sex and region
parest <- read.csv( paste0(getwd(),"/results/parEst_gam_",Sys.Date(),".csv")) %>%
  filter(variable != "Sigma")
## exponentiate logk
parest[parest$variable == 'log_k','value'] <- exp(parest[parest$variable == 'log_k','value'] )
parest$variable <- ifelse(parest$variable=='log_k',"k",paste(parest$variable))

  ggplot(parest, aes(x = st, y = value, col = Sex))+
    theme_bw() +
    theme(panel.grid = element_blank(), 
          legend.position = c(0.25,0.9),
          legend.text = element_text(size = 14),
          strip.text = element_text(size=14),
          axis.text = element_text(angle = 45))+
    scale_color_manual(values = c("#d8b365","#5ab4ac"))+
    geom_point() +
    geom_errorbar(aes(ymin = value - 1.96*sd, ymax = value + 1.96*sd)) +
    labs(x = 'Spatial Stratum', y = "") +
    facet_wrap(~variable, scales = "free_y") 
  
  ggsave(file = paste0(getwd(),"/plots/parplot_gam",Sys.Date(),".png"), 
         plot = last_plot(), height = 6, width = 8, unit = 'in', dpi = 520)

