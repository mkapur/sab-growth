## Updated GAM method for SAB
## kapurm spring 2019
## source from gist
source(paste0(getwd(),"/Deriv.R"))
compname <- c("Maia Kapur","mkapur")[1]
require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2);library(TMB)

## load data
load(paste0("C:/Users/",compname,"/Dropbox/UW/sab-growth//data/gam_data_sab_227.rda")) ## all_data -- made using gam_dataprep

dat <- all_data %>% filter(Age == 4 & Sex == 'F')
# with(dat, plot(Length_cm ~ Latitude_dd))

mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd), data = dat)

## get & eval derivatives ----
llsmooth <- mod$smooth[2][[1]]
# pdat <- sample_n(dat,100)
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
# pdat <- transform(pdat,
#                   upper = predLen + (crit.t * (se2_spt)),
#                   lower = predLen - (crit.t * (se2_spt)))



newD <- data.frame(Year = seq(1981,2017,length = 100),
                   Longitude_dd = seq(-186,-116, length = 100),
                   Latitude_dd = seq(0,64,length = 100))

# dr <- Deriv(mod,newdata = newD)
# confint.Deriv(dr)
breaksdf <- list()
Terms <- c("Year","Latitude_dd","Longitude_dd")
for(t in 1:3){
Term <- Terms[t]

m2.d <- Deriv(mod, newdata = newD)
m2.dci <- confint(m2.d, term = Term)

# crit.eval = mean(m2.d[[Term]]$deriv) ## use mean
crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
crit.eval.se = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$se.deriv) ## use tails

## identify where falls outside CI (NAs where inside CI)
# m2.dsig.outliers <- signifD(m2.d$eval[[Term]],
#                    d = m2.d[[Term]]$deriv,
#                    m2.dci[[Term]]$upper,
#                    m2.dci[[Term]]$lower,
#                    eval = crit.eval)

## identify where mean crosses zero or falls out of bounds (NAs where it does)
m2.dsig.zeros <- signifMK(x = m2.d$eval[[Term]],
                          d = m2.d[[Term]]$deriv,
                          upper = m2.dci[[Term]]$upper,
                          lower = m2.dci[[Term]]$lower,
                          crit.eval = crit.eval,
                          eval = 0)

## of passing ones, find outliers

## identify where CI crosses zero or falls out of bounds (NAs where it does)
# m2.dsig.CIs <- signifCI(x = m2.d$eval[[Term]],
#                         upper = m2.dci[[Term]]$upper,
#                         lower = m2.dci[[Term]]$lower,
#                          eval = 0)


## only retain vals that pass both crit
# outkeep <- !is.na(m2.dsig.outliers$incr)[!is.na(m2.dsig.outliers$decr)]
# zkeep <- !is.na(m2.dsig.outliers$incr)[!is.na(m2.dsig.outliers$decr)]
# ikeep <- !is.na(m2.dsig.zeros$incr)[!is.na(m2.dsig.CIs$incr)]
# dkeep <- !is.na(m2.dsig.zeros$decr)[!is.na(m2.dsig.CIs$decr)]
# 

pix <- !is.na(m2.dsig.zeros)
# pix <- outkeep[zkeep]
# pix <- zkeep[cikeep]
vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
breaksdf[[t]] <- sort(c(unique(round(vals)))) ## get rounded unique
# cat(idx)
# idx <- idx+1

}

# png(file = paste0(outdir,"/",b,"_gam_smooths.png"), height = 6, width = 8, units = 'in', res = 500)
layout(matrix(1:6, ncol = 2, byrow = T))
for(t in 1:3){
  plot(mod,   select = t, scheme  =2, lwd  =2, main = paste0(Terms[t],' Smoother'), cex.axis = 2, ylim = c())
  plot.Deriv(m2.d, term = Terms[t],cex.axis = 2,  main = paste0('derivative of ', Terms[t]))
  abline(v = breaksdf[[t]], col = 'red')
}

graphics.off()

# png(file = paste0(outdir,"/",b,"_gamAll.png"), height = 6, width = 8, units = 'in', res = 500)
# layout(matrix(1:6, ncol = 2,byrow = T))
# gam.check(mod)

plot(mod,   select = 2, scheme  =2, lwd  =2, main = 'Latitude Smoother', cex.axis = 2)
plot.Deriv(m2.d, term = 'Latitude_dd',cex.axis = 2,  main = 'derivative of latitude')
abline(v = breaksdf[[1]], col = 'red')
graphics.off()
return(breaksdf)
  