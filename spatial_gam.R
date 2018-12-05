# https://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/

require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)
library(ggsidekick)
## load and cbind data
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")

load( paste0(getwd(),"/data/gam_data.rda")) ## all_data

# all_data$Year <- as.factor(all_data$Year)
all_data$Year <- as.numeric(as.character(all_data$Year))
all_data$Sex <- as.factor(all_data$Sex)
all_data$Longitude_dd <- with(all_data, ifelse(Longitude_dd> 0, -1*Longitude_dd,Longitude_dd))
## change point gam in par est

## first fit with year only and check ACF
mod <- gam(Length_cm ~ s(Year, bs = "cc"), data = all_data)
# mod <- lm(Length_cm ~ Year, data = all_data)
acf(resid(mod),  main = "ACF")

mod <- gam(Length_cm ~ Age + Sex + s(Year, bs = "cc") + s(Latitude_dd,Longitude_dd), 
           data = all_data)
summary(mod)

## plotting model
# layout(matrix(1:2, ncol = 2))
# plot(mod$gam, scale = 0)
# layout(1)
# layout(matrix(1:2, ncol = 2))
# pacf(resid(mod$lme), lag.max = 36, main = "pACF")layout(matrix(1:2, ncol = 2))

## now try with some AR structures and check AIC
# mod1 <- gam(Length_cm ~ Age + Sex +
#               s(Longitude_dd,Latitude_dd), 
#             correlation = corAR1(form = ~ 1|Year, p = 1),
#             data = all_data)
# mod2 <- gam(Length_cm ~ Age + Sex + s(Year, bs = "cc") +
#               s(Longitude_dd,Latitude_dd), 
#             correlation = corAR1(form = ~ 1|Year, p = 2),
#             data = all_data)
# mod3 <- gam(Length_cm ~ Age + Sex + s(Year, bs = "cc") + 
#               s(Longitude_dd,Latitude_dd), 
#             correlation = corAR1(form = ~ 1|Year, p = 3),
#            data = all_data)


# MuMIn::model.sel(mod,mod1,mod2,mod3)
# save(mod2, file = paste0(getwd(),'/spatial_gam.rda'))
load(paste0(getwd(),'/spatial_gam.rda'))
## support for 2-year lag, not 3

## check that new model reduced error structure
layout(matrix(1:2, ncol = 2))
res <- resid(mod2, type = "deviance")
acf(res, main = "ACF - AR(2) errors")
pacf(res, main = "pACF- AR(2) errors")

png( file = paste0(getwd(),"/plots/gam_check.png"), height = 6, width = 8, units = 'in', res = 500)
layout(matrix(1:4, ncol = 2))
gam.check(mod2)
dev.off()

llsmooth <- mod2$smooth[2][[1]]
llsmooth$knots
## calc first derivatives
pdat <- sample_n(all_data,1000)
pTerm <- predict(mod2, newdata = pdat, type = "terms", se.fit = TRUE)
p2 <- predict(mod2, newdata = pdat) ## raw predicts
pdat <- transform(pdat, predLen = p2, se2_spt = pTerm$se.fit[,4], se2_yr = pTerm$se.fit[,3])
df.res <- df.residual(mod2)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
## variances are additive just FYI
pdat <- transform(pdat,
                  upper = predLen + (crit.t * (se2_spt+se2_yr)),
                  lower = predLen - (crit.t * (se2_spt+se2_yr)))

## source from gist
source(paste0(getwd(),"/Deriv.R"))

##iterate sex
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

dr <- Deriv(mod2,newdata = newD)
confint.Deriv(dr)

Term <- c("Year","Latitude_dd")[t]

m2.d <- Deriv(mod2, newdata = newD)
m2.dci <- confint(m2.d, term = Term)

# crit.eval = mean(m2.d[[Term]]$deriv) ## use mean
crit.eval = quantile(probs = c(0.05, 0.95), x=  m2.d[[Term]]$deriv) ## use tails
## identify where CI does NOT include some crit value (return NA where it does)
m2.dsig <- signifD(m2.d$eval[[Term]], 
                   d = m2.d[[Term]]$deriv,
                   m2.dci[[Term]]$upper, 
                   m2.dci[[Term]]$lower, eval = crit.eval)
pix <- c(which(!is.na(m2.dsig$incr)), which(!is.na(m2.dsig$decr)))
vals <- m2.d$eval[[Term]][pix]
breaksdf[[idx]] <- sort(c(unique(round(vals))))
# cat(idx)
idx <- idx+1

  }
}

png( file = paste0(getwd(),"/plots/gam_smooths.png"), height = 6, width = 8, units = 'in', res = 500)

layout(matrix(1:4, ncol = 2))
plot(mod2,  select  =1,  scheme  =2,  lwd  =2, main = 'Year Smoother', cex.axis = 2)
plot(mod2,  select  =2,  scheme  =2,  lwd  =2, main = '2d Spatial Smoother', cex.axis = 2)
plot.Deriv(m2.d, term = 'Year', cex.axis = 2, main = 'derivative of year')
# abline(v = breaksdf[[1]], col = 'red')
plot.Deriv(m2.d, term = 'Latitude_dd', cex.axis = 2, main = 'derivative of latitude')
# abline(v = breaksdf[[2]], col = 'red')

dev.off()

## I think it co-evaluates the long lat derivs... check if geom_tile matches
data.frame(lat = newD$Latitude_dd, lon = newD$Longitude_dd, val = m2.d$Latitude_dd$deriv) %>%
  ggplot(., aes(x = lon, y = lat, fill = val)) +
  geom_tile()
}
