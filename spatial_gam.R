# https://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/

require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)
library(ggsidekick)
## load and cbind data
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")

## WEST COAST ----
load(paste0(getwd(), "/data/raw/WC/Bio__NWFSC.Combo_2018-09-25.rda"), verbose = T) ## loads as "Data"
wcsurv0 <- Data; rm(Data)
wcsurv1 <- wcsurv0 %>% filter(!is.na(Age) & !is.na(Length_cm) & Depth_m < 549 & Depth_m >= 55 & Sex != 'U')
wcsurv <- wcsurv1 %>%
  select(Year, Length_cm, Age, Sex, Latitude_dd, Longitude_dd) %>%   mutate(REG = "WC")
rm(wcsurv0); rm(wcsurv1)

## British Columbia ----
bcsurv <- read.csv(paste0(getwd(),"/data/raw/BC/LWMSO.w_lat_long.csv")) %>%
  filter(!is.na(SPECIMEN_AGE) & !is.na(Fork_Length) & 
           SPECIMEN_SEX_CODE %in% c("1","2") &
           NS_AREA != "" & SABLE_AREA_GROUP != "" & slat != 0) %>%
  select(SPECIMEN_AGE, Fork_Length,SPECIMEN_SEX_CODE,YEAR,slat,slon) %>%
  mutate(Sex = ifelse(SPECIMEN_SEX_CODE == "2", 'F', "M"), 
         Fork_Length = Fork_Length/10) %>%
  plyr::rename(c("slat" = "Latitude_dd",
                 "slon" = "Longitude_dd",
                 "SPECIMEN_AGE" = "Age","Fork_Length" = "Length_cm","YEAR" = "Year")) %>%
  select(Year, Length_cm, Age, Sex, Latitude_dd, Longitude_dd) %>%
  mutate(REG = "BC")


# ## ALASKA ----
aksurv <- read.csv(paste0(getwd(),"/data/raw/ak/AK_age_view_2018.csv")) %>%
  ## drop period before 1995 and filter for top 6 as in Echave
  filter(., grepl(paste0(c("Southeast",'Kodiak',"Chirikof","Shumagin","Bering","Aleutian"), collapse="|"), GEOGRAPHIC_AREA_NAME)) %>%
  filter(SEX != 3 & !is.na(AGE) & !is.na(LENGTH) ) %>%
  mutate(SEX = ifelse(SEX == 2, 'F', "M")) %>%
  select(YEAR, LENGTH, AGE, SEX, STARTLAT, STARTLONG) %>%
  plyr::rename(c('YEAR' = 'Year', 'SEX' = 'Sex','AGE' = 'Age',
                 'LENGTH' = 'Length_cm', "STARTLAT" = "Latitude_dd","STARTLONG" = "Longitude_dd")) %>%
  mutate(REG = "AK")

## combine ---
all_data <- rbind(wcsurv,bcsurv,aksurv)
# save(all_data, file = paste0(getwd(),"/data/gam_data.rda"))
load( paste0(getwd(),"/data/gam_data.rda")) ## all_data

# all_data$Year <- as.factor(all_data$Year)
all_data$Year <- as.numeric(as.character(all_data$Year))
all_data$Sex <- as.factor(all_data$Sex)

## some exploratory plots ----
ggplot(all_data, aes(x = Length_cm, fill = Sex)) +
  theme_sleek() +
  theme(legend.position = c(0.05,0.9))+
  geom_histogram(stat = 'bin', position = 'stack', 
                 binwidth = 1) +
  scale_fill_manual(values = c("#d8b365","#5ab4ac"))+
  facet_wrap(~REG) +
  labs(x = 'Length (cm)')

ggplot(all_data, aes(x = Age, y = Length_cm, color = Sex)) +
  theme_sleek() +
  theme(legend.position = c(0.05,0.9))+
  geom_point( aes( alpha = 0.2) )+
  scale_alpha(guide = 'none')+
  scale_color_manual(values = c("#d8b365","#5ab4ac"))+
  facet_wrap(~REG) +
  labs(y = 'Length (cm)', x = 'Age (years)')

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
mod1 <- gam(Length_cm ~ Age + Sex +
              s(Longitude_dd,Latitude_dd), 
            correlation = corAR1(form = ~ 1|Year, p = 1),
            data = all_data)
mod2 <- gam(Length_cm ~ Age + Sex + s(Year, bs = "cc") +
              s(Longitude_dd,Latitude_dd), 
            correlation = corAR1(form = ~ 1|Year, p = 2),
            data = all_data)
mod3 <- gam(Length_cm ~ Age + Sex + s(Year, bs = "cc") + 
              s(Longitude_dd,Latitude_dd), 
            correlation = corAR1(form = ~ 1|Year, p = 3),
           data = all_data)

MuMIn::model.sel(mod,mod1,mod2,mod3)
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

png( file = paste0(getwd(),"/plots/gam_smooths.png"), height = 6, width = 8, units = 'in', res = 500)
layout(matrix(1:2, ncol = 2))
plot(mod2,  select  =1,  scheme  =2,  lwd  =2)
plot(mod2,  select  =2,  scheme  =2,  lwd  =2)
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
## newDF is data frame of points along X we want to evaluate the derivative, and eps is step
## separate by sex
sapply(subset(model.frame(mod2),Sex == 'M')[,c('Year')],
       function(x) seq(min(x), max(x), length = n))

newDF <- data.frame(Year=seq(1,33,0.5),Longitude_dd= seq(-180,-116,1),Latitude_dd = 0:64)

X0 <- predict(mod2, newDF, type = "lpmatrix")
newDF <- newDF + eps
X1 <- predict(mod, newDF, type = "lpmatrix")
Xp <- (X1 - X0) / eps
## ID breakpoints
## segregate data at breakpoints and fit VonB model using TMB or NLS
