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