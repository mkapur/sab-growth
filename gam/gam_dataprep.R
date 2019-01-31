# https://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/

require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)
library(ggsidekick)
## load and cbind data
# setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")

## WEST COAST ----
load("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/data/raw/WC/Bio__NWFSC.Combo_2018-09-25.rda") ## loads as "Data"
load("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/data/warehouse.RData") ## sent from melissa, no SAB after 2014 -- need to stitch
data0 <- WareHouse.All.Ages.Env %>% filter(
  common_name == 'sablefish' &
    !is.na(age_years) &
    !is.na(length_cm) & sex != 'U'
) %>% 
  select(
    Year = year,
    Length_cm = length_cm,
    Age = age_years,
    Sex = sex,
    Latitude_dd = latitude_dd,
    Longitude_dd = longitude_dd,
    REG = Salinity_at_Gear_psu,
    SPECIES_CODE = pass,
    GEAR_DEPTH = depth_ftm,
    Temp = Temperature_at_Surface_c,
    GEAR_TEMPERATURE = Temperature_at_Gear_c
  ) %>% bind_rows(.,
                  Data %>% filter(
                    Year > 2014,
                    !is.na(Age) &
                      !is.na(Length_cm) & Depth_m < 549 & Depth_m >= 55 & Sex != 'U'
                  ))


wcsurv0 <- data0; rm(Data);rm(WareHouse.All.Ages.Env)
wcsurv1 <- wcsurv0 
wcsurv <- wcsurv1 %>%
  select(Year, Length_cm, Age, Sex, Latitude_dd, Longitude_dd) %>%   mutate(REG = "WC")
rm(wcsurv0); rm(wcsurv1)



## British Columbia ----
bcsurv <- read.csv("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/data/raw/BC/LWMSO.w_lat_long.csv") %>%
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
  sample_n(.,8239) %>%
  mutate(REG = "BC")

# ## ALASKA ----
aksurv <- read.csv("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/data/raw/ak/AK_age_view_2018.csv") %>%
  ## drop period before 1995 and filter for top 6 as in Echave
  filter(., grepl(paste0(c("Southeast",'Kodiak',"Chirikof","Shumagin","Bering","Aleutian"), collapse="|"), GEOGRAPHIC_AREA_NAME)) %>%
  filter(SEX != 3 & !is.na(AGE) & !is.na(LENGTH) ) %>%
  mutate(SEX = ifelse(SEX == 2, 'F', "M")) %>%
  select(YEAR, LENGTH, AGE, SEX, STARTLAT, STARTLONG) %>%
  plyr::rename(c('YEAR' = 'Year', 'SEX' = 'Sex','AGE' = 'Age',
                 'LENGTH' = 'Length_cm', "STARTLAT" = "Latitude_dd","STARTLONG" = "Longitude_dd")) %>%
  sample_n(.,8239) %>%
  mutate(REG = "AK")

## combine ---
all_data <- rbind(wcsurv,bcsurv,aksurv)
save(all_data, file = "C:/Users/Maia Kapur/Dropbox/UW/sab-growth/data/gam_data.rda")

## some exploratory plots ----

all_data$REG <- as.factor(all_data$REG)
levels(all_data$REG) <- c("Alaska","Canada","California Current")
p <- ggplot(all_data, aes(x = Length_cm, fill = Sex)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.05,0.9),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14))+
  scale_fill_manual(values = c("#d8b365","#5ab4ac"))+
  scale_alpha(guide = 'none') +
  labs(x = 'Length (cm)', y = "") +
  geom_histogram(stat = 'bin', position = 'stack', 
                 binwidth = 1) +  facet_wrap(~ REG)
ggsave(file =  "C:/Users/mkapur/Dropbox/UW/sab-growth/plots/gam_rawHist.jpg",
       plot = p, height = 8, width = 12, unit = 'in', dpi = 520)

ggplot(all_data, aes(x = Age, y = Length_cm, color = Sex)) +
  # theme_sleek() +
  theme(legend.position = c(0.05,0.9))+
  geom_point( aes( alpha = 0.2) )+
  scale_alpha(guide = 'none')+
  scale_color_manual(values = c("#d8b365","#5ab4ac"))+
  facet_wrap(~REG) +
  labs(y = 'Length (cm)', x = 'Age (years)')
