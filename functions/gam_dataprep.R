# https://www.fromthebottomoftheheap.net/2014/05/15/identifying-periods-of-change-with-gams/

require(mgcv); require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)
library(ggsidekick)
## load and cbind data
# setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")

## WEST COAST ----
# load("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/raw_data/WC/Bio__NWFSC.Combo_2018-09-25.rda") ## loads as "Data"
# load("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/raw_data/WC/warehouse.RData") ## sent from melissa, no SAB after 2014 -- need to stitch
load(paste0("C:/Users/",compname,"/Dropbox/UW/sab-growth/raw_data/WC/WC_lengths_031519.rda")) ## extracted
data0 <- bio %>% filter(
  Common_name == 'sablefish' &
    !is.na(Age) &
    !is.na(Length_cm) & Sex != 'U' & !is.na(Sex)
) %>% 
  select(
    Year,
    Length_cm ,
    Age ,
    Sex ,
    Latitude_dd ,
    Longitude_dd
  ) %>%
mutate(REG = "WC") #%>%
  # sample_n(.,15000) 
  

# data0 %>%  ggplot(., aes(x = Year)) + geom_histogram() + scale_x_continuous(limits = c(1980,2020),breaks = seq(1983,2020,2))


wcsurv <- data0
wcsurv %>% group_by(Year) %>% summarise(n = n())
# 
# WareHouse.All.Ages.Env %>% filter(common_name == 'sablefish' & year < 2003 & 
#                                     project == "Groundfish Triennial Shelf Survey"
#                                   & !is.na(age_years)) %>% head()
## British Columbia ----
bcsurv <- read.csv("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/raw_data/BC/LWMSO.w_lat_long.csv") %>%
  filter(!is.na(SPECIMEN_AGE) & !is.na(Fork_Length) & 
           SPECIMEN_SEX_CODE %in% c("1","2") &
           NS_AREA != "" & SABLE_AREA_GROUP != "" & slat != 0 &  SABLE_SET_TYPE %in% c('StRS','OFFSHORE STANDARDIZED')) %>%
  select(SPECIMEN_AGE, Fork_Length,SPECIMEN_SEX_CODE,YEAR,slat,slon) %>%
  mutate(Sex = ifelse(SPECIMEN_SEX_CODE == "2", 'F', "M"), 
         Fork_Length = Fork_Length/10) %>%
  plyr::rename(c("slat" = "Latitude_dd",
                 "slon" = "Longitude_dd",
                 "SPECIMEN_AGE" = "Age","Fork_Length" = "Length_cm","YEAR" = "Year")) %>%
  select(Year, Length_cm, Age, Sex, Latitude_dd, Longitude_dd) %>%
  # sample_n(.,15000) %>%
  mutate(REG = "BC")

# ## ALASKA ----
aksurv <- read.csv("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/raw_data/ak/AK_age_view_2018.csv") %>%
  ## drop period before 1995 and filter for top 6 as in Echave
  # filter(., grepl(paste0(c("Southeast",'Kodiak',"Chirikof","Shumagin","Bering","Aleutian"), collapse="|"), GEOGRAPHIC_AREA_NAME)) %>%
  filter(SEX != 3 & !is.na(AGE) & !is.na(LENGTH) ) %>%
  mutate(SEX = ifelse(SEX == 2, 'F', "M")) %>%
  select(YEAR, LENGTH, AGE, SEX, STARTLAT, STARTLONG) %>%
  plyr::rename(c('YEAR' = 'Year', 'SEX' = 'Sex','AGE' = 'Age',
                 'LENGTH' = 'Length_cm', "STARTLAT" = "Latitude_dd","STARTLONG" = "Longitude_dd")) %>%
  # sample_n(.,15000) %>%
  mutate(REG = "AK")

## combine ---
all_data <- rbind(wcsurv,bcsurv,aksurv)
full_data <- rbind(wcsurv,bcsurv,aksurv)
# save(all_data, file = "C:/Users/Maia Kapur/Dropbox/UW/sab-growth/input_data/gam_data_sab_0315.rda")
# save(full_data, file = "C:/Users/Maia Kapur/Dropbox/UW/sab-growth/input_data/gam_data_sab_0415.rda")
# save(all_data, file = "C:/Users/Maia Kapur/Dropbox/UW/coursework/STAT-554/project_code_data/all_data.csv", row.names = F)
## some exploratory plots ----

all_data$REG <- as.factor(all_data$REG)
levels(all_data$REG) <- c("Alaska","Canada","California Current")

ggplot(all_data, aes(x = Year, fill = REG)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.9),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=14))+
  scale_alpha(guide = 'none') +
  labs(x = "Year", fill = 'Region of Origin', y = "# Records") +
  geom_histogram(stat = 'bin', position = 'stack', 
                 binwidth = 1, alpha = 0.6) 


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
