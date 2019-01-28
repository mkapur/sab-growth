## Generation of data vectors for implementation in TMB growth estimation model.
## Each age or length data point corresponds to the DES (design) matrix, which has a 0 - indexed code
## for the spatial grouping factor under that given hypothesis
## Every data frame should have the following columns:
## H2: for hypothesis 2, 1 of 3 regions
## H3: pooling as suggeH4ed in literature
## H4: pooling as available in survey strata

require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)
setwd("C:/Users/mkapur/Dropbox/UW/sab-growth"

## WEST COAST ----
load(paste0("C:/Users/mkapur/Dropbox/UW/sab-growth/data/raw/WC/Bio__NWFSC.Combo_2018-09-25.rda"), verbose = T) ## loads as "Data"
wcsurv0 <- Data; rm(Data)

## drop wcsurv/age NAs and depths outside of strata
wcsurv1 <- wcsurv0 %>% filter(!is.na(Age) & !is.na(Length_cm) & Depth_m < 549 & Depth_m >= 55 & Sex != 'U')

## merge on strata (based these on the strata_fn)
wcsurv1$H4 <- NULL
for(i in 1:nrow(wcsurv1)){
  if(wcsurv1$Depth_m[i] > 55 & wcsurv1$Depth_m[i] < 200 & wcsurv1$Latitude_dd[i] > 32 & wcsurv1$Latitude_dd[i] < 42){
    wcsurv1[i,'H4'] <- 'shallow_s'
  } else if(wcsurv1$Depth_m[i] >= 200 & wcsurv1$Depth_m[i] < 300 & wcsurv1$Latitude_dd[i] > 32 & wcsurv1$Latitude_dd[i] < 42){
    wcsurv1[i,'H4']<- 'mid_s'
    
  } else if(wcsurv1$Depth_m[i] >= 300 & wcsurv1$Depth_m[i] < 549 & wcsurv1$Latitude_dd[i] > 32 & wcsurv1$Latitude_dd[i] < 42){
    wcsurv1[i,'H4'] <- 'deep_s'
    
  } else if(wcsurv1$Depth_m[i] > 55 & wcsurv1$Depth_m[i] < 200 & wcsurv1$Latitude_dd[i] >= 42 & wcsurv1$Latitude_dd[i] < 49){
    wcsurv1[i,'H4'] <- 'shallow_n'
    
  } else if(wcsurv1$Depth_m[i] >= 200 & wcsurv1$Depth_m[i] < 300 & wcsurv1$Latitude_dd[i] >= 42 & wcsurv1$Latitude_dd[i] < 49){
    wcsurv1[i,'H4'] <- 'mid_n'
    
  } else if(wcsurv1$Depth_m[i] >= 300 & wcsurv1$Depth_m[i] < 549 & wcsurv1$Latitude_dd[i] >= 42 & wcsurv1$Latitude_dd[i] < 49){
    wcsurv1[i,'H4'] <- 'deep_n'
  } 
## split at MTY-CPN per GerH4eva et al. 
  wcsurv1[i,'H3'] <- ifelse(wcsurv1$Latitude_dd[i] >= 36, "North MTY-CPN","South MTY-CPN")
}

wcsurv1$H4 <- factor(wcsurv1$H4, levels=c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s"))
wcsurv <- wcsurv1 %>% mutate(H2 = 'WC', GEAR = NA) %>% select(Age, Length_cm, Sex, H2, H3, H4,GEAR)
rm(wcsurv1); rm(wcsurv0)

## British Columbia ----
bcsurv <- read.csv("C:/Users/mkapur/Dropbox/UW/sab-growth/data/raw/BC/BC_LWMSO_1970-present.csv") %>% 
  filter(!is.na(SPECIMEN_AGE) & !is.na(Fork_Length) & SPECIMEN_SEX_CODE %in% c("1","2") & NS_AREA != "" & SABLE_AREA_GROUP != "") %>%
  select(SPECIMEN_AGE, Fork_Length,SPECIMEN_SEX_CODE,SABLE_AREA_GROUP,NS_AREA,GEAR_CODE) %>%
  mutate(H2 = 'BC', H3 = NS_AREA,  Sex = ifelse(SPECIMEN_SEX_CODE == "2", 'F', "M"), Fork_Length = Fork_Length/10)  %>%
  plyr::rename(c("SPECIMEN_AGE" = "Age","Fork_Length" = "Length_cm", "SABLE_AREA_GROUP" = "H4", "GEAR_CODE" = 'GEAR')) %>%
  select(Age, Length_cm, Sex, H2, H3, H4, GEAR)

## beta -- subsample
bcsurv <- sample_n(bcsurv, mean( nrow(aksurv),nrow(wcsurv)))
  

# ALASKA ----
aksurv <- read.csv("C:/Users/mkapur/Dropbox/UW/sab-growth/data/raw/ak/AK_age_view_2018.csv") %>%
  ## drop period before 1995 and filter for top 6 as in Echave
  filter(., grepl(paste0(c("Southeast",'Kodiak',"Chirikof","Shumagin","Bering","Aleutian"), collapse="|"), GEOGRAPHIC_AREA_NAME)) %>%
  filter(SEX != 3 & !is.na(AGE) & !is.na(LENGTH) & YEAR > 1995) %>% 
  mutate(H2 = 'AK',
         H3 = ifelse(COUNCIL_SABLEFISH_MGMT_AREA == 'Bering Sea', 'Bering',
                     
                     ifelse(COUNCIL_SABLEFISH_MGMT_AREA == "Aleutians", "Aleutians",
                            ifelse(COUNCIL_SABLEFISH_MGMT_AREA == "East Yakutat/Southeast", "Southeast",
                                   as.character(GEOGRAPHIC_AREA_NAME )))), 
         SEX = ifelse(SEX == 2, 'F', "M"),
         GEAR = NA) %>%
  select(AGE, LENGTH, SEX, H2, H3, GEOGRAPHIC_AREA_NAME, GEAR) %>%
  plyr::rename(c('SEX' = 'Sex','AGE' = 'Age', 'LENGTH' = 'Length_cm', "GEOGRAPHIC_AREA_NAME" = "H4")) 

## combine ---
all_data <- rbind(wcsurv,bcsurv,aksurv) ## includes gear

# all_data <- all_data %>% select(-GEAR); save(all_data, file = paste0(getwd(),"/data/all_data.rda"))



## generate DES matrix of vectors and a KEY for later comparison
DES <- KEY <-  matrix(0, ncol = 4, nrow = nrow(all_data)) ## H1 is all zeros
for(d in 1:4){
  if(d == 1) temp <- paste0(all_data$Sex,"_ALL")
  else if(d != 1) temp <- paste0(all_data$Sex,"_",all_data[,d+2])
  DES[,d] <- as.numeric(factor(temp))-1 ## consider sex too, zero index
  KEY[,d] <- temp
}
REG<-KEY[,2] %>% data.frame("H2" = .) %>% mutate('REGION' = as.numeric(factor(sub(".*_ *(._?)", "\\1",.$H2)))) %>% select(REGION) 


# save(DES, file = paste0(getwd(),"/data/DES.rda"))
# save(KEY, file = paste0(getwd(),"/data/KEY.rda"))
# save(REG, file = paste0(getwd(),"/data/REG.rda"))

## generate selectivity vector for any BC data


## matrix of params by gear
selmat <- data.frame('gear' =1:8, 'L50min' = sample(20:30,8), 'L50max' = sample(35:42,8), "L95min" = sample(45:52,8), "L95max" = sample(50:60,8))
selfun <- function(Li, gear){
  t1 <- 1 + exp(-log(19)*(Li - selmat[gear,'L50min'])/(selmat[gear,'L95min']-selmat[gear,'L50min']))
  t2 <- 1 + exp(-log(19)*(Li - selmat[gear,'L50max'])/(selmat[gear,'L95max']-selmat[gear,'L50max']))
  return(t1/t2)
}

SEL <- all_data %>% mutate(Sel = ifelse(H2 == 'BC', selfun(Li = Length_cm,gear =  as.numeric(GEAR)), NA)) %>% select(Sel)
# save(SEL, file = paste0(getwd(),"/data/Sel.rda"))

plot( SEL[!is.na(SEL)] ~ all_data$Length_cm[all_data$H2 == 'BC'], type = 'p') ## looks logisitic...


