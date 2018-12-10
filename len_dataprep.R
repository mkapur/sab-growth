require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)

setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")


# ALASKA ----
aksurv <- read.csv(paste0(getwd(),"/data/raw/ak/AK_age_view_2018.csv")) %>%
  ## drop period before 1995 and filter for top 6 as in Echave
  filter(., grepl(paste(c("Southeast",'Kodiak',"Chirikof","Shumagin","Bering","Aleutian"), collapse="|"), GEOGRAPHIC_AREA_NAME)) %>%
  filter(SEX != 3 & !is.na(AGE) & !is.na(LENGTH) & YEAR > 1995) %>% 
  mutate(st = ifelse(COUNCIL_SABLEFISH_MGMT_AREA == 'Bering Sea', 'Bering',
                     
                     ifelse(COUNCIL_SABLEFISH_MGMT_AREA == "Aleutians", "Aleutians",
                            ifelse(COUNCIL_SABLEFISH_MGMT_AREA == "East Yakutat/Southeast", "Southeast",
                                                                                   as.character(GEOGRAPHIC_AREA_NAME )))), 
         SEX = ifelse(SEX == 2, 'F', "M")) %>%
  select(AGE, LENGTH, SEX, st, GEOGRAPHIC_AREA_NAME) %>%
  plyr::rename(c('SEX' = 'Sex','AGE' = 'Age', 'LENGTH' = 'Length_cm', "GEOGRAPHIC_AREA_NAME" = "st2"))



# from fishery 
aksurv$st <- factor(aksurv$st)
aksurv$st2 <- factor(aksurv$st2)
# aksurv$st <- factor(aksurv$st, levels = c("AK 0","AK 1","AK 2","AK 3","AK 4","AK 5","AK 6","AK 7","AK 8","AK 9"))
save(aksurv[,c(1:4)], file = paste0(getwd(),"/filtered_SAB_AK_1.rda")) ## raw
save(aksurv[,c(1:3,5)], file = paste0(getwd(),"/filtered_SAB_AK_2.rda")) ## raw


load( paste0(getwd(),"/filtered_SAB_AK_1.rda")) ## aksurv

nStrata <- length(unique(aksurv$st))
aka0 <-
  aksurv %>% select(Age, st, Sex) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID',"Sex"),'Age') %>% 
  reshape::cast(., ID + value ~ st + Sex, fill = NA, drop= T) %>% select(-ID, -value) 
akl0 <- 
  aksurv %>% select(Length_cm, Sex, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID',"Sex"),'Length_cm') %>% 
  reshape::cast(., ID + value ~ st + Sex, fill = NA, drop= T) %>% select(-ID, -value) 

aka <- akl <-  array(NA, dim = c(34069, length(unique(aksurv$st)), 2))

# matrix(NA, nrow = 8481, ncol = length(unique(len$st)))
aknmat <- matrix(NA, nrow = nStrata, ncol = 2)
for(s in 1:2){
  Atemp <- aka0[,grep(c("M","F")[s],names(aka0))]
  Ltemp <- akl0[,grep(c("M","F")[s],names(akl0))]
  for(i in 1:nStrata){
    # idx <- (s-1)*6+i
    nr <- length(na.omit(Atemp[,i]))
    # aka[1:nr,i] <- data.frame(na.omit(aka0[,i])) ##if there is a value in one, should be in the other as well
    # aka[,i] <- bind_cols(data.frame(aka), data.frame(na.omit(aka0[,i])))
    aka[1:nr,i,s] <- matrix(na.omit(Atemp[,i]),ncol = 1)
    akl[1:nr,i,s] <- matrix(na.omit(Ltemp[,i]),ncol = 1)
    aknmat[i,s] <- nr ## store lengths of each
    # idx <- sample(1:nrow(temp),500) ## pick 500 rows and store index
    # aka[,i] <- temp[idx,] %>% as.matrix()
    # akl[,i] <- data.frame(na.omit(akl0[,i]))[idx,] %>% as.matrix()
  }
}


save(aka, file = paste0(getwd(),"/agearray_AK.rda")); rm(aka0); rm(Atemp)
save(akl, file = paste0(getwd(),"/lenarray_AK.rda")); rm(akl0); rm(Ltemp)
save(aknmat, file = paste0(getwd(),"/nmat_AK.rda"))

## WEST COAST ----
load("C:/Users/mkapur/Dropbox/UW/sab-growth/data/raw/WC/Bio__NWFSC.Combo_2018-09-25.rda") ## loads as "Data"
len0 <- Data; rm(Data)

## drop len/age NAs and depths outside of strata
len1 <- len0 %>% filter(!is.na(Age) & !is.na(Length_cm) & Depth_m < 549 & Depth_m >= 55 & Sex != 'U')

## merge on strata (based these on the strata_fn)
len1$st <- NULL
for(i in 1:nrow(len1)){
       if(len1$Depth_m[i] > 55 & len1$Depth_m[i] < 200 & len1$Latitude_dd[i] > 32 & len1$Latitude_dd[i] < 42){
         len1[i,'st'] <- 'shallow_s'
       } else if(len1$Depth_m[i] >= 200 & len1$Depth_m[i] < 300 & len1$Latitude_dd[i] > 32 & len1$Latitude_dd[i] < 42){
         len1[i,'st']<- 'mid_s'

       } else if(len1$Depth_m[i] >= 300 & len1$Depth_m[i] < 549 & len1$Latitude_dd[i] > 32 & len1$Latitude_dd[i] < 42){
         len1[i,'st'] <- 'deep_s'

       } else if(len1$Depth_m[i] > 55 & len1$Depth_m[i] < 200 & len1$Latitude_dd[i] >= 42 & len1$Latitude_dd[i] < 49){
         len1[i,'st'] <- 'shallow_n'
         
       } else if(len1$Depth_m[i] >= 200 & len1$Depth_m[i] < 300 & len1$Latitude_dd[i] >= 42 & len1$Latitude_dd[i] < 49){
         len1[i,'st'] <- 'mid_n'
         
       } else if(len1$Depth_m[i] >= 300 & len1$Depth_m[i] < 549 & len1$Latitude_dd[i] >= 42 & len1$Latitude_dd[i] < 49){
         len1[i,'st'] <- 'deep_n'
       }
}
len <- len1; rm(len1)
len$st_f <- factor(len$st, levels=c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s"))
# save(len, file = paste0(getwd(),"/filtered_SAB_WC.rda")) ## raw


## some reshaping of the data, see https://stackoverflow.com/questions/28036294/collapsing-rows-where-some-are-all-na-others-are-disjoint-with-some-nas
## don't forget this is a ragged array, looks like TMB can't deal with NAs -- currently subsampling 500 data points from each strata :/
load( paste0(getwd(),"/filtered_SAB_WC.rda")) ## len

nStrata <- length(unique(len$st_f))
wca0 <-
  len %>% filter(Sex != 'U') %>% select(Age, st, Sex) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID',"Sex"),'Age') %>% 
  reshape::cast(., ID + value ~ st + Sex, fill = NA, drop= T) %>% select(-ID, -value) 
wcl0 <- 
  len %>% filter(Sex != 'U') %>% select(Length_cm, Sex, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID',"Sex"),'Length_cm') %>% 
  reshape::cast(., ID + value ~ st + Sex, fill = NA, drop= T) %>% select(-ID, -value) 

## omit NAs and select 500 from each strata
wca <- wcl <-  array(NA, dim = c(8481, length(unique(len$st)), 2))
  
# matrix(NA, nrow = 8481, ncol = length(unique(len$st)))
wcnmat <- matrix(NA, nrow = nStrata, ncol = 2)
for(s in 1:2){
  Atemp <- wca0[,grep(c("M","F")[s],names(wca0))]
  Ltemp <- wcl0[,grep(c("M","F")[s],names(wcl0))]
  for(i in 1:nStrata){
    # idx <- (s-1)*6+i
    nr <- length(na.omit(Atemp[,i]))
    # wca[1:nr,i] <- data.frame(na.omit(wca0[,i])) ##if there is a value in one, should be in the other as well
    # wca[,i] <- bind_cols(data.frame(wca), data.frame(na.omit(wca0[,i])))
    wca[1:nr,i,s] <- matrix(na.omit(Atemp[,i]),ncol = 1)
    wcl[1:nr,i,s] <- matrix(na.omit(Ltemp[,i]),ncol = 1)
    wcnmat[i,s] <- nr ## store lengths of each
    # idx <- sample(1:nrow(temp),500) ## pick 500 rows and store index
    # wca[,i] <- temp[idx,] %>% as.matrix()
    # wcl[,i] <- data.frame(na.omit(wcl0[,i]))[idx,] %>% as.matrix()
  }
}





save(wca, file = paste0(getwd(),"/agearray_WC.rda")); rm(wca0); rm(Atemp)
save(wcl, file = paste0(getwd(),"/lenarray_WC.rda")); rm(wcl0); rm(Ltemp)
save(wcnmat, file = paste0(getwd(),"/nmat_WC.rda"))



## combine ----
load( paste0(getwd(),"/nmat_WC.rda"))
load(paste0(getwd(),"/lenarray_WC.rda"))
load(paste0(getwd(),"/agearray_WC.rda")) 

load( paste0(getwd(),"/nmat_AK.rda"))
load(paste0(getwd(),"/lenarray_AK.rda"))
load(paste0(getwd(),"/agearray_AK.rda")) 

lenmat <- agemat <-  array(NA, dim = c(34069, ncol(wcl) +ncol(akl), 2))
for(i in 1:2){
  lenmat[,,i] <-  as.matrix(rowr::cbind.fill(wcl[,,i],akl[,,i]))
  agemat[,,i] <-  as.matrix(rowr::cbind.fill(wca[,,i],aka[,,i]))
}
nmat <- rbind(wcnmat,aknmat)


## create bin assignments (may vary by region)
Bin <- lenmat


## create selectivity matrix
Sel <- array(NA, dim = c(34069, 15, 2))

for(s in 1:2){
  for(j in 1:nStrata){
    L50 <- runif(1,30,55) ## random inflection point @ strata
    for(i in 1:nmat[j,s]){ ## only do non NAs
      Sel[i,j,s] <-  1/(1+exp(L50 - lenmat[i,j,s])) ## selectivity @ obs given l50; will be same for each strata
    }
  }
}
  
  


plot(Sel[,1,1] ~ lenmat[,1,1], pch = 19, main = 'simulated selectivity curves, n = 15')
for(i in 1:ncol(Sel)) points(Sel[,i,1]~ lenmat[,i,1], col =rainbow(30)[i], pch = 19)

save(agemat, file = paste0(getwd(),"/agearray.rda")); rm(wca0); rm(Atemp)
save(lenmat, file = paste0(getwd(),"/lenarray.rda")); rm(wcl0); rm(Ltemp)
save(nmat, file = paste0(getwd(),"/nmat.rda"))
save(Sel, file = paste0(getwd(),"/selarray.rda"))


# len %>% group_by(st) %>% summarise(minL = min(Length_cm), maxL = max(Length_cm))
# len %>% group_by(st) %>% summarise(quantile(Length_cm, probs = 0.025),quantile(Length_cm, probs = 0.5),quantile(Length_cm, probs = 0.975))

## Plotting ---
## raw len@age, by strata and sex
## change order for order
ggplot(len, aes(x = Age, y = Length_cm, color = Sex)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = c(0.9,0.15))+
  scale_y_continuous(limits = c(0,100)) +
  scale_x_continuous(limits = c(0,50)) +
  scale_color_brewer(palette = 'Dark2') +
  geom_point(alpha = 0.2) +
  facet_wrap(~ st_f) +
  labs(title = "Raw Data", y = 'Length (cm)', x= 'Age (yr)', subtitle = 'NWFSC Groundfish Survey')
ggsave(file = paste0(getwd(),"/plots/raw_data.png"), plot = last_plot(), height = 5, width = 7, unit = 'in', dpi = 520)

len.bins = 22:90 ## from Kelli's assesssment

# Calculate the effN
n = GetN.fn(dir=getwd(), dat = len, type = "length", species = "shelfrock", printfolder = "forSS")

# The GetN.fn calculated input sample sizes based on Hamel & Stewart bootstrap approach.

# Expand and format length composition data for SS
LFs <- SurveyLFs.fn(dir = getwd(), datL = len, datTows = catch,  
                    strat.df = strata, lgthBins = len.bins, gender = 3, 
                    sexRatioStage = 2, sexRatioUnsexed = 0.5, maxSizeUnsexed = 26, 
                    nSamps = n)

# The code offers two options for applying the sex ratio based on expansion stage. The sex ratio will be
# applied based on a tow basis first if sexRatioStage = 1. The other option applies the sex ratio to the
# expanded numbers of fish across a whole strata (sexRatioStage = 2, this was the option applied to the
# NWFSC combo survey data in the past).


PlotFreqData.fn(dir = getwd(), LFs, survey = "NWFSCBT", ylim=c(0, max(len.bins) + 4), yaxs="i", ylab="Length (cm)", dopng = TRUE)
PlotSexRatio.fn(dir = getwd(), dat = len, data.type = "length", survey = "NWFSCBT", dopng = TRUE, main = "NWFSCBT")

#============================================================================================
#Length Biological Data 
#============================================================================================
age = len
age.bins = 1:40
n = GetN.fn(dir = getwd(), dat = age, type = "age", species = "shelfrock", printfolder = "forSS")
# Exand and format the marginal age composition data for SS
Ages <- SurveyAFs.fn(dir = getwd(), datA = age, datTows = catch,  
                     strat.df = strata, ageBins = age.bins, 
                     sexRatioStage = 2, sexRatioUnsexed = 0.50, maxSizeUnsexed = 5, 
                     gender = 3, nSamps = n)



PlotFreqData.fn(dir = getwd(), dat = Ages, survey = "NWFSCBT", ylim=c(0, max(age.bins) + 2), yaxs="i", ylab="Age (yr)", dopng=TRUE)
PlotVarLengthAtAge.fn(dir = getwd(), dat = age, survey ="NWFSCBT", dopng = TRUE) 
PlotSexRatio.fn(dir = getwd(), dat = age, data.type = "age", survey = "NWFSCBT", dopng = TRUE, main = "NWFSCBT")



## dataprep for mapping
load( paste0(getwd(),"/data/filtered_SAB_WC.rda")) ## len

wcsurv <- len %>% select(Latitude_dd,Longitude_dd,st_f) %>%
  plyr::rename(c("Latitude_dd" = "LAT", "Longitude_dd" = "LON", "st_f" = "STRATA")) %>%
  mutate(REGION = 'WC')

aksurv <- read.csv(paste0(getwd(),"/data/raw/AK/AK_age_view_2018.csv")) %>% 
  select(STARTLAT,STARTLONG, STRATUM) %>% 
  plyr::rename(c("STARTLAT" = "LAT", "STARTLONG" = "LON","STRATUM" = "STRATA")) %>%
  mutate(REGION = 'AK', STRATA = as.character('STRATA'))

## no lat and longs in BC data just area code -- maybe do runif within limits
bcsurv <- read.csv(paste0(getwd(),"/data/raw/BC/BC_LWMSO_1970-present.csv")) %>% 
  select(SABLE_AREA_GROUP) %>%
  plyr::rename(c("SABLE_AREA_GROUP" = "STRATA")) %>%
  mutate(LAT = runif(nrow(.),48,58), LON = runif(nrow(.),-135,-128), REGION = 'BC')  %>%
  select(LAT,LON,STRATA,REGION)

mapdf <- bind_rows(wcsurv,bcsurv,aksurv)
# write.csv(mapdf, file = paste0(getwd(),'/data/mapdf.csv'),row.names = F)
