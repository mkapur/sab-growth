require(dplyr); require(reshape); require(RColorBrewer); require(ggplot2)

setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")

load(paste0(getwd(), "/Bio__NWFSC.Combo_2018-09-25.rda"), verbose = T) ## loads as "Data"
len0 <- Data; rm(Data)

## drop len/age NAs and depths outside of strata
len1 <- len0 %>% filter(!is.na(Age) & !is.na(Length_cm) & Depth_m < 549 & Depth_m >= 55)

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
len$st_f = factor(len$st, levels=c('deep_n','mid_n','shallow_n','deep_s', "mid_s","shallow_s"))

# save(len, file = paste0(getwd(),"/filtered_SAB_WC.rda")) ## raw
## some reshaping of the data, see https://stackoverflow.com/questions/28036294/collapsing-rows-where-some-are-all-na-others-are-disjoint-with-some-nas
## don't forget this is a ragged array, looks like TMB can't deal with NAs -- currently subsampling 500 data points from each strata :/
load( paste0(getwd(),"/filtered_SAB_WC.rda"))


agemat0 <-
  len %>% select(Age, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID'),'Age') %>% 
  reshape::cast(., ID + value ~ st, fill = NA, drop= T) %>% select(-ID, -value) 
lenmat0 <- 
  len %>% select(Length_cm, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID'),'Length_cm') %>% 
  reshape::cast(., ID + value ~ st, fill = NA, drop= T) %>% select(-ID, -value) 

## omit NAs and select 500 from each strata
agemat <- lenmat <-  matrix(NA, nrow = 500, ncol = length(unique(len$st)))
for(i in 1:nStrata){
  temp <- data.frame(na.omit(agemat0[,i])) ##if there is a value in one, should be in the other as well 
  idx <- sample(1:nrow(temp),500) ## pick 500 rows and store index
  agemat[,i] <- temp[idx,] %>% as.matrix()
  lenmat[,i] <- data.frame(na.omit(lenmat0[,i]))[idx,] %>% as.matrix()
}

save(agemat, file = paste0(getwd(),"/agemat_WC.rda")); rm(agemat0)
save(lenmat, file = paste0(getwd(),"/lenmat_WC.rda")); rm(lenmat0)

len %>% group_by(st) %>% summarise(minL = min(Length_cm), maxL = max(Length_cm))
len %>% group_by(st) %>% summarise(quantile(Length_cm, probs = 0.025),quantile(Length_cm, probs = 0.5),quantile(Length_cm, probs = 0.975))

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

len.bins = 11:47

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



