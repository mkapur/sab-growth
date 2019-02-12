## not needed for GAM analysis - moved here Feb 11 2019

## Generate Catch Matrix
CATCH_MATRIX = subset(IBM_DEAD_MATRIX, dtype=="FISHED") ## all records of fishing death
write.table(CATCH_MATRIX,paste0(path_name, "/CATCH_MATRIX.csv",sep=""),row.names = FALSE,col.names = TRUE,sep=",")

IBM_EXP_MATRIX = read.table(paste0(path_name,"/IBM_EXP_MATRIX.txt"),sep=",",header=T) 
CATCH_MATRIX = read.table(paste0(path_name,"/CATCH_MATRIX.csv"),sep=",",header=T) 

SS3_Obs_Size_COM_MinBin = 5
SS3_Obs_Size_COM_MaxBin = 1005
Obs_Size_COM_bin_interval = 10
Effect_N = 100

OBS_size_bin = seq(SS3_Obs_Size_COM_MinBin,SS3_Obs_Size_COM_MaxBin,Obs_Size_COM_bin_interval) 
lower_OBS_size_bin = OBS_size_bin[-length(OBS_size_bin)]
N_obs_size_bin = length(lower_OBS_size_bin)

OBS_SIZE_COM = NULL
por_TRUE_SIZE_COM = NULL
por_OBS_SIZE_COM = NULL

for(y in (start_yr+M_only_yr):(start_yr+simu_year-1)){
  temp_Year = subset(CATCH_MATRIX, Year==y)
  
  if(length(temp_Year[,1])>=1){
    # lower bin
    hist_data = hist(temp_Year$fish_size,right=F,breaks=OBS_size_bin,plot=FALSE)  
    temp_true_catch_com = hist_data$counts/sum(hist_data$counts)
    por_TRUE_SIZE_COM = rbind(por_TRUE_SIZE_COM,temp_true_catch_com)
    
    catch_com_new = NULL
    temp_por_OBS_SIZE_COM = NULL
    for(i in 1:length(temp_true_catch_com)){
      n_count = 0
      for(j in 1:Effect_N){
        r4=runif(1, min=0, max=1)
        if(r4<=temp_true_catch_com[i]){
          n_count=n_count+1
        }
      }
      catch_com_new[i]=n_count
      temp_por_OBS_SIZE_COM[i] = n_count/Effect_N
    }
    por_OBS_SIZE_COM = rbind(por_OBS_SIZE_COM,temp_por_OBS_SIZE_COM)
    #Yr","Seas","Flt","Gender","Part","Nsamp"
    OBS_SIZE_COM = rbind(OBS_SIZE_COM,c(y,   1,      1,    0,      0,      Effect_N,catch_com_new))
  } # end of if 
}
colnames(OBS_SIZE_COM) <-  c("#Yr","Seas","Flt","Gender","Part","Nsamp",lower_OBS_size_bin)
colnames(por_OBS_SIZE_COM) <-  c(lower_OBS_size_bin)


write.table("#_N_LengthBins",paste0(path_name, "/SS3_OBS_Size_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(N_obs_size_bin,paste0(path_name, "/SS3_OBS_Size_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)
write.table(t(lower_OBS_size_bin),paste0(path_name, "/SS3_OBS_Size_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)

write.table("#_N_Length_obs",paste0(path_name, "/SS3_OBS_Size_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  
N_size_observation = length(OBS_SIZE_COM[,1])
write.table(t(N_size_observation),paste0(path_name, "/SS3_OBS_Size_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 

# data-vector
write.table(t(colnames(OBS_SIZE_COM)),paste0(path_name, "/SS3_OBS_Size_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(OBS_SIZE_COM,paste0(path_name, "/SS3_OBS_Size_COM.csv"),row.names=FALSE,col.names=FALSE,sep=",",append = T)

# print(paste0('wrote SS3 OBS size comp, sim ',sim_number," FLEVEL ",level))

## SS3 AGE COMPOSITION ----
# file.remove(paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""))

max_IBM_age =  max(max(CATCH_MATRIX$AgeE),max(CATCH_MATRIX$Age),max(IBM_SAA_MATRIX$Age))
if(max_IBM_age>a2){
  full_age_breaks = c(c(0:a2),(max_IBM_age+0.5)) # a2:max_IBM_age as plus group
}else{
  full_age_breaks = c(c(0:(a2+1)))               # a2 as plus group
}


lower_age_bin = full_age_breaks[-length(full_age_breaks)]

OBS_full_AgeE_COM = NULL   
por_full_Agetrue_COM = NULL
por_full_AgeE_COM = NULL
por_OBS_AgeE_COM = NULL


for(y in unique(CATCH_MATRIX$Year)){
  temp_Year = subset(CATCH_MATRIX, Year==y) ## only works where a fish died in that year
  
  for(k in 1:length(OBS_size_bin-1)){
    temp_conditional_size = subset(temp_Year,temp_Year$fish_size >= OBS_size_bin[k] & temp_Year$fish_size < OBS_size_bin[k+1])
    
    
    if(length(temp_conditional_size[, 1])>=1) {
      
      hist_data_Agetrue = hist(temp_conditional_size$Age,right=F,breaks=full_age_breaks,plot=FALSE) 
      hist_data_AgeE = hist(temp_conditional_size$AgeE,right=F,breaks=full_age_breaks,plot=FALSE)
      
      temp_data_Agetrue = hist_data_Agetrue$counts/sum(hist_data_Agetrue$counts) 
      temp_data_AgeE = hist_data_AgeE$counts/sum(hist_data_AgeE$counts)
      
      por_full_Agetrue_COM = rbind(por_full_Agetrue_COM,temp_data_Agetrue)
      por_full_AgeE_COM = rbind(por_full_AgeE_COM,temp_data_AgeE)
      
      obs_Age_COM = NULL
      temp_por_obs_Age_COM = NULL
      for(i in 1:length(temp_data_AgeE)){
        n_count = 0
        for(j in 1:Effect_N_Age_COM){
          r5=runif(1, min=0, max=1)
          if(r5<=temp_data_AgeE[i]){
            n_count=n_count+1
          }
        }
        obs_Age_COM[i]=n_count
        temp_por_obs_Age_COM[i]=n_count/Effect_N_Age_COM
      }
      por_OBS_AgeE_COM = rbind(por_OBS_AgeE_COM,temp_por_obs_Age_COM)
      
      #Yr","Seas","Flt","Gender","Part", "Ageerr", "Lbin_lo", "Lbin_hi", "Nsamp",           "obs_Age_COM"
      OBS_full_AgeE_COM = rbind(OBS_full_AgeE_COM,c(y,   1,      1,    0,      0,      1,         OBS_size_bin[k],         OBS_size_bin[k],        Effect_N_Age_COM, obs_Age_COM))
      
    } # end of if 
    
  } # end of OBS_size_bin
  
} # end of year

colnames(por_full_Agetrue_COM) <-  c(hist_data_AgeE$breaks[-length(hist_data_AgeE$breaks)]) 
colnames(por_full_AgeE_COM) <-  c(hist_data_AgeE$breaks[-length(hist_data_AgeE$breaks)])   
colnames(por_OBS_AgeE_COM) <-  c(hist_data_AgeE$breaks[-length(hist_data_AgeE$breaks)])   
colnames(OBS_full_AgeE_COM) <-  c("#Yr","Seas","Flt","Gender","Part","Ageerr","Lbin_lo","Lbin_hi","Nsamp",hist_data_AgeE$breaks[-length(hist_data_AgeE$breaks)])

OBS_full_AgeE_COM = as.data.frame(OBS_full_AgeE_COM)
por_full_Agetrue_COM = as.data.frame(por_full_Agetrue_COM)
por_full_AgeE_COM = as.data.frame(por_full_AgeE_COM)

OBS_AgeE_COM = subset(OBS_full_AgeE_COM, select = c("#Yr","Seas","Flt","Gender","Part","Ageerr","Lbin_lo","Lbin_hi","Nsamp",OBS_age_bin))
por_full_Agetrue_COM = subset(por_full_Agetrue_COM, select = paste(OBS_age_bin))
por_full_AgeE_COM = subset(por_full_AgeE_COM, select = paste(OBS_age_bin))
por_OBS_AgeE_COM = subset(por_OBS_AgeE_COM, select = paste(OBS_age_bin))


# Number of OBS_age_bin (age' bin)
N_obs_age_bin = length(OBS_age_bin)
write.table("#_N_age_bins",paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(t(N_obs_age_bin),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(t(OBS_age_bin),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 

# Ageing error
write.table("#_N_ageerror_definitions",paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  
write.table("1",paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  
# vector of mean age' for each true age  
write.table(t(mean_age_true_age[1:(a2+1)]),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  
# vector of SD of age' for each true age
write.table(t(mean_SD_age_true_age[1:(a2+1)]),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 

N_age_observation = length(OBS_AgeE_COM[,1])
write.table("#_N_Agecomp_obs",paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(t(N_age_observation),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 

Age_COM_length_bin = 1
write.table("#_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths",paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  
write.table(t(Age_COM_length_bin),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  

Combine_M_into_F = 1
write.table("#_combine males into females at or below this bin number",paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  
write.table(t(Combine_M_into_F),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T)  

# data-vector
write.table(t(colnames(OBS_AgeE_COM)),paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(OBS_AgeE_COM,paste0(path_name, "/SS3_OBS_AgeE_COM.csv",sep=""),row.names=FALSE)
# print(paste0('wrote SS3 OBS age comp, sim ',sim_number, " FLEVEL ",level))
## SS3 NUMBER/FISHERY BASED CPUE ----

Exp_POP = NULL
for(i in start_yr:(start_yr+simu_year-1)){
  temp_Year = subset(IBM_EXP_MATRIX,IBM_EXP_MATRIX$Year==i)
  EXP_temp = data.frame(Year=c(i),Exp_number=c(length(temp_Year$fish_size)),Exp_biomass=c(sum(lw(temp_Year$fish_size))))
  Exp_POP = rbind(Exp_POP,EXP_temp)
}

CPUE_error = rnorm(length(Exp_POP[,1]),0,sigma_CPUE)                 
OBS_CPUE = data.frame(Year=c(Exp_POP[,1]),season=c(1),fleet=c(2),CPUE=c(Exp_POP[,2]*fishery_q*exp(CPUE_error-0.5*sigma_CPUE^2)),sd=c(sigma_CPUE))

write.table(OBS_CPUE,paste0(path_name, "/SS3_Fishery_CPUE.csv",sep=""),row.names = FALSE,col.names = TRUE,sep=",") 

## SS3 CATCH DATA----

F_MATRIX = subset(IBM_DEAD_MATRIX,IBM_DEAD_MATRIX$dtype=="F")
M_MATRIX = subset(IBM_DEAD_MATRIX,IBM_DEAD_MATRIX$dtype=="M")

CATCH_wt = NULL
CATCH = NULL
DEATH_wt = NULL 

for(i in start_yr:(start_yr+simu_year-1)){
  temp_Year = subset(F_MATRIX,F_MATRIX$Year==i)
  temp_Year_M = subset(M_MATRIX,M_MATRIX$Year==i)
  CATCH = rbind(CATCH,length(temp_Year[,1]))
  CATCH_wt = rbind(CATCH_wt,sum(lw(temp_Year$fish_size)))
  DEATH_wt = rbind(DEATH_wt,sum(lw(temp_Year_M$fish_size)))
}
catch_out = data.frame(CATCH_num=c(CATCH),CATCH_Kg=c(CATCH_wt),CATCH_ton=c(CATCH_wt)/1000,Year=c(start_yr:(start_yr+simu_year-1)),season=c(1))
write.table(catch_out,paste0(path_name, "/SS3_Catch.csv",sep=""),row.names = FALSE,col.names = TRUE,sep=",") 
# print(paste0('wrote SS3 catch data, sim ',sim_number," FLEVEL ",level))

## SS3 MEAN LENGTH AT AGE ----

# file.remove(paste0(path_name, "/SS3_OBS_Mean_Size_At_AgeE.csv",sep=""))

OBS_Mean_Length_At_Age_RESULT = NULL

for(y in (start_yr+M_only_yr):(start_yr+simu_year-1)){
  
  temp_Year = subset(CATCH_MATRIX,CATCH_MATRIX$Year==y)
  
  temp_mean_size_at_age = NULL
  RESULT_fish_N = NULL
  for(age in 1:(length(full_age_breaks)-1)){
    temp_age = subset(temp_Year,temp_Year$AgeE >= full_age_breaks[age] & temp_Year$AgeE < full_age_breaks[age+1])
    if(length(temp_age[,1])>=Nsamp_Mean_length_at_age){
      index_all = c(1:length(temp_age[,1]))
      sample_index = sample(index_all,Nsamp_Mean_length_at_age,replace = FALSE)
      sample_temp_age = temp_age[sample_index,]
      mean_size_at_age = mean(sample_temp_age$fish_size)
      fish_N = Nsamp_Mean_length_at_age
      temp_mean_size_at_age=rbind(temp_mean_size_at_age,mean_size_at_age)
      RESULT_fish_N = rbind(RESULT_fish_N,fish_N)
    }
    else if(length(temp_age[,1])==0){
      mean_size_at_age = c(-999)
      fish_N = 0 
      temp_mean_size_at_age = rbind(temp_mean_size_at_age,mean_size_at_age) 
      RESULT_fish_N = rbind(RESULT_fish_N,fish_N)
    }
    else{
      mean_size_at_age = mean(temp_age$fish_size)
      fish_N = length(temp_age[,1]) 
      temp_mean_size_at_age = rbind(temp_mean_size_at_age,mean_size_at_age) 
      RESULT_fish_N = rbind(RESULT_fish_N,fish_N) 
    } 
  } # age loop
  #Yr","Seas","Flt","Gender","Part", "Ageerr", "Nsamp", "Mean-size-at-age for each age* " , "N_fish"
  OBS_Mean_Length_At_Age_RESULT = rbind(OBS_Mean_Length_At_Age_RESULT,c(y,   1,      1,    0,      0,      1,        2,     temp_mean_size_at_age,            RESULT_fish_N))
  
} # year loop

colnames(OBS_Mean_Length_At_Age_RESULT) <- c("#Yr","Seas","Flt","Gender","Part","Ageerr","Nsamp",lower_age_bin, paste("N_fish_a",lower_age_bin,sep=""))
OBS_Mean_Length_At_Age_RESULT = as.data.frame(OBS_Mean_Length_At_Age_RESULT)  
OBS_Mean_Length_At_Age_RESULT = subset(OBS_Mean_Length_At_Age_RESULT, select = c("#Yr","Seas","Flt","Gender","Part","Ageerr","Nsamp",OBS_age_bin, paste("N_fish_a",OBS_age_bin,sep="")))

# Number of OBS_age_bin (age' bin)
N_obs_mean_size_at_age = length(OBS_Mean_Length_At_Age_RESULT[,1])
write.table("#_N_MeanSize-at-Age_obs",paste0(path_name, "/SS3_OBS_Mean_Size_At_AgeE.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(t(N_obs_mean_size_at_age),paste0(path_name, "/SS3_OBS_Mean_Size_At_AgeE.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 

# data-vector
write.table(t(colnames(OBS_Mean_Length_At_Age_RESULT)),paste0(path_name, "/SS3_OBS_Mean_Size_At_AgeE.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
write.table(OBS_Mean_Length_At_Age_RESULT,paste0(path_name, "/SS3_OBS_Mean_Size_At_AgeE.csv",sep=""),row.names=FALSE,col.names=FALSE,sep=",",append = T) 
# print(paste0('wrote SS3 mean size at age, sim ',sim_number," FLEVEL ",level))