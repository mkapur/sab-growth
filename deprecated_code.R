## deprecated code


## a function to get 500 samples within selectivity bounds
# data.sub <- function(len, s){
#   ## reshape raw length and age data
#   agemat0 <-
#     len %>% select(Age, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID'),'Age') %>% 
#     reshape::cast(., ID + value ~ st, fill = NA, drop= T) %>% select(-ID, -value) 
#   lenmat0 <- 
#     len %>% select(Length_cm, st) %>% mutate(ID = 1:nrow(.)) %>% melt(., id = c('st', 'ID'),'Length_cm') %>% 
#     reshape::cast(., ID + value ~ st, fill = NA, drop= T) %>% select(-ID, -value) 
#   
#   agemat <- lenmat <- age_unif <- len_unif <-  matrix(NA, nrow = 500, ncol = length(unique(len$st)))
#   
#   
# 
#   for(i in 1:nStrata){
#     temp0 <- matrix(na.omit(lenmat0[,i]),ncol = 1) ##if there is a value in one, should be in the other as well 
#     
#     if(s == 1) temp <- temp0 
#     if(s == 2) { temp <- temp0[temp0 >= minSel[i]]}
#     # if(s == 3) { temp <- temp[temp <= maxSel[i]]}  ## sample size issues
#     # if(s == 4) { temp <- temp[temp >= minSel[i] & temp <= maxSel[i]]} ## sample size issues
#     idx <- sample(1:length(temp),500) ## pick 500 rows and store index
#     lenmat[,i] <- temp[idx] %>% as.matrix()
#     agemat[,i] <- data.frame(na.omit(agemat0[,i]))[idx,] %>% as.matrix() ## select same rows
#     
#     len_unif[,i] <- temp0[idx] %>% as.matrix() ## store raw data (used in all sims)
#     age_unif[,i] <- data.frame(na.omit(agemat0[,i]))[idx,] %>% as.matrix()  
#     
#     rm(temp0)
#   }
#   return(list(lenmat,agemat,len_unif,age_unif))
# }