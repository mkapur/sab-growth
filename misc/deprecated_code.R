## deprecated code

Type unif_obj_fun = 0.0; // used as numerator in all simulations
Type trunc_fun = 0.0; // denominator

for(int s = 0; s < 2; s++) // iterate sexes
for(int j = 0; j < nStrata; j++){     // iterate strata columns
  for(int i = 0; i < nmat(j,s); i++){ // iterate rows, unique to strata
    int idx = (s)*nStrata + j;
    yfit = Linf(idx)*(1-exp(-k(idx)*(Age(i,j,s) - t0(idx))));
    Sigma = sigma0*pow(yfit,sigma1) // Francis 1988
    unif_obj_fun = dnorm(Length_cm(i,j,s), yfit, Sigma, false);
    
    switch( selType ){
      case 1 : // uniform selectivity (Sel = 1)
      trunc_fun = 1*pnorm(yfit, Length_cm(i,j,s), Sigma);
      obj_fun -= log( (Sel(i,j,s)*unif_obj_fun)/trunc_fun );
      ypreds(i,j,s) = yfit;
      break;
      case 2 : // selectivity correction
      trunc_fun = Sel(i,j,s)*pnorm(yfit, Length_cm(i,j,s), Sigma);
      obj_fun -= log( (Sel(i,j,s)*unif_obj_fun)/trunc_fun );
      
      ypreds(i,j,s) = yfit;
      break;
      // case 3 : // no correction
      //   obj_fun -= log(unif_obj_fun);
      //   ypreds(i,j,s) = yfit;
      //   break;
      // case 4 : // selectivity co
      //   trunc_fun = pnorm(maxSel(j), yfit, Sigma) - pnorm(minSel(j), yfit, Sigma) ;
      //   obj_fun = unif_obj_fun/log(trunc_fun);
      //   ypreds(i,j,s) = yfit;
      //   break;
      
    } // end selType
  } // end rows
}

REPORT(ypreds);
REPORT(yfit);
REPORT(unif_obj_fun);
REPORT(trunc_fun);
ADREPORT(Linf);
ADREPORT(k);
ADREPORT(t0);


## trunc fun switch
# switch( selType ){
#   case 1 : // uniform selectivity (Sel = 1)
#   trunc_fun = 1*pnorm(yfit, Length_cm(i,j,s), Sigma);
#   obj_fun -= log( (Sel(i,j,s)*unif_obj_fun)/trunc_fun );
#   ypreds(i,j,s) = yfit;
#   break;
#   case 2 : // selectivity correction
#   trunc_fun = Sel(i,j,s)*pnorm(yfit, Length_cm(i,j,s), Sigma);
#   obj_fun -= log( (Sel(i,j,s)*unif_obj_fun)/trunc_fun );
#   
#   ypreds(i,j,s) = yfit;
#   break;

# ggplot(subset(mdf, model == 'Length'), aes(x = age, y = Length)) +
#   theme_minimal() +
#   theme(panel.grid = element_blank(), legend.position = c(0.9,0.15))+
#   # scale_y_continuous(limits = c(0,100)) +
#   # scale_x_continuous(limits = c(0,50)) +
#   geom_point(alpha = 0.2) +
#   geom_line(data = subset(mdf, model != 'Length'), aes(x = age, y = Length, color = model), lwd = 1.1)+
#   facet_wrap(~ st_f) +
#   labs(
#     title = "Predicted Model Fits and Subsampled Data",
#     y = 'Length (cm)',
#     x = 'Age (yr)',
#     subtitle = 'Points are actual subsampled data used in parameter fitting',
#     color = 'selectivity model'
#   )
# ggsave(file = paste0(getwd(),"/plots/dome_uni_fits1.png"), plot = last_plot(), height = 5, width = 7, unit = 'in', dpi = 520)


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