## Generic GAM & growth estimation script
## M Kapur kapurm@uw.edu

require(mgcv);require(dplyr);require(ggplot2); require(TMB); library(reshape2)
library(gridExtra); library(grid); library(lattice)

## source some functions from github
source("https://raw.githubusercontent.com/mkapur/sab-growth/master/functions/Deriv.R")

## load in your data
## columns are Year, Length_cm, Age, Sex, Latitude_dd, Longitude_dd, and REG
load(paste0("./input_data/gam_data_sab_0415.rda")) 
## ensure correct wrapping around dateline
full_data$Longitude_dd[full_data$Longitude_dd > 0] <- full_data$Longitude_dd[full_data$Longitude_dd > 0]*-1


## Fit GAMS to ID breakpoints ----
breaksdf <- list(); idx = 1## storage for breakpoints;
for(AA in c(4,6,30)){ ## ages to define breakpoints; often roughly correspond with L1/L2
  for(SS in c("F","M")){
    dat <- full_data %>% filter( Age == AA  & Sex == SS)
    mod <- gam(Length_cm ~ s(Year, bs = "cc") + s(Latitude_dd) + s(Longitude_dd),data = dat)
    
    # png(paste0("./figures/sab_gam_diagnostics_",AA,SS,".png"), width = 7, height = 5, units = 'in', res = 420)
    layout(matrix(1:4, ncol = 2))
    gam.check(mod)
    # graphics.off()

    ## get & eval derivatives ----
    pdat <- dat
    pTerm <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE) ## predict on original data
    p2 <- predict(mod, newdata = pdat) ## raw predicts
    pdat <- transform(pdat, predLen = p2, se2_lat = pTerm$se.fit[,2],se2_lon = pTerm$se.fit[,3], se2_yr = pTerm$se.fit[,1])

    df.res <- df.residual(mod)
    crit.t <- qt(0.025, df.res, lower.tail = FALSE)
    ## variances are additive just FYI
    pdat <- transform(pdat,
                      upper = predLen + (crit.t * (se2_lat+se2_lon+se2_yr)),
                      lower = predLen - (crit.t * (se2_lat+se2_lon+se2_yr)))


    Terms <- c("Year","Latitude_dd","Longitude_dd") 
    ## predict over parameter space
    for(t in 1:length(Terms)){
      Term <- Terms[t]


      newD <- data.frame( Year = seq(1981,2017,length = 100),
                          Longitude_dd = seq(-186,-116, length = 100),
                          Latitude_dd = seq(0,64,length = 100)) ## whatever isn't modeled will just get ignored


      m2.d <- Deriv(mod, newdata = newD)
      m2.dci <- confint(m2.d, term = Term)

      crit.eval = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$deriv) ## use tails
      crit.eval.se = quantile(probs = c(0.025, 0.975), x =  m2.d[[Term]]$se.deriv) ## use tails


      ## identify where mean crosses zero or falls out of bounds (NAs where it does)
      m2.dsig.zeros <- signifMK(x = m2.d$eval[[Term]],
                                d = m2.d[[Term]]$deriv,
                                upper = m2.dci[[Term]]$upper,
                                lower = m2.dci[[Term]]$lower,
                                crit.eval = crit.eval,
                                eval = 0)

      pix <- !is.na(m2.dsig.zeros)

      vals <- m2.d$eval[[Term]][pix] ## what test vals did these correspond to
      breaksdf[[idx]] <- sort(c(unique(round(vals)))) ## get rounded unique
      idx = idx+1
      ## fill NAs in bdf for binding
      for(i in 1:length(breaksdf)){
        if (length(breaksdf[[i]]) == 0){ ## fill NA for empty
          breaksdf[[i]] <- NA
        }
      }## end breaksdf

    } ## end terms
  } ## end sexes
} ## end key ages

## reformat breaks DF for interpredation
breaks_df <- data.frame(Age = rep(c(4,6,30),2), Sex = c(rep("F",3),rep("M",3)), 
                        Year = NA, Lat = NA, Lon = NA)
vec <- c(seq(1,length(breaksdf),3))
for(i in 1:length(vec)){
  breaks_df$Year[i] <- breaksdf[[vec[i]]]
  breaks_df$Lat[i] <- breaksdf[[vec[i]+1]]
  breaks_df$Lon[i] <- breaksdf[[vec[i]+2]]
}




# Estimate growth parameters at detected breaks ----
## dwnld CPP from https://github.com/mkapur/sab-growth/blob/master/sptlVB_Sel_Sigma.cpp 
compile("sptlVB_Sel_Sigma.cpp") 
dyn.load(dynlib("sptlVB_Sel_Sigma"))