## Build and run VAST Analysis on LENGTH data for Sablefish
## M Kapur sourced from C Stawitz and J Thorson Winter 2018
## kapurm@uw.edu
library(VAST); library(compiler); library(dplyr); library(TMB)
Version = "VAST_v5_3_0"
mapply(source, list.files(paste0(getwd(),"/R/"), pattern = ".R", full.names=TRUE))

#Vector to rename columns to
renames <- c('Year',
             'Lat','Lon',
             'length', 'depth', 'temp')
#Build Data_Geostat
load("./data/vast_data.rda")

vast_data <-
  vast_data %>% mutate(SPECIES_CODE = 20510,
                      GEAR_DEPTH = NA,
                      Temp = NA,
                      GEAR_TEMPERATURE = NA) %>%
  filter(Longitude_dd <0) %>%
  plyr::rename(
    c(
      'Length_cm' = 'LENGTH..cm.',
      'Age' = 'AGE',
      'Latitude_dd' = 'START_LATITUDE',
      'Longitude_dd' = 'START_LONGITUDE',
      'Year' = 'YEAR'
    )
  )

# with(vast_data, plot(START_LATITUDE ~ START_LONGITUDE, pch = 19))
p <- proc.time() ## about 9 mins
sab <- build_corrected_df(
  vast_data,
  species_code = 20510,
  sex = 'F',
  age = "4",
  renames
)
proc.time() - p
save(sab, file = "./data/corrected_sab.rda")
run_one_spp(sab, 
            config_file = "vast_config_sab.R",
            folder_name = "sab_spatiot")

run_one_sppOG(pcod, 
            config_file = "vast_config_pcod.R",
            folder_name = "pcod_spatiot")

pos <- filter(pcod, Lon<0)
for(i in unique(pcod$Year)){
  png(paste("pcod",i,".png"))
plot(Lat~Lon, cex = 0.01, main = i, data=filter(pos, Year==i))
map("world", add=T)
dev.off()
}


#Format covariates for pollock model
covsperknot <- suppressMessages(FishStatsUtils::format_covariates(
  Lat_e = yellowfin$Lat,
  Lon_e = yellowfin$Lon,
  t_e = yellowfin$Year,
  Cov_ep = yellowfin[,"depth"],
  Extrapolation_List = Extrapolation_List,
  Spatial_List = Spatial_List,
  FUN = mean,
  Year_Set = sort(unique(yellowfin$Year)),
  na.omit = "time-average"))
#X_xtp <- array(data=NA, dim=c(100,33,2))
X_xtp <- apply(covsperknot$Cov_xtp, 2:3, scale)
#X_xtp[, , 2] <- scale(exp(Dens_xt))
dimnames(X_xtp)[[1]] <- dimnames(covsperknot$Cov_xtp)[[1]]

run_one_spp(Data_Geostat, config_file="vast_config_pollock",
            folder_name="Pollock_Spatial",
            covars=X_xtp)


pcod = build_corrected_df(clean_data, species_code =21720, 
                          sex=2, age="4", renames)
pcod = filter(pcod, Year>=1988)

#Format covariates for pollock model

run_one_spp(Data_Geostat=pcod, config_file="vast_config_pcod.R",
            folder_name="Pcod_Spatiotemp_Depth",
            covar_columns=c("depth"))


bestmods <- c("pcod_spatiotemp", "VAST_output_pollock_spatiotemp_depth", "yellowfin_spatio_depth", "nsole_spatio_depth")
filedirs<-paste0(here(),"/",bestmods, "/parameter_estimates.RData")

sd_vect <- paramest <- ts <- vector("list")
ind2 <- c(32, 39, 39, 25)
ind1 <- c(3,4,4,4)
for(i in 1:4){
 load(filedirs[i])
  sd_vect[[i]] <- parameter_estimates$SD
  paramest[[i]] <- parameter_estimates$par
  ts[[i]] <- paramest[[i]][ind1[i]:ind2[i]]
}

results_tx<-rbind(c(rep(NA,6),
                    unname(ts[[1]])),
                  unname(ts[[2]]),
                  unname(ts[[3]]),
                  c(rep(NA,14),
                    unname(ts[[4]])))

png("pollockvyellowfin.png")
plot(results_tx[2,]~years, type="l", ylim=c(3,4.1))
lines(results_tx[3,]~years, col="blue")
dev.off()

strsplit(outputs[[1]]$X.par, " ")
unlist(outputs[[2]]$X.par)[1:18]
unlist(outputs[[3]]$X.par)[1:22]
unlist(outputs[[4]]$X.par)[1:12]
?read.delim
