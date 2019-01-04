## test script for SAB
##load data
library(VAST);library(TMB)
load("C:/users/mkapur/dropbox/uw/sab-growth/sab_VAST/data/corrected_sab.rda")
run_one_spp(sab, 
                config_file = "vast_config_sab.R",
                folder_name = "sab_spatiot")