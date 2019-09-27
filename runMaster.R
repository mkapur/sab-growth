# runMaster
## run everything...with confidence!

source("IBM_master_reboot.R") ## regen if needed. Builds IBM
cat("################################# \n 
    IBM COMPLETE \n
    ################################## \n")

source("GAM_master.R")
cat("################################# \n 
    GAM COMPLETE \n
    ################################## \n")

source("STARS_master.R")
cat("################################# \n 
    STARS COMPLETE \n
    ################################## \n")


source("n_sensitivity.R")
cat("################################# \n 
    NSENS COMPLETE \n
    ################################## \n")


source("plot_master.R")
cat("################################# \n 
    PLOTS COMPLETE \n
    ################################## \n")