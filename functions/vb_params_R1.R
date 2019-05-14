## GROWTH MODULE ----
# length-weight
lw_a = 1.3528e-006
lw_b = 3.4297

# SS3 growth model (vonBert with L1&L2, SS3 3.34f mannual, Page 42)
L1 = 62
L2 = 215
a1 = 0
a2 = 15
growth_lognormal_SD = 0.025
VBGF_K = 0.25
# 
# source("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/misc/Predict_Length_Fn.R")
# 
# 
# ldat <- read.csv("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/iPopSim/F0L1S_25/R2/MED_MED_MED/MED_MED_MED1sim_Comp.csv") %>% 
#   filter(REG == 'R2') %>% select(length = Length_cm,age= Age)
# vals <- sample_fit_vbgf(
#   length.data = ldat,
#   start.L1 = L1,
#   start.L2 = L2,
#   start.k = VBGF_K,
#   start.cv.young = 0.02,
#   start.cv.old = 0.11,
#   lo.L1 = L1-2,
#   lo.L2 = L2-2,
#   lo.k = VBGF_K*0.8,
#   lo.cv.young = 0.01,
#   lo.cv.old = 0.03,
#   hi.L1 = L1+2,
#   hi.L2 = L2+2,
#   hi.k = VBGF_K*1.1,
#   hi.cv.young = 0.15,
#   hi.cv.old =0.20,
#   a3 = 0,
#   A =  a2
# )
# 
# vals$CV_young_Fem_GP_1*vals$L_at_Amax_Fem_GP_1 ## cv = sd/mean, multiply to get sd. note labs are reversed
