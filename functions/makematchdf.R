# scens <- unique(scenarios$DESC)
# vars <- paste0('R',1:3)
# vars2 <- paste0('R',1:4)
# vis <- c("pooled",'early','late')
# parms <- c('L1','L2')
# 
# matchdf <- data.frame(expand.grid(scens,vars,vars2,vis,parms))
# names(matchdf) <- c("nscen","trueREG","guessREG","guessCREG","variable")
# matchdf$value[matchdf$trueREG == 'R1' & matchdf$variable == 'L1'] <- 62.69
# matchdf$value[matchdf$trueREG == 'R1'& matchdf$variable == 'L2'] <- 216.72
# 
# matchdf$value[matchdf$trueREG == 'R2' & matchdf$variable == 'L1'] <- 50
# matchdf$value[matchdf$trueREG == 'R2'& matchdf$variable == 'L2'] <- 350
# 
# matchdf$value[matchdf$trueREG == 'R3' & matchdf$variable == 'L1'] <- 62.69
# matchdf$value[matchdf$trueREG == 'R3'& matchdf$variable == 'L2'] <- 248
# matchdf<-matchdf[!duplicated(matchdf),]
# ## FORCE DELETE
# write.csv(matchdf,"./input_data/matchdf3.csv")
# 
# 
# 
# 
