## viz of some distributions, using WC data
# https://www.statmethods.net/advgraphs/probability.html
require(dplyr)
require(ggplot2)

load(paste0(getwd(),"/lenarray.rda")) ## length data

lendf<- lenmat[,,1][!is.na(lenmat[,,1])]
lendf <- data.frame(as.vector(lendf))
names(lendf) <- 'len'
ggplot(lendf, aes(x = len)) + 
  theme_bw() +
  # geom_histogram(fill = alpha('red',0.2)) +
  geom_density(aes(x = rnorm(32728,54.54,2)), fill = alpha('blue',0.2)) +
  geom_vline(xintercept = 54.54) +
  geom_density(aes(x = rnorm(32728,54.54,2)), fill = alpha('blue',0.2)) +
  
  
mean(lendf)

dnorm(63,54.54,0.01)
