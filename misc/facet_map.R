## Figure X facet map with hypotheses

setwd("C:/Users/mkapur/Dropbox/UW/sab-growth")


## load mappable data, made in len_dataprep


# some standard map packages.
# https://eriqande.github.io/rep-res-eeb-2017/map-making-in-R.html
# the ggmap package.  Might as well get the bleeding edge version from GitHub
# devtools::install_github("dkahle/ggmap")
# install.packages(c("maps", "mapdata"))
require(dplyr);require(maps);require(mapdata);require(ggplot2)
require(Rmisc)
# mapdf <- sample_n(read.csv(paste0(getwd(),"/data/mapdf.csv")),10000)

usa <- map_data("world") # we already did this, but we can do it again
map <- ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
  coord_quickmap() +
  scale_x_continuous(expand = c(0,0), limits = c(-180,-115), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
  scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels =  paste(seq(30,75,10), "°N"))  +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), axis.title =element_blank()) #+

map2 <- ggplot() + 
  geom_polygon(data = usa, aes(x = long, y = lat, group = group)) + 
  coord_quickmap() +
  scale_x_continuous(expand = c(0,0), limits = c(-180,-115), breaks = seq(-180,-120,10), labels = paste(seq(-180,-120,10), "°W")) +
  scale_y_continuous(expand = c(0,0), limits = c(30,75), breaks = seq(30,75,10), labels = paste(seq(30,75,10), "°N")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        axis.title.x.bottom =  element_blank()) 

multiplot(map,map,map,map,cols = 2)


