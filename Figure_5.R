## Figure 5
# Out of desperation.

install.packages('oceanmap')
library(oceanmap)
# https://rstudio-pubs-static.s3.amazonaws.com/301644_61297506548c4b3cb9e7cc8cbc8578ce.html
# https://www.arcgis.com/home/item.html?id=24bfd85e97b042948e6ed4928dc45a8b
devtools::install_github("yonghah/esri2sf")
library("esri2sf")
url <- "http://services1.arcgis.com/VAI453sU9tG9rSmh/arcgis/rest/services/WorldGeo_Physical_Climate_features/FeatureServer"
df <- esri2sf(url)
plot(df)
# https://data.amerigeoss.org/en_AU/dataset/major-ocean-currents-arrowpolys-30m-85


require(rgdal)
shape <- readOGR(dsn = "C:/Users/Maia Kapur/Dropbox/UW/sab-growth/raw_data/Major_Ocean_Currents_arrowPolys_30m_8", 
                 layer = 'Major_Ocean_Currents_arrowPolys_30m_8')

png("C:/Users/Maia Kapur/Dropbox/UW/sab-growth/figures/Figure5.png",height = 10, width = 8, unit = 'in', res = 520)
lon <- c (-180, -110)
lat <- c (26, 74)
# figure ( width =9.75, height =5.28)
plotmap ( lon=lon ,  lat=lat , main ="Northeast Pacific", grid = T, col.bg = 'grey90', col.land = 'white')
grid()
segments(lon[1], 50, x1 = lon[2], y1 = 50, col = 'red', lty = 'dashed')
segments(lon[1], 36, x1 = lon[2], y1 = 36, col = 'red', lty = 'dashed')
segments(-130, lat[1], x1 = -130, y1 = lat[2], col = 'red', lty = 'dashed')
segments(-145, lat[1], x1 = -145, y1 = lat[2], col = 'blue', lty = 'dashed')

rect(xleft = lon[1], ybottom = lat[1], xright = -130, ytop = 49.5, col = 'grey88', border = NA, add = T)
text(x = c(-155,-135,-117,-117,-115), y = c(65,65,55,45,35), labels = paste('Region',5:1))
plot(shape,  ylim = lat, col = c( "gold", "dodgerblue","#ABDDA4", "#ABDDA4" ,"grey22"), border = 'white', add = T)
legend(x = -175, y = 40, legend = c('Alaskan Current','N. Pacific Current','California Current', 'S. California Bight', 'GAM-detected Regions', 'Ecosystem break'),
       col = c('dodgerblue','#ABDDA4','grey22','gold','red','blue'), lty = c(rep(1,4),2,2), lwd = c(rep(5,4),1,1))

dev.off()