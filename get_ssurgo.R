library(rgdal)
library(soilDB)
huc_wc <- rgdal::readOGR("geodata/HUC10_WillowCreek.shp")

# query mukeys using HUC boundary
res <- SDA_spatialQuery(huc_wc)

# fetch full extent of mukeys used within HUC
# add musym and muname to base result
s <- fetchSDA_spatial(res$mukey, chunk.size = 1, 
                          add.fields = c('mapunit.musym',
                                         'mapunit.muname'))

# convert to projected CRS (boundary is in UTM Zone 11)
s.t <- spTransform(s, CRS(proj4string(huc_wc)))

# display SSURGO lines as backdrop
plot(s.t, col='blue')

# plot watershed boundary on top
plot(huc_wc, lwd=2, border='red', add=T)

# inspect SPDF
s.t

# write result to shapefile in geodata folder
writeOGR(s.t, dsn="geodata", layer="SSURGO", 
         driver="ESRI Shapefile", overwrite_layer = TRUE)

