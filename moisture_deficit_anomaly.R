# moisture deficit anomaly for sierra national forest - 1981-2015
# @author: andrew brown; 2020/01/09

# load rgdal and raster libraries
library(rgdal)
library(raster)

# load multiband raster for 1981-2015 as a RasterStack
deficit.ts <- stack('def_SNF_1981-2015.tif')

# load boundary shapefile
boundary.shp <- readOGR(dsn=".", layer = "HUC10_WillowCreek")

# project boundary shapefile to CRS of raster
boundary.shp <- spTransform(boundary.shp, CRS(proj4string(deficit.ts)))
            
# define year intervals of interest as numeric vectors
test.period <- 2012:2015
full.period <- 1981:2015

# give nice names to each year/data column in stack
names(deficit.ts) <- paste0("def_SNF_", full.period)

# inspect
deficit.ts

# view plots of individual years
#plot(deficit.ts)

# use names to subset the test period
deficit.ts.test <- deficit.ts[[paste0("def_SNF_", test.period)]]

# calculate mean (all layers assigned index 1) from a stack of 35 
r.baseline <- stackApply(deficit.ts, indices=1, fun=mean)

# calculate mean (all layers assigned index 1) from subset (test) stack 
r.test <- stackApply(deficit.ts.test, indices=1, fun=mean)

r.anomaly <- r.baseline - r.test

writeRaster(r.anomaly, filename = "snf_mdeficit_anomaly_1981-2015.tif")

# make a plot of full extent(write to PNG file)
png('snf_mdeficit_anomaly_1981-2015.png')
plot(r.anomaly, xlab="Longitude", ylab="Latitude",
     main=paste0("Sierra National Forest - Drought Anomaly\n Mean Moisture Deficit ",
     min(test.period),"-",max(test.period)," versus ",min(full.period),"-",max(full.period)))
# add project area  lines to plot
lines(boundary.shp)
dev.off()

# crop r.anomaly to project boundary
r.anomaly.crop <- crop(r.anomaly, boundary.shp)

# mask out portions within extent rectangle, but outside boundary
r.anomaly.crop <- mask(r.anomaly.crop, mask = rasterize(boundary.shp, r.anomaly.crop))

# make a plot of project extent(write to PNG file)
png('WillowCreek_mdeficit_anomaly_1981-2015.png')
plot(r.anomaly.crop, xlab="Longitude", ylab="Latitude",
     main=paste0("Sierra National Forest - Drought Anomaly\nMean Moisture Deficit ",
                 min(test.period),"-",max(test.period)," versus ",min(full.period),"-",max(full.period)))

# convert project area boundary to r.anomaly CRS and add lines to plot
lines(spTransform(boundary.shp, CRS(proj4string(r.anomaly.crop))))
dev.off()
