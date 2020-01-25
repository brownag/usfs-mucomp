library(raster)
july <- raster('LC08_L1TP_042034_20180705_20180717_01_T1_sr_ndvi_utm.tif')
febr <- raster('LC08_L1TP_042034_20180211_20180222_01_T1_sr_ndvi_utm.tif')
july <- resample(july, febr, "bilinear")
july <- crop(july, extent(febr))
delta_NDVI <- july - febr
plot(density(getValues(delta_NDVI), na.rm=T))
