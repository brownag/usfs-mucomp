library(knitr, quietly=TRUE)
library(soilReports, quietly=TRUE)

packz <- c("MASS","rgdal","rgeos","raster","plyr","latticeExtra","cluster", 
           "clhs","randomForest","spdep","reshape2","aqp","soilDB","sharpshootR")
newpackz <- packz[!(packz %in% installed.packages()[,"Package"])]

loaded <- lapply(packz, FUN=require, character.only=TRUE, quietly=TRUE)
if(sum(as.numeric(loaded)) != length(packz)) {
  stop("Failed to load one or more required packages! Be sure you have the latest version of the soilReports package from GitHub. Then run `soilReports::reportSetup('region2/mu-comparison')` to install all required packages. Use reportUpdate() to ensure you have the latest version of report.", call. = FALSE)
  geterrmessage()
}

source('custom.R')
source('config.R')

report.extent <- try(readOGR(dsn=extent.dsn, layer=extent.layer, stringsAsFactors=FALSE))
report.legend <- soilDB::SDA_query_features(report.extent, id=extent.col)

chunk_SDA_spatial <- function(mukey.list, nchunk = 10) {
  mukey.chunk <- (1:length(mukey.list) %% nchunk) + 1
  s <- NULL
  
  for(i in 1:max(mukey.chunk)) {
    idx <- which(mukey.chunk == i)
    
    q <- paste0("SELECT G.MupolygonWktWgs84 as geom, mapunit.mukey, mapunit.nationalmusym FROM mapunit 
                CROSS APPLY SDA_Get_MupolygonWktWgs84_from_Mukey(mapunit.mukey) as G 
                WHERE mukey IN ", format_SQL_in_statement(mukey.list[idx]))
    
    #NB: FedData has a different (much simpler, but not equivalent) definition of SDA_query
    #    it also uses the post/rest interface
    sp.res.sub <- suppressMessages(soilDB::SDA_query(q))
    s.sub <- soilDB::processSDA_WKT(sp.res.sub)
    
    if(is.null(s)) {
      s <- s.sub
    } else {
      s <- rbind(s, s.sub)
    }
  }
  return(s)
}
mu.all <- chunk_SDA_spatial(report.legend$mukey)
legall <- suppressMessages(get_mapunit_from_SDA(WHERE=paste0('mukey IN ', format_SQL_in_statement(report.legend$mukey))))
mu.all <- merge(mu.all, legall, by="mukey", sort=F)

mu <- spTransform(mu.all, CRS(proj4string(report.extent))) # skip cropping at this point -- derive samples from full extent of MUs

mu.crop <- crop(mu, report.extent)
writeOGR(mu.crop, dsn='cleanreport', layer = 'WillowCreek_SSURGO', driver="ESRI Shapefile")
if(class(mu) == 'try-error')
  stop(paste0('Cannot read extent polygon/feature file: "', extent.dsn, ' / ', extent.layer, '"'), call. = FALSE)
if(!(mu.col %in% names(mu)))
  stop(paste0('Cannot find mapunit symbol column (',mu.col,') in column names from Soil Data Access result.'), call. = FALSE)

accessible.inputs <- file.access(as.character(unlist(raster.list)), mode = 4) + 1 
if(any( accessible.inputs == 0 )) {
  unreadable.files <- names(which(accessible.inputs == 0))
  stop(paste0("The following input files either do not exist or are unreadable: \n", paste0(unreadable.files, collapse=", ")))
}

mu[[mu.col]] <- as.character(mu[[mu.col]])

if(exists('mu.set')) { 
  # coerce mu.set to character just in integers were specified
  mu.set <- as.character(mu.set)
  
  # check to see if the mu.set specifies musyms that are not in the spatial layer
  mu.set.bak <- mu.set
  mu_nosp <- !(mu.set %in% mu[[mu.col]])
  
  # check if any of the predefined musyms are absent from the spatial
  if(any(mu_nosp))
    mu.set <- mu.set[-which(mu_nosp)] #if so, remove them
  
  # if we removed everything (e.g. no musyms match due to wrong config file?), fail gracefully
  if(length(mu.set) == 0) 
    stop(paste0("Cannot find map unit polygons with symbol: ",paste0(mu.set.bak,collapse=", ")))
  
  # if mu.set defined in config.R, only keep the features with musym matching set
  mu <- mu[which(mu[[mu.col]] %in% mu.set), ] 
} else {
  # if mu.set is not predefined, ordering is determined by the order of musyms in the source file
  mu.set <- sort(unique(mu[[mu.col]]))
}

if(!dir.exists('output')) 
  dir.create('./output')

if(!exists('shp.unsampled.fname')) shp.unsampled.fname <- paste0('un-sampled-', paste(mu.set, collapse='_'))
if(!exists('shp.stats.fname')) shp.stats.fname <- paste0('polygons-with-stats-', paste(mu.set, collapse='_'))
if(!exists('shp.qc.fname')) shp.qc.fname <- paste0('poly-qc-', paste(mu.set, collapse='_'))

if(!exists('csv.qc.fname')) csv.qc.fname <- paste0('poly-qc-', paste(mu.set, collapse='_'),'.csv')
if(!exists('csv.stats.fname')) csv.stats.fname <- paste0('poly-stats-', paste(mu.set, collapse='_'), '.csv')

outputfiles <- c(paste0(c(shp.unsampled.fname,shp.stats.fname,shp.qc.fname),".shp"), csv.qc.fname, csv.stats.fname)
if(any(grepl(basename(outputfiles), pattern='([/\\|<>:\\*?\"])'))) {
  stop("Map unit set or output file name contains invalid characters for filename. Either override default output filenames in config.R or remove [/\\|<>:\\*?\"] characters from your map unit symbols.", call. = FALSE)
}

mu <- mu[, mu.col, drop=FALSE]

mu$pID <- seq(from=1, to=length(mu))

if(cache.samples & file.exists('cached-samples.Rda')) {
  message('Using cached raster samples...')
  .sampling.time <- 'using cached samples'
  load('cached-samples.Rda')
} else {
  .timer.start <- Sys.time()
  sampling.res <- suppressWarnings(sampleRasterStackByMU(mu, mu.set, mu.col, raster.list, pts.per.acre, estimateEffectiveSampleSize = correct.sample.size))
  .timer.stop <- Sys.time()
  
  .sampling.time <- format(difftime(.timer.stop, .timer.start, units='mins'), digits=2)
  print(paste0("Completed sampling of raster stack for poly symbols (",paste(mu.set,collapse=","),") in ",.sampling.time,"."))
  
  sampling.res$raster.samples$.id <- factor(sampling.res$raster.samples$.id, levels=mu.set)
  print("DONE!")
  
  if(cache.samples)
    save(mu, sampling.res,file = 'cached-samples.Rda')
}

## todo: truncate training data by quantiles -- use central 90% for each .id*variable combo

samples <- dcast(sampling.res$raster.samples, .id + pID + sid ~ variable, value.var = "value")
names(samples) <- c(".id", "pID", "sid", "abr", "elev", "gdd",
                    "moistdef", 'ndvi_july', 'nlcd', 'gradient')

# convert moistdeficit to integer
samples$moistdef <- round(samples$moistdef)

# filter training data to samples collected from consociations + Chawanakee-Rockoutcrop cplx,Tollhouse-Rockoutcorp cplx
training.idx <- which(samples$.id %in% c(unique(legall$musym[legall$mukind == "Consociation"]),'126','166'))
train <- samples[training.idx, ]
train <- train[complete.cases(train),]

# omit dams, refactor
train <- train[-which(train$.id == "DAM"),]
train$.id <- factor(train$.id)
train$nlcd <- factor(train$nlcd)

# do mapping of musym (.id) to series name
train.legend <- legall[legall$musym %in% unique(train$.id),]
series.lut <- (stringr::word(train.legend$muname, 1))
names(series.lut) <-  train.legend$musym

# correlation hack to align FS Auberry concept with Holland from Madera
series.lut[names(series.lut) == "HoD"] <- "Auberry"

# correlate Neuns to Chaix (same depth class), lacking geo info
series.lut[series.lut == 'Neuns'] <- "Chaix"

# fix for RO unit
series.lut[names(series.lut) == '147'] <- "Rock outcrop"

train$series <- as.factor(as.character(series.lut[match(train$.id, names(series.lut))]))

library(randomForest)
rf <- randomForest(data=train, series ~ abr+gdd+moistdef+ndvi_july+gradient, ntree=250)
rf

predict.ras <- as.character(unlist(raster.list))

moistdef <- raster(predict.ras[1])
moistdef <- projectRaster(moistdef, crs = CRS("+proj=utm +zone=11 +datum=NAD83"), res=30)
project.extent <- spTransform(readOGR(".", "HUC10_WillowCreek"), CRSobj = CRS(proj4string(moistdef)))
moistdef <- crop(moistdef, project.extent)

gdd <- raster(predict.ras[2])
elev <- raster(predict.ras[3])
gradient <- raster(predict.ras[4])
ndvi_july <- raster(predict.ras[5])
abr <- raster(predict.ras[6])
nlcd <- raster(predict.ras[7])

alignRasters <- function(rasters, categorical) {
  target <- rasters[[1]]
  res <- list(target)
  for(i in 2:length(rasters)) {
    tmp <- crop(rasters[[i]], spTransform(project.extent, CRS(proj4string(rasters[[i]]))))
    if(!categorical[i])
      tmp <- projectRaster(from = rasters[[i]], to = target, method = "bilinear")
    else
      tmp <- projectRaster(from = rasters[[i]], to = target, method = "ngb")
    tmp <- crop(tmp, extent(target))
    res[[i]] <- tmp
    message(paste0("Aligned ", names(rasters[[i]]), " with ", names(target),"."))
  }
  return(res)
}

predictors <- stack(alignRasters(rasters=list(moistdef,abr,elev,gdd,ndvi_july,nlcd,gradient), categorical = c(F,F,F,F,F,T,F)))
names(predictors) <- c("moistdef","abr","elev","gdd","ndvi_july","nlcd","gradient")

predictors <- crop(predictors, extent(project.extent))
predictors <- mask(predictors, mask = fasterize::fasterize(sf::st_as_sf(project.extent), raster = predictors[[1]]))
plot(is.na(predictors))
writeRaster(predictors, filename = "WillowCreek_MUCOMP_predictors.tif", overwrite=T)

the.map <- raster::predict(object=predictors, model=rf, filename="prediction6.tif", overwrite=T)
rasterVis::levelplot(the.map, margin = FALSE)

save(rf, file = "prediction6_rf_model.Rda")

