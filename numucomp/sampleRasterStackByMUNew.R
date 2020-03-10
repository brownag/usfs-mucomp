#sampleRasterStackByMUNew
sampleRasterStackByMUNew <- function (mu, mu.set, mu.col, raster.list, pts.per.acre, 
                                      p = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), 
                                      progress = TRUE, estimateEffectiveSampleSize = TRUE) 
{
  if (!requireNamespace("rgdal") | !requireNamespace("rgeos") | 
      !requireNamespace("raster") | !requireNamespace("spdep")) 
    stop("please install the packages: rgdal, rgeos, raster, spdep", 
         call. = FALSE)
  
  if (!is.projected(mu)) 
    stop("map unit polygons must be in a projected CRS", 
         call. = FALSE)
  
  # extract components of raster list
  raster.type <- lapply(raster.list, function(l) return(l$type))
  raster.files <- lapply(raster.list, function(l) return(l$file))
  raster.pnames <- lapply(raster.list, function(l) return(l$name))
  raster.labels <- lapply(raster.list, function(l) return(l$label))
  
  validity.res <- data.frame(id = mu[[mu.col]], 
                             Polygon.Validity = rgeos::gIsValid(mu, byid = TRUE, reason = TRUE), 
                             stringsAsFactors = FALSE)
  l.mu <- list()
  l.unsampled <- list()
  a.mu <- list()
  
  raster.layers <- rapply(raster.files, how = "replace", f = function(i) {
    i <- try(raster::raster(i))
    if (class(i) == "try-error") 
      stop(paste0("Cannot find raster file: ", i), call. = FALSE)
    else return(i)
  })
  
  message("Loading raster data...")
  nm <- unlist(lapply(raster.layers, names))
  raster.layers <- rapply(raster.layers, how = "replace", f = function(r) {
    r.mem <- try(raster::readAll(r), silent = TRUE)
    if (class(r.mem) == "RasterLayer") {
      return(r.mem)
    }
    else {
      return(r)
    }
  })
  e.mu <- as(raster::extent(mu), "SpatialPolygons")
  proj4string(e.mu) <- proj4string(mu)
  
  # raster containment test (rasters cover full polygon extent)
  message("Checking raster/MU extents...")
  raster.containment.test <- rapply(raster.layers, f = function(r) {
    e.r <- as(raster::extent(r), "SpatialPolygons")
    proj4string(e.r) <- proj4string(r)
    e.mu.r <- spTransform(e.mu, CRS(proj4string(e.r)))
    return(rgeos::gContainsProperly(e.r, e.mu.r))
  })
  if (any(!raster.containment.test)) 
    warning("Raster extent does not completly cover MU extent")
  
  # estimate effective sampling size? (account for global autocorrlation)
  if (estimateEffectiveSampleSize) {
    message("Estimating effective sample size...")
    MI <- ldply(raster.list$continuous, Moran_I_ByRaster, 
                mu.extent = e.mu, .progress = ifelse(progress, "text", 
                                                     NULL))
    names(MI) <- c("Variable", "Moran.I")
  } else {
    MI <- data.frame(Variable = nm, Moran.I = 0)
  }
  
  message("Sampling polygons, and extracting raster values...")
  if (progress) 
    pb <- txtProgressBar(min = 0, max = length(mu.set), 
                         style = 3)
  for (mu.i in mu.set) {
    mu.i.sp <- mu[which(mu[[mu.col]] == mu.i), ]
    suppressMessages(s <- constantDensitySampling(mu.i.sp, 
                                                  n.pts.per.ac = pts.per.acre, 
                                                  min.samples = 1, 
                                                  polygon.id = "pID", 
                                                  iterations = 10))
    l.unsampled[[mu.i]] <- setdiff(mu.i.sp$pID, unique(s$pID))
    
    if (is.null(s)) {
      a.mu[[mu.i]] <- NULL
      l.mu[[mu.i]] <- NULL
      if (progress) 
        setTxtProgressBar(pb, match(mu.i, mu.set))
      next
    } else {
      s$sid <- 1:nrow(s)
      l.mu[[mu.i]] <- rapply(raster.layers, how = "replace", 
                             f = function(r) {
                               res <- data.frame(value = raster::extract(r, s), 
                                                 pID = s$pID, sid = s$sid)
                               return(res)
                             })
      a <- sapply(slot(mu.i.sp, "polygons"), slot, "area") * 0.000247
      
      .quantiles <- quantile(a, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
      .total.area <- sum(a)
      .samples <- nrow(s)
      .mean.sample.density <- round(.samples/.total.area,  2)
      .polygons <- length(a)
      .unsampled.polygons <- length(l.unsampled[[mu.i]])
      
      a.stats <- c(round(c(.quantiles, .total.area, .samples, 
                           .polygons, .unsampled.polygons)), .mean.sample.density)
      names(a.stats) <- c("Min", "Q5", "Q25", "Median", 
                          "Q75", "Q95", "Max", "Total Area", "Samples", 
                          "Polygons", "Polygons Not Sampled", "Mean Sample Dens.")
      a.mu[[mu.i]] <- a.stats
      
      if (progress) 
        setTxtProgressBar(pb, match(mu.i, mu.set))
    }
  }
  
  if (progress) 
    close(pb)
  
  mu.area <- do.call('rbind', a.mu)
  df.mu <- data.frame(musym=rownames(mu.area))
  colnames(df.mu) <- mu.col
  mu.area <- cbind(df.mu, mu.area)
  
  d.mu <- lapply(l.mu, function(i) {
    df.sid <- i[[1]][,c("pID","sid")]
    res <- cbind(df.sid, as.data.frame(do.call('cbind', lapply(i, function(j) return(j$value)))))
    colnames(res) <- c("pID","sid",raster.pnames)
    return(res)
  })
  d.mu.musym <- do.call('c', lapply(1:length(d.mu), function(i) return(rep(names(d.mu)[i], nrow(d.mu[[i]])))))
  
  d.mu <- do.call('rbind', d.mu)
  df.musym <- data.frame(mysym = d.mu.musym)
  colnames(df.musym) <- mu.col
  d.mu <- cbind(df.musym, d.mu)
  rownames(d.mu) <- NULL
  
  unsampled.idx <- unlist(l.unsampled)
  rs <- rapply(raster.layers, f = raster::filename, how = "unlist")
  rs <- gsub("\\\\", "/", rs)
  
  rs.df <- data.frame(Variable = nm, File = rs, 
                      inMemory = as.character(rapply(raster.layers, f = raster::inMemory, how = "unlist")), 
                      ContainsMU = raster.containment.test)
  
  rs.df <- merge(rs.df, MI, by = "Variable", sort=FALSE)
  
  rs.df$Moran.I[is.na(rs.df$Moran.I)] <- 0
  return(list(raster.samples = d.mu, area.stats = mu.area, 
              unsampled.ids = unsampled.idx, raster.summary = rs.df, 
              mu.validity.check = validity.res, Moran_I = MI))
}
