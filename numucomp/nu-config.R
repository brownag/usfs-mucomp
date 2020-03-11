# config.R
# nu-mucomp report

# new list format for rasters
raster.list <- list(
  # DEM meters
    list(label="Elevation (m)", 
         name="elevation", 
         file='../geodata/WillowCreek_elev.tif', 
         type="continuous", 
         precision=0),
  # Slope Gradient %
    list(label="Slope (%)", 
         name="slope", 
         file='../geodata/WillowCreek_gradient.tif', 
         type="continuous",
         precision=0)
)

threshold.list <- list(
    list(pattern='^13[67]$', name="elevation", vmin=1095, vmax=1770),
    list(pattern='^13[8]$', name="elevation", vmin=820, vmax=1675),
    list(pattern='^13[68]$', name="slope", vmin=5, vmax=35),
    list(pattern='^137$', name="slope", vmin=65, vmax=100)
)

# path to parent folder of SHP, no trailing forward slash (/)
mu.dsn <- '../geodata'
# SHP name, without file extension
mu.layer <- 'SSURGO'

############################################
### column with map unit ID / key / symbol #
############################################

# could be 'MUKEY', 'MUSYM', or any valid column name
mu.col <- 'musym'
mu.set <- c('136','137','138','139','140','141','142')

#########################################################
### polygon sampling density (samples / acre / polygon) #
#########################################################

# consider using a sampling density between 1-2 points / ac.
# increase if there are un-sampled polygons
# delineations smaller than 5 ac. may require up to 5 points / ac.
# values > 6-7 points / ac. will only slow things down
pts.per.acre <- 0.1

###########################
### quantiles of interest #
###########################

# the most important quantiles (percentiles / 100) are: 0.1, 0.5 (median), and 0.9
# optionally reduce the number of quantiles for narrower tables
p.quantiles <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)

# cache samples?
cache.samples <- TRUE

# correct sample size?
correct.sample.size <- FALSE
