# Load libraries ---------------------------------------------------------------
library(raster)
library(sp)
library(sf)
library(rgdal)

# Tell R where the data is -----------------------------------------------------
here::here() # make sure you are already in the folder where the data is

# Data cleaning ----------------------------------------------------------------

# Load raster and shapefile data 

accessibility_raster <- raster::raster(here::here("data", "2015_friction_surface_v1.geotiff"))
# The friction surface, which is available here:  http://www.map.ox.ac.uk/accessibility_to_cities/

accessibility_raster

kasungu <- shapefile(here::here("data", "Kasungu_excluding_national_park"))
# Kasungu district polygon having national park clipped out

# Check CRS projection
proj4string(accessibility_raster)

proj4string(kasungu)

# Match projection
kasungu <- spTransform(kasungu, proj4string(accessibility_raster)) 

# Crop and mask raster 
clip_raster <- mask(crop(accessibility_raster, extent(kasungu)), kasungu)

# Check that clipping worked
plot(clip_raster)
plot(kasungu, add=TRUE, lwd=1)

# Save clipped raster
writeRaster(clip_raster, "data/friction_surface_clip.tif", overwrite=TRUE)

# Accessibility Mapping in R ---------------------------------------------------
# 
# Dan Weiss, Malaria Atlas Project, University of Oxford
# 2017-11-06
#
# This script requires the gdistance package (van Etten, J. R Package gdistance: 
# Distances and Routes on Geographical Grids. Journal of Statistical Software 76, 1-21)
#
# This script requires the two user supplied datasets:
# (a) The friction surface, which is available here:  http://www.map.ox.ac.uk/accessibility_to_cities/
# (b) A user-supplied .csv of points (i.e., known geographic coordinates) 
#
# Notes:
# (a) All file paths and names should be changed as needed.
# (b) Important runtime details can be found in the comments.
# (c) This script is suitable only for analyses of moderately sized areas for most (e.g., up to 10 million km^2 in lower latitude settings - GLOBAL RUNS WILL NOT WORK).
#     We recommend using Google Earth Engine for larger areas, with the exception of high-latitude areas where custom approaches are typically required.
#
# Citation: D.J. Weiss, A. Nelson, H.S. Gibson, W. Temperley, S. Peedell, A. Lieber, M. Hancher, E. Poyart, S. Belchior, N. Fullman, B. Mappin, U. Dalrymple, J. Rozier, 
# T.C.D. Lucas, R.E. Howes, L.S. Tusting, S.Y. Kang, E. Cameron, D. Bisanzio, K.E. Battle, S. Bhatt, and P.W. Gething. A global map of travel time to cities to assess 
# inequalities in accessibility in 2015. (2018). Nature. doi:10.1038/nature25181.
# 

## Required Packages
require(gdistance)

# User Defined Variables - used if clipping from the global layer, if no clipping is needed, see lines 54-55 (currently commented out).
# This could also be accomplished by importing a shapefile (for example) 
# Geographic Coordinates (WGS84)
# left   <- -2.0
# right  <- 0.0
# bottom <- 50.0
# top    <- 52.0
transition.matrix.exists.flag <- 0 # if the geo-corrected graph has already been made, this can save time.  Uses the same T.GC.filename as specified using the T.GC.filename variable.

# Input Files
friction.surface.filename <- here::here('data/friction_surface_clip.tif') # Clipped Kasungu friction surface raster

point.filename <- here::here('data/health_facilities_aggregated.csv') #  Use a header.

# Output Files
T.filename <- here::here('data/study.area.T.rds')
T.GC.filename <- here::here('data/study.area.T.GC.rds')
output.filename <- here::here('data/study.area.accessibility_clip.tif')

# Read in the points table
points <- read.csv(file = point.filename) 

points <- points[,c("LONGITU", "LATITUD")] # Just 2 columns.  Structured as [LONGITU, LATITUD] aka [LONG, LAT].

points <- points[!is.na(points$LONGITU),] # remove NA values

# Fetch the number of points
temp <- dim(points)
n.points <- temp[1]

#  Define the spatial template
friction <- raster(friction.surface.filename)
# fs1 <- crop(friction, extent(left, right, bottom, top))
# Use the following line instead of the preceding 2 if clipping is not needed (i.e., to run globally), but be warned that trying this will far exceed the computational capacity available to most users.
fs1 <- raster(friction.surface.filename) 

# Make the graph and the geocorrected version of the graph (or read in the latter).
if (transition.matrix.exists.flag == 1) {
  # Read in the transition matrix object if it has been pre-computed
  T.GC <- readRDS(T.GC.filename)
} else {
  # Make and geocorrect the transition matrix (i.e., the graph)
  T <- transition(fs1, function(x) 1/mean(x), 8) # RAM intensive, can be very slow for large areas
  saveRDS(T, T.filename)
  T.GC <- geoCorrection(T)                    
  saveRDS(T.GC, T.GC.filename)
}

# Convert the points into a matrix
xy.data.frame <- data.frame()
xy.data.frame[1:n.points,1] <- points[,1]
xy.data.frame[1:n.points,2] <- points[,2]
xy.matrix <- as.matrix(xy.data.frame)

# Run the accumulated cost algorithm to make the final output map. This can be quite slow (potentially hours).
temp.raster <- accCost(T.GC, xy.matrix)

# Write the resulting raster
writeRaster(temp.raster, output.filename) 

# Run the r.cost algorithm in QGIS or GRASS GIS using this raster
# as the unit cost layer  and the health facility points as the 'start points'
# to create a raster layer of catchment areas which can then be vectorised

# Plot the resulting raster
plot(temp.raster)

tm_shape(temp.raster)+
  tm_raster(palette = "Greens", style = "fisher", n = 5)+
  tm_layout(legend.position = c("left","bottom"),
            frame = FALSE)



