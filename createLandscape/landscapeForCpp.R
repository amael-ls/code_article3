
#### Aim of prog: Create the landscape for C++ program
## Temperature:
#		- annual_mean_temperature (growth)
#		- min_temperature_of_coldest_month (mortality)
#
## Precipitation:
#		- annual_precipitation (growth)
#		- precipitation_of_driest_quarter (mortality)
#
## Remarks:
# The climatic data are stored in the article 1 folder on Beluga (Compute Canada)
#

#### Load packages
library(data.table)
library(stringi)
library(raster)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

#### Tool functions
## Write climate values for C++ program, with cppNames for the C++ names
writeCppClimate = function(climate, cppNames, id, crs, sep = " = ", rm = FALSE)
{
	outfileName = paste0("./climate_", id, ".txt")
	if (rm & file.exists(outfileName))
		file.remove(outfileName)

	sink(file = outfileName, append = TRUE)
	for (i in 1:length(cppNames))
		cat(paste0(cppNames[i], sep, climate[, cppNames[i], with = FALSE][[1]]), sep = "\n")

	cat(paste0("plotArea", sep, 400), sep = "\n") # To modify!!!
	cat(paste0("proj4string", sep, crs), sep = "\n") # To modify!!!

	sink(file = NULL)
}

## Tu use multi threading from data.table
writeCppClimate_DT = function(climateNamedVector, cppNames, id, crs, sep = " = ", rm = FALSE)
{
	outfileName = paste0("./climate_", id, ".txt")
	if (rm & file.exists(outfileName))
		file.remove(outfileName)

	sink(file = outfileName, append = TRUE)
	for (i in 1:length(cppNames))
		cat(paste0(cppNames[i], sep, climateNamedVector[cppNames[i]]), sep = "\n")

	cat(paste0("plotArea", sep, 400), sep = "\n") # To modify!!!
	cat(paste0("proj4string", sep, crs), sep = "\n") # To modify!!!

	sink(file = NULL)
}

#### Common variables
## Folders
loadPath = "~/projects/def-dgravel/amael/article1/progToSendToReview/"

## C++ names
cppNames = c("annual_mean_temperature", "min_temperature_of_coldest_month",
	"annual_precipitation", "precipitation_of_driest_quarter", "longitude", "latitude")

## Load data
# Climate
climate = readRDS(paste0(loadPath, "createData/clim_2010.rds"))

# Projection
crs = crs(raster(paste0(loadPath, "clim60sec/clim_2010.grd")), asText = TRUE)

## Rename columns climate data
climVar = readRDS(paste0(loadPath, "createData/climaticVariables.rds"))
setnames(climate, old = names(climate), new = c("latitude", "longitude", climVar))

## Set (x, y) coordinates in NYC (A tree grows in Brooklyn...)
x = 2150631.67
y = -131511.69

climate[, dist := sqrt((latitude - y)^2 + (longitude - x)^2)]
# Apparently East river state park, fortunately there is uniqueness of the minimum!
clim = climate[dist == min(dist), ]

#### Write file for Cpp
writeCppClimate(clim, cppNames, id = 1, crs, sep = " = ", rm = TRUE)

#### Get climate within a disc
## Transform dist in m to km
climate[, dist := dist/1000]

clim = climate[dist < 5]

for (id in 1:clim[, .N])
	writeCppClimate(clim[id], cppNames, id, crs, sep = " = ", rm = TRUE)

#### Get climate from raster
## Load raster
climate_rs = stack(paste0(loadPath, "clim60sec/clim_2010.grd"))

# ## Transform raster if not done yet
# if (!file.exists("./reprojRast/clim_2010.grd"))
# {
# 	#Load raster
# 	climate_rs = raster(paste0(loadPath, "clim60sec/clim_2010.grd"))

# 	# Transform
# 	climate_rs = projectRaster(climate_rs, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# 	dir.create("./reprojRast")
# 	writeRaster(climate_rs, filename = "reprojRast/clim_2010.grd", bandorder = "BIL")
# }

## Define crop extent
lonMin = -74
lonMax = -73.9
latMin = 41
latMax = 41.1

crop_extent = extent(c(lonMin, lonMax, latMin, latMax))
crop_extent = projectExtent(raster(crop_extent, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
	crs(climate_rs))

croppedClimate = crop(climate_rs[[c("bio60_01", "bio60_06", "bio60_12", "bio60_17")]], crop_extent)
names(croppedClimate) = c("annual_mean_temperature", "min_temperature_of_coldest_month",
	"annual_precipitation", "precipitation_of_driest_quarter")

vals = as.data.table(rasterToPoints(croppedClimate))
setnames(vals, old = c("x", "y"), new = c("latitude", "longitude"))

for (id in 1:vals[, .N])
	writeCppClimate(vals[id], cppNames, id, crs, sep = " = ", rm = TRUE)

vals[, id := 1:.N]
vals[, writeCppClimate_DT(unlist(vals[id]), cppNames, id, crs, sep = " = ", rm = TRUE), by = id]