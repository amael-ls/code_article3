
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

# ####
# ## Load raster
# climate_rs = raster(paste0(loadPath, "clim60sec/clim_2010.grd"))

# ## Transform
# climate_rs = projectRaster(climate_rs, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# dir.create("./reprojRast")
# writeRaster(climate_rs, filename = "reprojRast/clim_2010.grd", bandorder = "BIL")


# ## Define crop extent
# bottomLeft
# topLeft
# topRight
# bottomRight

# crop_extent = SpatialPoints()
