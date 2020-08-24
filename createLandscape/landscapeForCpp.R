
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
cppNames = c("annual_mean_temperature", "min_temperature_of_coldest_month", "annual_precipitation",
	"precipitation_of_driest_quarter", "isPopulated", "longitude", "latitude", "patch_id")

#### Get climate from raster
## Load raster
climate_rs = stack(paste0(loadPath, "clim60sec/clim_2010.grd"))

# Projection
crs = crs(climate_rs, asText = TRUE)

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

## Change resolution to a 20 x 20 m
# downScale_factors = floor(res(croppedClimate)/20)
# croppedClimate_downscale = disaggregate(x = croppedClimate, fact = downScale_factors)

# nrow(croppedClimate)
# ncol(croppedClimate)

## Get centroid and in which cell it belongs
centroid = colMeans(coordinates(croppedClimate))

## Coerce to data table
vals = as.data.table(rasterToPoints(croppedClimate))
setnames(vals, old = c("x", "y"), new = c("longitude", "latitude"))

## Add id and isPopulated
vals[, patch_id := 0:(.N - 1)] # C++ starts at 0, not 1
vals[, rowNumber := 1:.N]
vals[, isPopulated := "false"]

## Compute distance to centroid
vals[, dist := sqrt((longitude - centroid["x"])^2 + (latitude - centroid["y"])^2)]

## Set patches that are populated (currently, only the closest to the centroid)
vals[dist == min(dist), isPopulated := "true"]

#### Save files
## Environment files for C++ prog
vals[, writeCppClimate_DT(unlist(vals[rowNumber]), cppNames, patch_id, crs, sep = " = ", rm = TRUE), by = rowNumber]

## Id files of populated patches (to create the Initial Condition in createIC folder)
saveRDS(vals[ifelse(isPopulated == "true", TRUE, FALSE), patch_id], "./populatedPatches.rds")

#### Plot to check centroid
# pdf("test.pdf", height = 8, width = 8)
# croppedClimate = setValues(x = croppedClimate, values = sample(1:ncell(croppedClimate), ncell(croppedClimate), replace = TRUE), layer = 1)
# plot(croppedClimate[["annual_mean_temperature"]])
# text(x = coordinates(croppedClimate)[, "x"],
# 	y = coordinates(croppedClimate)[, "y"], label = 1:ncell(croppedClimate), pos = 3)
# text(x = coordinates(croppedClimate)[, "x"],
# 	y = coordinates(croppedClimate)[, "y"], label = paste0("(", values(croppedClimate[["annual_mean_temperature"]]), ")"), pos = 1)
# points(x = centroid["x"], y = centroid["y"], pch = 15, col = "black")
# points(x = vals[dist == min(dist), longitude], y = vals[dist == min(dist), latitude], pch = 19, col = "blue")
# dev.off()
