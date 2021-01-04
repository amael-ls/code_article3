
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
## Write climate values for C++ program, with cppNames for the C++ names using multi threading from data.table
writeCppClimate_DT = function(climateNamedVector, cppNames, id, crs, deltaX, deltaY, path = "./", sep = " = ", rm = FALSE)
{
	if (!(stri_detect(str = path, regex = "/")) | (stri_locate_last(str = path, regex = "/")[, "end"] != stri_length(path)))
		path = paste0(path, "/")
	
	outfileName = paste0(path, "climate_", id, ".txt")
	if (rm & file.exists(outfileName))
		file.remove(outfileName)

	ofstream = file(outfileName)
	line = ""
	for (i in 1:length(cppNames))
		line = paste0(line, cppNames[i], sep, climateNamedVector[cppNames[i]], "\n")
	
	line = paste0(line, "plotArea", sep, deltaX*deltaY, sep = "\n")
	line = paste0(line, "proj4string", sep, crs)

	writeLines(line, ofstream)
	close(ofstream)
}

#### Common variables
## Folders
loadPath = "~/projects/def-dgravel/amael/article1/progToSendToReview/"
outputPath = "./climate_200x7_refugia/"

if (!dir.exists(outputPath))
	dir.create(outputPath)

if (length(list.files(outputPath, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0)
	unlink(paste0(outputPath, "*"))

## C++ names
cppNames = c("annual_mean_temperature", "min_temperature_of_coldest_month", "annual_precipitation",
	"precipitation_of_driest_quarter", "isPopulated", "longitude", "latitude", "patch_id", "row", "col")

## Options
# Averaged conditions on landscape
homogeneousClimate = TRUE

# Refugia beyond southern distribution
refugiaOption = TRUE
nbRefugia = 10

#### Get climate from raster
## Load raster
climate_rs = stack(paste0(loadPath, "clim60sec/clim_2010.grd"))

# Projection
crs = crs(climate_rs, asText = TRUE)

## Define crop extent
lonMin = -74
lonMax = -73.9
latMin = 41
latMax = 41.2
# lonMin = -74
# lonMax = -73.99
# latMin = 41.09
# latMax = 41.1

crop_extent = extent(c(lonMin, lonMax, latMin, latMax))
crop_extent = projectExtent(raster(crop_extent, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
	crs(climate_rs))

croppedClimate = crop(climate_rs[[c("bio60_01", "bio60_06", "bio60_12", "bio60_17")]], crop_extent)
names(croppedClimate) = c("annual_mean_temperature", "min_temperature_of_coldest_month",
	"annual_precipitation", "precipitation_of_driest_quarter")

nrows = nrow(croppedClimate)
ncols = ncol(croppedClimate)

deltaX = xres(croppedClimate) # 1060
deltaY = yres(croppedClimate) # 1830

## Change resolution to a 20 x 20 m
downScale_factors = floor(res(croppedClimate)/20)
croppedClimate = disaggregate(x = croppedClimate, fact = downScale_factors, method = "bilinear")

#! --- Crash test zone, to have a smaller landscape
# For a 100 x 100 landscape
crop_extent = extent(c(2135078, 2137078, -93158, -91146)) # order = xmin, xmax, ymin, ymax

# For a 200 x 7 landscape
crop_extent = extent(c(2135078, 2135218, -93158, -89136))

# For a 100 x 5 landscape
crop_extent = extent(c(2135078, 2135178, -93158, -91146))

# For a 10 x 100 landscape
crop_extent = extent(c(2135078, 2137078, -93158, -92957))

# For a 20 x 20 landscape
crop_extent = extent(c(2135078, 2135478, -93158, -92758))

# For a 5 x 5 landscape
crop_extent = extent(c(2135078, 2135178, -93158, -93057))
#! --- End crash test zone, to have a smaller landscape

croppedClimate = crop(croppedClimate, crop_extent)

nrows = nrow(croppedClimate)
ncols = ncol(croppedClimate)

deltaX = xres(croppedClimate) # 20
deltaY = yres(croppedClimate) # 20.10989

## Get centroid and in which cell it belongs
centroid = colMeans(coordinates(croppedClimate))

## Coerce to data table
vals = as.data.table(rasterToPoints(croppedClimate))
setnames(vals, old = c("x", "y"), new = c("longitude", "latitude"))

## If homogeneous climate, compute average
if (homogeneousClimate)
{
	cols = c("annual_mean_temperature", "min_temperature_of_coldest_month",
		"annual_precipitation", "precipitation_of_driest_quarter")
	vals[ , (cols) := lapply(.SD, "mean"), .SDcols = cols]
}

## Add id, row, col, and isPopulated
vals[, patch_id := 0:(.N - 1)] # C++ starts at 0, not 1
vals[, col := patch_id %% ncols]
vals[, row := (patch_id - col)/ncols]
vals[, rowNumber := 1:.N]
vals[, isPopulated := "false"]

#! Select one or more option below
## Compute distance to centroid
vals[, dist := sqrt((longitude - centroid["x"])^2 + (latitude - centroid["y"])^2)]

## Set patches that are populated
# Populations are the closest to the centroid
vals[dist == min(dist), isPopulated := "true"]

# Populations are at the bottom of the landscape
vals[row == max(row), isPopulated := "true"]

# Populations are at the 10 bottom lines of the landscape
vals[row %in% seq(max(row) - 9, max(row)), isPopulated := "true"]

# Populations are at the 20 bottom lines of the landscape
vals[row %in% seq(max(row) - 19, max(row)), isPopulated := "true"]

# Add refugia
if (refugiaOption)
{

	#! --- Crash test zone, to have refugia
	# Randomly assignated refugia
	set.seed(1969-08-18) # Woodstock seed
	if (vals[isPopulated != "true", .N] < nbRefugia)
	{
		nbRefugia = vals[!isPopulated, .N]
		print("*** The number of refugia was above the number of free patches. All the landscape is colonised! ***")
	}
	
	refugia = sample(vals[isPopulated != "true", patch_id], nbRefugia, replace = FALSE)

	# Refugia assigned by user
	#* This is for the 200x7 landscape: (row 65, 66 and cols 0, 1); (row 12-14 and cols 4, 5)
	refugia = c(64*7, 64*7 + 1, 65*7, 65*7 + 1,
		11*7 + 3, 12*7 + 4, 12*7 + 3, 13*7 + 4, 13*7 + 3, 11*7 + 4)
	#! --- End crash test zone, to have refugia

	# Assigned refugia
	vals[patch_id %in% refugia, isPopulated := "true"]
}

#### Save files
## Environment files for C++ prog
vals[, writeCppClimate_DT(unlist(vals[rowNumber]), cppNames, patch_id, crs, deltaX, deltaY, path = outputPath, sep = " = ", rm = TRUE), by = rowNumber]

## Id files of populated patches (to create the Initial Condition in createIC folder)
saveRDS(vals[ifelse(isPopulated == "true", TRUE, FALSE), patch_id], paste0(outputPath, "populatedPatches.rds"))

## Write data to create Matlab's data
saveRDS(vals, paste0(outputPath, "climate.rds"))

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
