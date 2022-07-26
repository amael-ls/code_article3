
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
#	1. The raw climatic data are stored in ../raw_data
#	2. I use the projection EPSG:3799 simply because it is easier to count with meters than degrees. It is quite precise for Quebec
#		For another region, it would worth taking something else
#

#### Load packages
library(data.table)
library(stringi)
library(terra)

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

## Write file landscape.txt
writeLandscape_txt = function(values, outfileName = "../run/landscape.txt", sep = "=", rm_dots = TRUE, rm_dots_vars = "path", rm = TRUE)
{
	if (!is.list(values))
		stop("Values must be a list")

	varNames = names(values)
	
	if (rm & file.exists(outfileName))
		file.remove(outfileName)

	if (rm_dots)
		values[[rm_dots_vars]] = stri_replace(str = values[[rm_dots_vars]], replacement = "./", regex = "^../")

	ofstream = file(outfileName)
	line = ""
	for (current_var in varNames)
		line = paste0(line, current_var, sep, values[[current_var]], "\n")

	writeLines(line, ofstream)
	close(ofstream)
}

#### Common variables
## Folders
loadPath = "../raw_data/"
outputPath = "../run/data/landscape_5x5_acsa/"

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

#### Get climate from raster
## Load raster
ls_variables = list.dirs(path = loadPath, recursive = FALSE, full.names = FALSE)
ls_rast = paste0(loadPath, ls_variables, "/2010_",
	rep(c("quarter_min.tif", "annual_min.tif", "annual_mean.tif"), each = length(ls_variables)))
climate_rs = rast(ls_rast)

names(climate_rs) = paste0(ls_variables, "_", rep(c("quarter_min", "annual_min", "annual_mean"), each = length(ls_variables)))

## Define crop extent number 1 (Mont-Orford National Park):
lonMin = -72.207962
lonMax = -72.107962
latMin = 45.345086
latMax = 45.545086

crop_extent = ext(c(lonMin, lonMax, latMin, latMax))

croppedClimate = crop(climate_rs[[c("tas_annual_mean", "tas_annual_min", "pr_annual_mean", "pr_quarter_min")]], crop_extent)

names(croppedClimate) = c("annual_mean_temperature", "min_temperature_of_coldest_month",
	"annual_precipitation", "precipitation_of_driest_quarter")

# Projection
croppedClimate = project(x = croppedClimate, "epsg:3799")

crs_str = crs(croppedClimate, describe = TRUE)[c("name", "authority", "code")]
crs_str = paste0(crs_str["name"], ", ", crs_str["authority"], ":", crs_str["code"])

nrows = nrow(croppedClimate)
ncols = ncol(croppedClimate)

deltaX = xres(croppedClimate)
deltaY = yres(croppedClimate)

## Fill the gaps created by the projection (no idea why...)
croppedClimate = focal(croppedClimate, w = matrix(data = c(1, 1, 1, 1, 0, 1, 1, 1, 1), nrow = 3),
	fun = mean, na.policy = "only", na.rm = TRUE)

## Change resolution to a 20 x 20 m
downScale_factors = floor(res(croppedClimate)/20)
croppedClimate = disagg(x = croppedClimate, fact = downScale_factors, method = "bilinear")

## Crop to a 5 x 5 landscape
crop_extent = ext(c(626881, 626983, 151133, 151235))
croppedClimate = crop(croppedClimate, crop_extent)

(nrows = nrow(croppedClimate))
(ncols = ncol(croppedClimate))

deltaX = xres(croppedClimate)
deltaY = yres(croppedClimate)

coords = crds(croppedClimate)

## Coerce to data table
vals = as.data.table(as.points(x = croppedClimate, values = TRUE, na.rm = TRUE))
if (vals[, .N] != nrows*ncols)
	warning("You probably have NAs in your raster. Check it, as this brings problems later on...")

vals[, c("longitude", "latitude") := split(coords, rep(1:ncol(coords), each = nrow(coords)))]

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

## Set patches that are populated
# Populations are at the bottom of the landscape
vals[row == max(row), isPopulated := "true"]

# Print number of colonised patches
print(paste0("Number of colonised patches = ", sum(vals$isPopulated == "true")))
print(paste0("Number of colonised rows = ", sum(vals$isPopulated == "true")/ncols))

#### Save files
## Environment files for C++ prog
vals[, writeCppClimate_DT(unlist(vals[rowNumber]), cppNames, patch_id, crs_str, deltaX, deltaY,
	path = outputPath, sep = " = ", rm = TRUE), by = rowNumber]

writeLandscape_txt(values = list(
	nRow = nrow(croppedClimate),
	nCol = ncol(croppedClimate),
	path = outputPath,
	delimiter = " = ",
	filenamePattern = "climate_",
	deltaLon = xres(croppedClimate),
	deltaLat = yres(croppedClimate),
	distance = "euclidean"))

## Id files of populated patches (to create the Initial Condition in createIC folder)
# List occupied patches
patch_id = vals[isPopulated == "true", patch_id]

# Acer saccharum
minRow_acsa = 5
min_index_limitRow = (minRow_acsa - 1)*ncols
ls_id_acsa = patch_id[patch_id >= min_index_limitRow]
patch_data = data.table(patch_id = ls_id_acsa, species = "Acer_saccharum")

# Save patch data
saveRDS(patch_data, paste0(outputPath, "populatedPatches.rds"))
