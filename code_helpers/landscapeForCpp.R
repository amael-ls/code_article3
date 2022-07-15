
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

#### Common variables
## Folders
loadPath = "../raw_data/"
outputPath = "../run/data/landscape_300x7/"

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
refugiaOption = FALSE
nbRefugia = 10

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

## Define crop extent number 2 (Old Tappan Borough, NJ 07675, USA):
# lonMin = -74
# lonMax = -73.9
# latMin = 41
# latMax = 41.2

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

#! --- Crash test zone, to have a smaller landscape
#? What follows is for the crop extent number 1 (Mont-Orford National Park)
# For a 300 x 11 landscape
crop_extent = ext(c(626881, 627106, 151133, 157253)) # order = xmin, xmax, ymin, ymax

#? What follows is for the crop extent number 2 (Old Tappan Borough, NJ 07675, USA)
# # For a 100 x 100 landscape
# crop_extent = extent(c(xxx, xxx, xxx, xxx))

# # For a 200 x 7 landscape
# crop_extent = extent(c(xxx, xxx, xxx, xxx))

# # For a 100 x 5 landscape
# crop_extent = extent(c(xxx, xxx, xxx, xxx))

# # For a 10 x 100 landscape
# crop_extent = extent(c(xxx, xxx, xxx, xxx))

# # For a 20 x 20 landscape
# crop_extent = extent(c(xxx, xxx, xxx, xxx))

# # For a 5 x 5 landscape
# crop_extent = extent(c(xxx, xxx, xxx, xxx))
#! --- End crash test zone, to have a smaller landscape

croppedClimate = crop(croppedClimate, crop_extent)

(nrows = nrow(croppedClimate))
(ncols = ncol(croppedClimate))

deltaX = xres(croppedClimate)
deltaY = yres(croppedClimate)

## Get centroid and in which cell it belongs
coords = crds(croppedClimate)
centroid = colMeans(coords)

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

# Populations are at the 30 bottom lines of the landscape #! i.e., 10% of the 300*11 landscape, starting from the bottom
vals[row %in% seq(max(row) - 29, max(row)), isPopulated := "true"]

# Populations are at the 100 top lines of the landscape
vals[row %in% seq(min(row), min(row) + 99), isPopulated := "true"]

# Populations are at the 180 top lines of the landscape
vals[row %in% seq(min(row), min(row) + 179), isPopulated := "true"]

# Populations are at the 100 bottom lines of the landscape
vals[row %in% seq(max(row) - 99, max(row)), isPopulated := "true"]

# Print number of colonised patches
print(paste0("Number of colonised patches = ", sum(vals$isPopulated == "true")))
print(paste0("Number of colonised rows = ", sum(vals$isPopulated == "true")/ncols))

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
vals[, writeCppClimate_DT(unlist(vals[rowNumber]), cppNames, patch_id, crs_str, deltaX, deltaY, path = outputPath, sep = " = ", rm = TRUE), by = rowNumber]

## Id files of populated patches (to create the Initial Condition in createIC folder)
# List occupied patches
patch_id = vals[isPopulated == "true", patch_id]

# Abies balsamea
maxRow_abba = 270
max_index_limitRow = (maxRow_abba - 1)*ncols + ncols - 1
ls_id_abba = patch_id[patch_id <= max_index_limitRow]
abba = data.table(patch_id = ls_id_abba, species = "Abies_balsamea")

# Acer saccharum
minRow_acsa = 271
min_index_limitRow = (minRow_acsa - 1)*ncols
ls_id_acsa = patch_id[patch_id >= min_index_limitRow]
acsa = data.table(patch_id = ls_id_acsa, species = "Acer_saccharum")

# Merge all the species in a data.table
patch_data = rbind(abba, acsa)
patch_data[, .N, by = "species"]

# Save patch data
saveRDS(patch_data, paste0(outputPath, "populatedPatches.rds"))

## Write data to create Matlab's data
saveRDS(vals, paste0(outputPath, "climate.rds"))

# #### Plot
# pdf("test.pdf")
# plot(croppedClimate, axes = FALSE, col = viridis::viridis(256))
# dev.off()

# #### Plot to check centroid
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
