
#### Aim of prog: Create the initial condition for the C++ prog
## Initial condition from real data
#
## Initial condition from a continuous function (here a lognormal)
#

#### Load packages
library(data.table)
library(stringi)
library(raster)
library(rgdal)
library(sf)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

######## PART I: CREATE DATA REAL DATA
#### Load data and keep columns of interest
loadPath = "~/projects/def-dgravel/amael/article1/progToSendToReview/"
trees = readRDS(paste0(loadPath, "createData/tree_sStar.rds"))
trees = trees[, .(plot_id, tree_id, year_measured, species_id, dbh, latitude, longitude, plot_size, s_star)]

## Compute number of individuals per plot/year/species
trees[, nbIndiv := .N, by = c("plot_id", "year_measured", "species_id")]

## Compute number of individuals per plot/year
trees[, nbIndivPlotYear := .N, by = c("plot_id", "year_measured")]

#### Create initial condition
## Keep the maximum individuals of one sp per year per plot
treesIC = trees[nbIndiv == max(nbIndiv)]

## Check whether there is only one species (should be the case)
print(treesIC[, unique(species_id)])

## Compute densities per dbh (we are already in a plot/species/year group)
treesIC[, nbIndiv := .N, by = dbh]
treesIC[, density := nbIndiv/plot_size]

coords = unique(treesIC[, .(longitude, latitude)])
year = treesIC[, unique(year_measured)]

treesIC = unique(treesIC[, .(density, dbh)])

## Sort individuals by decreasing order (although Cpp prog sorts if necessary)
setorder(treesIC, -dbh) # Sort by decreasing dbh

#### Get climate (averaged over 5 years) associated to the year, plot_id
## Read data
climate = readRDS(paste0(loadPath, "createData/averagedClim_5yearsAllVar.rds"))

## Select the year
climate = climate[year_measured == year]

## Select coordinates
climate = climate[coords, nomatch = 0, on = names(coords)]

## Select the four variables of interest
climate = climate[, .(annual_mean_temperature, min_temperature_of_coldest_month,
	annual_precipitation, precipitation_of_driest_quarter)]

print(climate) # To write them manually if required

#### Write the initial condition
write.table(treesIC, file = "./ic_1.txt", append = FALSE, quote = FALSE, sep = " ",
	eol = "\n", na = "NA", dec = ".", row.names = FALSE,
	col.names = TRUE)

#### Get a landscape of climate
## Load raster
climate_rs = stack(paste0(loadPath, "clim60sec/clim_2010.grd"))
crs_ref = crs(climate_rs)

## Transform coords
coords_sp = SpatialPoints(coords, proj4str = CRS("+init=epsg:4326"))
coords_sp = spTransform(coords_sp, CRSobj = crs_ref)
coords_sp = coordinates(coords_sp)

## Define the crop extent (square of 40km around coords_sp)
bottomLeft_lon = coords_sp[1, "longitude"] - 20000
bottomLeft_lat = coords_sp[1, "latitude"] - 20000

topLeft_lon = coords_sp[1, "longitude"] - 20000
topLeft_lat = coords_sp[1, "latitude"] + 20000

topRight_lon = coords_sp[1, "longitude"] + 20000
topRight_lat = coords_sp[1, "latitude"] + 20000

bottomRight_lon = coords_sp[1, "longitude"] + 20000
bottomRight_lat = coords_sp[1, "latitude"] - 20000

coords_mx = matrix(data = c(
		bottomLeft_lon, bottomLeft_lat,
		topLeft_lon, topLeft_lat,
		topRight_lon, topRight_lat,
		bottomRight_lon, bottomRight_lat,
		bottomLeft_lon, bottomLeft_lat), ncol = 2, byrow = TRUE)

crop_extent = st_polygon(list(coords_mx))
crop_extent = as(st_sfc(crop_extent), "Spatial")

## Crop the raster
croppedClimate_rs = crop(climate_rs, crop_extent)

## Keep only the 4 variables of interest (BIO1, BIO6, BIO12, BIO17)
# BIO1 = Annual Mean Temperature
# BIO6 = Min Temperature of Coldest Month
# BIO12 = Annual Precipitation
# BIO17 = Precipitation of Driest Quarter

croppedClimate_rs = subset(croppedClimate_rs, subset = paste0("bio60_", c("01", "06", "12", "17")))

## Get the spatial dimension
m_nRow = nrow(croppedClimate_rs)
m_nCol = ncol(croppedClimate_rs)

## Coerce to data table
croppedClimate_rs = rasterToPoints(croppedClimate_rs)
croppedClimate_rs = as.data.table(croppedClimate_rs)

setnames(croppedClimate_rs, old = names(croppedClimate_rs),
	new = c("longitude", "latitude",
	"annual_mean_temperature", "min_temperature_of_coldest_month",
	"annual_precipitation", "precipitation_of_driest_quarter"))

#### Write the metadata for C++
## Where to save the data
metadata_fileName = "clim_metadata.txt"
savingPath = "./landscapeData"
sep = " = "

## Metadata
savingPath_meta = savingPath
if (stri_detect(savingPath, regex = "^./"))
	savingPath_meta = stri_sub(savingPath, from = stri_locate(savingPath, regex = "./")[2] + 1)

m_path = paste0("../createIC/", savingPath_meta, "/")
m_delimiter= "\" = \""

if (!dir.exists(savingPath))
	dir.create(savingPath)

if (file.exists(paste0(savingPath, metadata_fileName)))
	file.remove(paste0(savingPath, metadata_fileName))

sink(file = paste0("./", metadata_fileName), append = TRUE)

cat(paste0("nRow", sep, m_nRow, sep = "\n"))
cat(paste0("nCol", sep, m_nCol, sep = "\n"))
cat(paste0("path", sep, m_path), sep = "\n")
cat(paste0("delimiter", sep, m_delimiter), sep = "\n")

sink(file = NULL)

#### Write the files
## Remarks:
#	1/ rasterToPoints goes line by line, i.e., for a raster (n, p), n rows and p cols
#		the first p lines of the data table is the first row of the raster and so on.
for (i in 1:nrow(croppedClimate_rs))
{
	sink(file = paste0(savingPath, metadata_fileName), append = TRUE)

	cat(paste0("nRow", sep, m_nRow, sep = "\n"))
	cat(paste0("nCol", sep, m_nCol, sep = "\n"))
	cat(paste0("path", sep, m_path), sep = "\n")
	cat(paste0("delimiter", sep, m_delimiter), sep = "\n")

	sink(file = NULL)
}


#### CRASH TEST ZONE
# jpeg("test.jpg")
# plot(croppedClimate_rs[[1]])
# plot(crop_extent, lwd = 2, col = "blue", add = TRUE)
# dev.off()

# jpeg("test2.jpg")
# plot(croppedClimate_rs[[1]])
# plot(crop_extent, lwd = 2, border = "blue", add = TRUE)
# points(coords_sp, col = "red", pch = 19, cex = 2)
# dev.off()

#### END CRASH TEST ZONE

######## PART II: CREATE DATA FROM FUNCTION
#### Create data
# To do if required
