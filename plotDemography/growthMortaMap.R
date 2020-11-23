
#### Aim of prog: Plot maps of growth and mortality, for s* = 10m

#### Aim of prog: Make maps from matlab results
library(data.table)
library(stringi)
library(raster)
library(sf)

options(max.print = 500)
rm(list = ls())

#### Tool function
## Get the neighbours values of a cell from a raster
getNeighboursVal = function(val_raster, id, ncell_r, ncol_r)
{
	# Get indices of the neighbours
	ind = c(id - ncol_r - 1, id - ncol_r, id - ncol_r + 1,
		id - 1, id + 1,
		id + ncol_r - 1, id + ncol_r, id + ncol_r + 1)

	# Remove indices outside of raster (happen when working on the edge)
	ind = ind[(ind > 0) & (ind < ncell_r)]
	return (val_raster[ind])
}

## Fill the gaps in raster, using neighbourhood. r = raster, p = polygon
fillGap = function(r, p)
{
	vals = getValues(r)
	vals_orig = vals # Needs it
	ind_na_orig = sort(which(is.na(vals)))

	# Get values in polygon p
	ind_p = cellFromPolygon(r, p) # list of indices for each polygon (if p is a multipolygon)
	ind_p = sort(unlist(ind_p))
	ind_p = unique(ind_p)

	vals = vals[ind_p]

	# Update ind_na_orig, and create a new ind_na within polygon
	ind_na_orig = ind_p[ind_p %in% ind_na_orig]
	ind_na = sort(which(is.na(vals)))
	ind_na_orig_withinPolygon = ind_na

	ncell_r = ncell(r)
	ncol_r = ncol(r)

	unstable = TRUE
	count = 0

	ind_dt = data.table(ind_val = ind_na, ind_orig = ind_na_orig)

	while (length(ind_na) != 0 & unstable)
	{
		prev_stage = vals
		print(paste0("Number of NAs: ", length(ind_na)))
		for (i in ind_na)
		{
			neighVals = getNeighboursVal(vals_orig, ind_dt[ind_val == i, ind_orig], ncell_r, ncol_r)
			if (sum(is.na(neighVals)) == 8)
				next;
			vals[i] = mean(neighVals, na.rm = TRUE)
		}

		vals_orig[ind_dt[ind_val %in% ind_na, ind_orig]] = vals[ind_na]

		ind_na = sort(which(is.na(vals)))

		ind_dt = ind_dt[ind_val %in% ind_na]

		if(identical(prev_stage, vals))
		{
			if (length (ind_na) == 1)
			{
				unstable = FALSE
			} else {
				count = count + 1
				ind_na = sample(ind_na)
				if (count == 10)
					unstable = FALSE
			}
		}
	}

	r[ind_na_orig] = vals[ind_na_orig_withinPolygon]
	return (r)
}

#### Common variables
## Path acer rubrum
acerub_g = "../growth/array_4/" 
acerub_m = "../mortality/array_4/"

height_canopy = "10m"

## Shapefiles Northern America
canada = readRDS("~/database/shapefiles/North_America_Shapefile/canada_union.rds")
usa = readRDS("~/database/shapefiles/North_America_Shapefile/usa_union.rds")

erie = st_read("~/database/shapefiles/North_America_Shapefile/lakes/erie")
huron = st_read("~/database/shapefiles/North_America_Shapefile/lakes/huron")
michigan = st_read("~/database/shapefiles/North_America_Shapefile/lakes/michigan")
ontario = st_read("~/database/shapefiles/North_America_Shapefile/lakes/ontario")
stClair = st_read("~/database/shapefiles/North_America_Shapefile/lakes/stClair")
superior = st_read("~/database/shapefiles/North_America_Shapefile/lakes/superior")

canada = st_transform(canada,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
usa = st_transform(usa,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
erie = st_transform(erie,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
huron = st_transform(huron,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
michigan = st_transform(michigan,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
ontario = st_transform(ontario,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
stClair = st_transform(stClair,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
superior = st_transform(superior,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")

canada = st_geometry(canada)
usa = st_geometry(usa)
erie = st_geometry(erie)
huron = st_geometry(huron)
michigan = st_geometry(michigan)
ontario = st_geometry(ontario)
stClair = st_geometry(stClair)
superior = st_geometry(superior)

## Load bounding box determined by the data. Will be used to crop the results
data_bbox = readRDS("../createData/bbox_treeData.rds")
data_extent = extent(matrix(data_bbox, nrow = 2, ncol = 2))

## Get species' name
species = list.files(path = acerub_g, pattern = "^[0-9]{,4}.*.txt")
species = stri_sub(str = species, to = stri_locate(str = species, regex = ".txt")[1] - 1)

## Load growth and mortality data (cf matlab program averageGrowthMortality.m)
averagedGrowth = fread(paste0("./results/", species, "/averageGrowth.csv"))[, V1]
averagedMortality = fread(paste0("./results/", species, "/averageMortality.csv"))[, V1]

## Load climate and merge with averagedG/M
clim_2010 = readRDS(paste0("../R0/Matlab_data/", species, "/", species, "_clim2010.rds"))
clim_2010 = clim_2010[, .(longitude, latitude)]
clim_2010[, c("G_ave", "M_ave") := .(averagedGrowth, averagedMortality)]

## Load the Little's 1971 shapefile, crop to data bbox
ls_little = list.files(path = paste0("../little1971/", species), pattern = ".shp$")
little1971 = st_read(paste0("../little1971/", species, "/", ls_little))
little1971 = st_transform(little1971,
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
little1971 = st_crop(little1971, data_bbox)
little1971 = st_geometry(little1971)
little1971_sp = as(little1971, "Spatial")

## Set raster of growth and mortality
growth_rs = rasterFromXYZ(clim_2010[, .(longitude, latitude, G_ave)],
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
growth_rs = crop(growth_rs, data_extent)

mortality_rs = rasterFromXYZ(clim_2010[, .(longitude, latitude, M_ave)],
	crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
mortality_rs = crop(mortality_rs, data_extent)

## Fill the NA within little1971
growth_rs = fillGap(growth_rs, little1971_sp)
mortality_rs = fillGap(mortality_rs, little1971_sp)

## Bounding box of plot
lonMin = unname(data_bbox["xmin"]) - 100
lonMax = unname(data_bbox["xmax"]) + 30000 # To get Nova Scotia completely
latMin = unname(data_bbox["ymin"]) - 60000 # To get Florida completely
latMax = unname(data_bbox["ymax"]) + 100

#### Plot maps
## Growth
tiff(paste0("./", species, "_averageG_h=", height_canopy, "_scaled.tiff"), width = 1000, height = 1000)
# jpeg(paste0("./", species, "_averageG_h=", height_canopy, "_scaled.jpg"), width = 1000, height = 1000, quality = 100)
plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
	xlab = "", ylab = "", bg = "transparent")
plot(growth_rs[["G_ave"]], cex = 0.5, pch = 19, add = TRUE, legend = FALSE, # , "#4269E2"
	col = c("#000000", "#2A3344", "#1122AA", "#2058DC", "#5A82EA", "#FDECBE", "#FDDB5B", "#FAB935", "#FD9859", "#B6343A"))
plot(little1971, col = NA, border = "#000000", add = TRUE, lwd = 2) # "#DD925C"
plot(canada, col = NA, add = TRUE, lwd = 2)
plot(usa, col = NA, add = TRUE, lwd = 2)
plot(erie, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(huron, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(michigan, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(ontario, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(stClair, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(superior, col = "#FFFFFF", border = "#000000", add = TRUE)
dev.off()

## Mortality
tiff(paste0("./", species, "_averageM_h=", height_canopy, "_scaled.tiff"), width = 1000, height = 1000)
# jpeg(paste0("./", species, "_averageM_h=", height_canopy, "_scaled.jpg"), width = 1000, height = 1000, quality = 100)
plot(0, pch = "", xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), axes = FALSE,
	xlab = "", ylab = "", bg = "transparent")
plot(mortality_rs[["M_ave"]], cex = 0.5, pch = 19, add = TRUE, legend = FALSE, # , "#4269E2"
	col = c("#000000", "#2A3344", "#1122AA", "#2058DC", "#5A82EA", "#FDECBE", "#FDDB5B", "#FAB935", "#FD9859", "#B6343A"))
plot(little1971, col = NA, border = "#000000", add = TRUE, lwd = 2) # "#DD925C"
plot(canada, col = NA, add = TRUE, lwd = 2)
plot(usa, col = NA, add = TRUE, lwd = 2)
plot(erie, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(huron, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(michigan, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(ontario, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(stClair, col = "#FFFFFF", border = "#000000", add = TRUE)
plot(superior, col = "#FFFFFF", border = "#000000", add = TRUE)
dev.off()

print(paste0("species: ", species, " done"))

