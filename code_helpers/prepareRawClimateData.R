
#### Aim of prog: compute yearly average for the four variables tmax, tmean, tmin, and prec
## *Variables and units:
# tmax			mean daily maximum air temperature					K/10
# tmean			mean daily air temperature							K/10
# tmin			mean daily minimum air temperature					K/10
# prec			precipitation amount								kg.m-2
#
# tas			near-surface (usually, 2 meter) air temperature		K/10
# tasmin		idem, but minimum									K/10
# tasmax		idem, but maximum									K/10
# pr			precipitation flux									kg.m-2/100
#
## *Comments
# Many important informations can be found on the official website of the data (especially the units!):
#	https://chelsa-climate.org
#
# Here is the pdf containing the informations for the version 2 (version read: 24 May 2021)
#	https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf
#
# The three temperature variables are the daily maximum/average/minimum air temperatures at
#	2 metres since previous post-processing from ERA interim averaged over 1 month
#
# For the precipiations, "Amount" means mass per unit area while "Precipitation" in the earth's atmosphere
#	means precipitation of water in all phases.
#
## *Shapefile
#	The shapefile was created with the script create_usa-canada_shapefile.R

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)

#### Tool functions
## Function to get the dates of each raster (uses the naming pattern of Chelsa)
getDates = function(str, var, pattern = "varYearMonth", isThereUnderscore = TRUE)
{
	recognisedPatterns = c("varYearMonth", "varMonthYear")
	if (!(pattern %in% recognisedPatterns))
		stop(paste0("Pattern <", pattern, "> unknown"))

	if (pattern == "varYearMonth")
	{
		indexStart = stri_locate(str, regex = var)[, "end"] + ifelse(isThereUnderscore, 2, 1)
		getYear = stri_sub(str = str, from = indexStart, to = indexStart + 3)
		getMonth = stri_sub(str = str, from = indexStart + 3 + ifelse(isThereUnderscore, 2, 1),
			to = indexStart + 4 + ifelse(isThereUnderscore, 2, 1))
		monthName = month.abb[as.integer(getMonth)]
	}

	if (pattern == "varMonthYear")
	{
		indexStart = stri_locate(str, regex = var)[, "end"] + ifelse(isThereUnderscore, 2, 1)
		getMonth = stri_sub(str = str, from = indexStart, to = indexStart + 1)
		getYear = stri_sub(str = str, from = indexStart + 1 + ifelse(isThereUnderscore, 2, 1),
			to = indexStart + 4 + ifelse(isThereUnderscore, 2, 1))
		monthName = month.abb[as.integer(getMonth)]
	}

	return (data.table(fileName = str, year = getYear, month = getMonth, monthName = monthName))
}

## Compute min over a quarter (i.e., the three contiguous min value of the year for a given location). Example: driest quarter of the year
min_quarter = function(stack_rs)
{
	combinations = data.table(l1 = 1:10, l2 = 2:11, l3 = 3:12)
	nComb = combinations[, .N]
	temporary_stack = vector(mode = "list", length = nComb)
	names(temporary_stack) = paste0("l", 1:nComb)
	for (i in 1:nComb)
		temporary_stack[[paste0("l", i)]] = app(x = subset(stack_rs, c(i, i + 1, i + 2)), fun = "sum")
	rs = rast(temporary_stack)
	rs = min(rs)
	return(rs)
}

#### Define common variables
## Folder containing the data
folderRawData = "../raw_data/"

## Load map for cropping
usa_canada = vect(paste0(folderRawData, "usa-canada.shp"))

ls_variables = list.dirs(path = folderRawData, recursive = FALSE, full.names = FALSE)
selectedVariables = c("tas", "tasmin", "tasmax", "pr")

if (any(!(selectedVariables %in% ls_variables)))
{
	warning(paste0("The following variables are not available and are removed:\n- ",
		paste0(selectedVariables[!(selectedVariables %in% ls_variables)], collapse = "\n- ")))

	selectedVariables = selectedVariables[selectedVariables %in% ls_variables]

	if (length(selectedVariables) == 0)
		stop("No selected variable recognised")
}

#### Compute average
for (var in selectedVariables)
{
	pathToData = paste0(folderRawData, var, "/")
	ls_data = dir(path = pathToData, pattern = paste0("^CHELSA_", var, ".*.tif$"))
	ls_dates = getDates(str = ls_data, var = var, pattern = "varMonthYear")
	
	ls_dates[, dataPerYear := .N, by = year]
	if (any(ls_dates[, unique(dataPerYear)] != 12))
	{
		stop(paste0("Some months are missing for the variable <", var, "> for the following years:\n- ",
			paste0(ls_dates[dataPerYear != 12, unique(year)], collapse = "\n- ")))
	}

	if (ls_dates[, .N] != 12)
		stop("There is either more than one year, or not a single complete year")

	setorder(ls_dates, year, month)
	year = ls_dates[1, year]

	# Load climate for the year
	climate_rs = rast(paste0(pathToData, ls_dates[, fileName]))

	# Reproject
	usa_canada = project(usa_canada, climate_rs)

	# Crop to USA/Canada
	climate_rs = crop(x = climate_rs, y = usa_canada)

	# Compute 'min' quarter (such as driest, coolest)
	quarter_min_rs = min_quarter(climate_rs)

	# Compute annual min
	annual_min_rs = min(climate_rs)

	# Compute annual mean
	annual_mean_rs = mean(climate_rs)

	# Compute temperature in celsius (°C), check at the beginning of this script for the x/10:
	if (var %in% c("tmax", "tmean", "tmin", "tas", "tasmin", "tasmax"))
	{
		quarter_min_rs = app(x = quarter_min_rs, fun = function(x){x/10 - 273.15})
		annual_min_rs = app(x = annual_min_rs, fun = function(x){x/10 - 273.15})
		annual_mean_rs = app(x = annual_mean_rs, fun = function(x){x/10 - 273.15})
	}

	if (var == "pr")
	{
		quarter_min_rs = quarter_min_rs/100 # See comments for /100
		annual_min_rs = 12/100*annual_min_rs # Times 12 is to get per year (rather than per month), and see comments for /100
		annual_mean_rs = 12/100*annual_mean_rs
	}

	savingPath = paste0(folderRawData, var, "/")

	mask(x = quarter_min_rs, mask = usa_canada, inverse = FALSE, updatevalue = NA, touches = TRUE,
		filename = paste0(savingPath, year, "_quarter_min.tif"), filetype = "GTiff", overwrite = TRUE)
	mask(x = annual_min_rs, mask = usa_canada, inverse = FALSE, updatevalue = NA, touches = TRUE,
		filename = paste0(savingPath, year, "_annual_min.tif"), filetype = "GTiff", overwrite = TRUE)
	mask(x = annual_mean_rs, mask = usa_canada, inverse = FALSE, updatevalue = NA, touches = TRUE,
		filename = paste0(savingPath, year, "_annual_mean.tif"), filetype = "GTiff", overwrite = TRUE)

	print(paste(var, " done"))
}
