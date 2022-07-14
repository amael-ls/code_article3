
#### Aim of prog: download Chelsa climate data
## Explanations
# 1. This script downloads climate data for the world.
# 2. The downloaded data are tif.
# 4. Documentation is in English
# 5. To cite the data, check in the doc:
# CHELSA data should be cited as:
#	- General citation: Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E.,
#		Linder, H.P. & Kessler, M. (2017) Climatologies at high resolution for the earth's land surface areas. Scientific Data 4, 170122.
#	- Data citations (Version 1.0): Karger, Dirk Nikolaus; Dabaghchian, Babek; Lange, Stefan; Thuiller, Wilfried; Zimmermann, Niklaus E.;
#		Graham, Catherine H. (2020). High resolution climate data for Europe. EnviDat. doi:10.16904/envidat.150.
# 6. Naming convention for files: CHELSA_<variable>_<month>_<Version>_land.tif
# 7. For European data, check https://chelsa-climate.org/high-resolution-climate-data-for-europe/
#
# Variables and units:
#	pr		= precipitation [mm/month]
#	tas		= monthly mean of daily mean temperature [°C*10]
#	tasmax	= monthly mean of daily maximum temperature [°C*10]
#	tasmin	= monthly mean of daily minimum temperature [°C*10]

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500, timeout = 200)

##################################################################
############            DATA FOR THE WORLD            ############
##################################################################

#### Tool function
## Get year and month from iteration. It is assumed that the years are complete (i.e., 12 months)
getYearMonth = function(id, years)
{
	month = ifelse(id %% 12 == 0, 12, id %% 12)
	year_ind = ifelse(id %% 12 == 0, id %/% 12, id %/% 12 + 1)
	return(list(year = years[year_ind], month = month))
}

#### Define common variables
## Folder where data will be downloaded
folderDownloadData = "../raw_data/"

if (!dir.exists(folderDownloadData))
	dir.create(folderDownloadData)

## Link to download the data
downloadLink = "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/"

## Selected variables
availableVariables = c("tas", "tasmin", "tasmax", "pr")

selectedVariables = c("tas", "tasmin", "tasmax", "pr")

if (any(!(selectedVariables %in% availableVariables)))
{
	warning(paste0("The following variables are not available and are removed:\n- ",
		paste0(selectedVariables[!(selectedVariables %in% availableVariables)], collapse = "\n- ")))

	selectedVariables = selectedVariables[selectedVariables %in% availableVariables]

	if (length(selectedVariables) == 0)
		stop("No selected variable recognised")
}

## Selected years
availableYears = 1980:2018 # January 1979 is missing for tas variable, and year 2019 is missing for pr variable

selectedYears = 2010

if (any(!(selectedYears %in% availableYears)))
{
	warning(paste0("The following years are not available and are removed:\n- ",
		paste0(selectedYears[!(selectedYears %in% availableYears)], collapse = "\n- ")))

	selectedYears = selectedYears[selectedYears %in% availableYears]

	if (length(selectedYears) == 0)
		stop("No selected year available")
}

#### Quick and dirty download for monthly timeseries for the whole world
n_years = length(selectedYears)
n_months = 12
n_iter = n_years*n_months

for (var in selectedVariables)
{
	url_file = paste0(downloadLink, var, "/")
	
	exitDir = paste0(folderDownloadData, var, "/")
	if (!dir.exists(exitDir))
		dir.create(path = exitDir)

	for (id in 1:n_iter)
	{
		info = getYearMonth(id = id, years = selectedYears)
		year = info[["year"]]
		month = info[["month"]]
		file = paste0("CHELSA_", var, "_", ifelse(month < 10, "0", ""), month, "_", year, "_V.2.1.tif")
		download.file(url = paste0(url_file, file), destfile = paste0(exitDir, file))
	}
	print(paste("Variable ", var, "done"))
}
