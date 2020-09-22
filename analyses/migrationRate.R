
#### Aim of prog: Plot and analyse the ouptuts of Cpp prog
## Plot densities (initial and last states)
#

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

#### Tool function
stringCleaner = function(str, fixed, skip = NULL)
{
	if (is.null(skip))
		return (stri_replace(str = str, replacement = "", fixed = fixed));
	
	skipPos = stri_locate_all(str, fixed = skip)[[1]]
	fixedPos = stri_locate_all(str, fixed = fixed)[[1]]
	fixedInSkip = stri_locate_all(skip, fixed = fixed)[[1]]
	
	if (anyNA(skipPos))
		return (stri_replace(str = str, replacement = "", fixed = fixed));
	
	if (fixedInSkip[, "start"] == 1)
		fixedPos = fixedPos[!(fixedPos[, "start"] %in% skipPos[, "start"]), ]

	if (fixedInSkip[, "end"] == stri_length(skip))
		fixedPos = fixedPos[!(fixedPos[, "end"] %in% skipPos[, "end"]), ]

	if (nrow(fixedPos) == 0) # All the fixed must be skipped
		return (str);

	strToReturn = stri_sub(str, from = 1, to = fixedPos[1, "start"] - 1) # if fixedPos[1, "start"] == 1, then it returns an empty string
	
	for (i in 1:nrow(fixedPos) - 1)
		strToReturn = paste0(strToReturn, stri_sub(str, from = fixedPos[i, "end"] + 1, to = fixedPos[i + 1, "start"] - 1))
	
	strToReturn = paste0(strToReturn, stri_sub(str, from = fixedPos[nrow(fixedPos), "end"] + 1, to = stri_length(str)))

	return (strToReturn);
}

######## Part I: Population dynamics
#### Load parameters c++
## Simulation parameters
simulationParameters = setDT(read.table(file = "../cpp/simulationParameters.txt", header = FALSE, sep = "=", comment.char = "#", blank.lines.skip = TRUE))
setnames(x = simulationParameters, new = c("parameters", "values"))

## Clean parameters
simulationParameters[, parameters := stringCleaner(parameters, " ")]
simulationParameters[, values := stringCleaner(values, " ")]
simulationParameters[, values := stringCleaner(values, fixed = "./", skip = "../"), by = parameters]

## Species
speciesList = simulationParameters[parameters == "species_filenames", values]

if (stri_detect(speciesList, regex = ","))
	speciesList = unlist(stri_split(str = speciesList, regex = ","))

speciesList = stringCleaner(speciesList, ".txt")

## Paths to files
pathCpp = "../cpp/"

pathSummary = paste0(pathCpp, simulationParameters[parameters == "summaryFilePath", values], speciesList, "/")
pathPopDyn = paste0(pathCpp, simulationParameters[parameters == "popDynFilePath", values], speciesList, "/")
initPath = paste0(simulationParameters[parameters == "initPath", values], speciesList, "/")

## Landscape metadata
landscape_metadata = setDT(read.table(file = paste0(pathCpp, simulationParameters[parameters == "climate_file", values]),
	header = FALSE, sep = "=", blank.lines.skip = TRUE, fill = TRUE))

if (ncol(landscape_metadata) == 3)
{
	landscape_metadata[, V3 := NULL]
	landscape_metadata[V1 == "delimiter", V2 := "="]
}

setnames(landscape_metadata, new = c("parameters", "values"))

nRow_land = as.integer(landscape_metadata[parameters == "nRow", values])
nCol_land = as.integer(landscape_metadata[parameters == "nCol", values])

## Time
t0 = as.numeric(simulationParameters[parameters == "t0", values])
tmax = as.numeric(simulationParameters[parameters == "tmax", values])
nIter = as.integer(simulationParameters[parameters == "nIter", values])
delta_t

#### Load results c++
## Initial condition
# List files
ls_init = list.files(initPath)

# Determine their c++ coordinates (starting from 0 to n-1 rather than 1 to n)
init_index = as.integer(stri_sub(ls_init, from = stri_locate_last(ls_init, fixed = "_")[, "end"] + 1,
	to = stri_locate_last(ls_init, fixed = ".txt")[, "start"] - 1))

init_col = init_index %% nCol_land
init_row = (init_index - init_col)/nCol_land

## Load files belonging to same transect (i.e., either same row = East-West or col = North-South)
nbData = nRow_land*
results = data.table()

#! RESTART HERE
# init = fread(paste0(initPath, "init.txt"))

# ## Last state
# end = fread(paste0(pathCpp, "end.txt"))

# ## Competition, reproduction, basal area and total density
# compReprod = fread(paste0(pathCpp, "compReprod.txt"))

# ## Dynamics
# # Load the density and dbh dynamics of cohorts
# dyn = fread(paste0(pathCpp, "popDyn.txt"))

# # Determine the number of cohorts and time steps
# nbTimeSteps = dim(compReprod)[1]
# nbCohorts = dim(dyn)[1]/nbTimeSteps

# #### Plots
# ## Initial state
# plot(init$dbh, init$density)

# ## Last state
# plot(end$dbh, end$density)

# ## Dynamic plot
# # Limits xlim and ylim
# max_dbh = dyn[, max(dbh)]
# max_density = dyn[, max(density)]

# # Plot
# plot(dyn[1:nbCohorts, dbh], dyn[1:nbCohorts, density], pch = '', # type = 'l',
# 	xlim = c(0, max_dbh), ylim = c(0, max_density))
# invisible(sapply(2:nbTimeSteps, function(x, nbCohorts) {
# 	Sys.sleep(0.005)
# 	ind_start = (x - 1)*nbCohorts + 1
# 	ind_end = x*nbCohorts
# 	plot(dyn[ind_start:ind_end, dbh], dyn[ind_start:ind_end, density], type = 'l', lwd = 1,
# 		xlim = c(0, max_dbh), ylim = c(0, max_density))
# }, nbCohorts = nbCohorts))

# ## Reproduction and competition
# par(mfrow = c(4,1))
# plot(compReprod$time, compReprod$reproduction, type = "l", lwd = 2,
# 	xlab = "Time", ylab = "Reproduction")
# plot(compReprod$time, compReprod$competition, type = "l", lwd = 2,
# 	xlab = "Time", ylab = "Competition")
# plot(compReprod$time, compReprod$basalArea, type = "l", lwd = 2,
# 	xlab = "Time", ylab = "Basal area")
# plot(compReprod$time, compReprod$totalDensity, type = "l", lwd = 2,
# 	xlab = "Time", ylab = "Total density")

# #### Tikz version
# tikz(paste0("./ba_plot.tex"), width = 3.1, height = 3.1) #, standAlone = TRUE)
# op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
# plot(compReprod$time, compReprod$basalArea, type = "l", lwd = 2,
# 	xlab = "Time", ylab = "Basal area")
# dev.off()
