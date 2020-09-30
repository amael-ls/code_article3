
#### Aim of prog: Plot and analyse the ouptuts of Cpp prog
## Plot densities (initial and last states)
#

#! I SHOULD CREATE A LOOP FOR SPECIES OR CALL THIS PROG IN A FOR LOOP

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

#### Tool functions
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

## Check cohorts are in decreasing order of dbh (within a species)
checkOrder = function(dbh)
{
	n = length(dbh)
	decreasingOrder = unique(dbh[1:(n-1)] - dbh[2:n] >= 0)

	if (length(decreasingOrder) == 1)
	{
		if (decreasingOrder)
			return (TRUE);
		
		return (FALSE);
	}
	return (FALSE);
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

## Files' patterns
summaryPattern = simulationParameters[parameters == "summaryFilePattern", values]
initPattern = simulationParameters[parameters == "initFilenamePattern", values] # Not used

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

deltaLon = as.numeric(landscape_metadata[parameters == "deltaLon", values])
deltaLat = as.numeric(landscape_metadata[parameters == "deltaLat", values])
plotArea = deltaLon*deltaLat

## Time
t0 = as.numeric(simulationParameters[parameters == "t0", values])
tmax = as.numeric(simulationParameters[parameters == "tmax", values])
nIter = as.integer(simulationParameters[parameters == "nIter", values])
delta_t = (tmax - t0)/(nIter - 1)

#### Load results c++
## Initial condition
# List files
ls_init = list.files(initPath)

# Determine their c++ coordinates (starting from 0 to n-1 rather than 1 to n)
init_index = as.integer(stri_sub(ls_init, from = stri_locate_last(ls_init, fixed = "_")[, "end"] + 1,
	to = stri_locate_last(ls_init, fixed = ".txt")[, "start"] - 1))

init_col = unique(init_index %% nCol_land)
init_row = unique((init_index - init_col)/nCol_land)

## Load files belonging to same transect (i.e., either same row = East-West or col = North-South)
# Prepare results data tables
nbData_ns = nRow_land*nIter*length(init_col) # nb of row x nIter x number of cols to cover
nbData_ew = nCol_land*nIter*length(init_row) # nb of row x nIter x number of rows to cover

transect_ns = data.table(patch_id = integer(length = nbData_ns), iteration = numeric(length = nbData_ns),
	localSeedProduced = numeric(length = nbData_ns), localSeedBank = numeric(length = nbData_ns),
	sumTrunkArea = numeric(length = nbData_ns), totalDensity = numeric(length = nbData_ns),
	distance = numeric(length = nbData_ns), transectOrigin = integer(length = nbData_ns))

transect_ew = data.table(patch_id = integer(length = nbData_ew), iteration = numeric(length = nbData_ew),
	localSeedProduced = numeric(length = nbData_ew), localSeedBank = numeric(length = nbData_ew),
	sumTrunkArea = numeric(length = nbData_ew), totalDensity = numeric(length = nbData_ew),
	distance = numeric(length = nbData_ew), transectOrigin = integer(length = nbData_ew))

# Loop on the North-South transects
for (i in 1:length(init_col))
{
	ind_start = (i - 1)*nRow_land*nIter + 1
	ind_end = (i - 1)*nRow_land*nIter + nIter
	currentOrigin = init_index[i] # This works only when length(init_col) == length(init_index)
	currentRow = (currentOrigin - init_col[i])/nCol_land

	for (row in 0:(nRow_land - 1))
	{
		patch_id = init_col[i] + row*nCol_land
		temporary = fread(paste0(pathSummary, summaryPattern, patch_id, ".txt"))
		transect_ns[ind_start:ind_end, patch_id := ..patch_id]
		transect_ns[ind_start:ind_end, c("iteration", "localSeedProduced", "localSeedBank", "sumTrunkArea", "totalDensity") := temporary]
		transect_ns[ind_start:ind_end, distance := abs(currentRow - row)*deltaLat] # This works only when there is one line in init_row
		transect_ns[ind_start:ind_end, transectOrigin := currentOrigin]
		ind_start = ind_start + nIter
		ind_end = ind_end + nIter
	}
}

## Compute basal area
transect_ns[, basalArea := sumTrunkArea/plotArea]

#### Plot density at different time
## Initiate plot
plot(transect_ns[(iteration == 0) & (transectOrigin == 54), distance], transect_ns[(iteration == 0) & (transectOrigin == 54), basalArea],
	type = "l", ylim = c(0, transect_ns[transectOrigin == 54, max(basalArea)]))

## For loop on time
for (i in (seq(100, nIter, 300) - 1))
	lines(transect_ns[(iteration == i) & (transectOrigin == 54), distance], transect_ns[(iteration == i) & (transectOrigin == 54), basalArea])


#### ! CRASH TEST ZONE
aa = fread("../cpp/popDyn/Acer_saccharum/pd_54.txt")
aa[, verifHeight := checkOrder(height), by = iteration]
aa[, verifDbh := checkOrder(dbh), by = iteration]

unique(aa[, verifHeight])
unique(aa[, verifDbh])
aa[, c("verifHeight", "verifDbh") := NULL]

aa[, nbCohorts := .N, by = iteration]
aa[, max(nbCohorts)]

aa[, sumTrunckArea := pi*sum(dbh*dbh*density)/(4*1000*1000), by = iteration] # trunk area in metres (hence the /1000*1000)
plotArea_ha = plotArea/1e4
aa[, basalArea := sumTrunckArea/plotArea_ha, by = iteration]
aa[, basalArea2 := pi*sum(dbh*dbh*density)/(4e2*plotArea), by = iteration]
bb = unique(aa[, .(iteration, basalArea)])
plot(bb$iteration, bb$basalArea)

plot(aa[iteration == 999, height], aa[iteration == 999, density], type = "l", lwd = 2, col = "#FF9933")
lines(aa[iteration == 500, height], aa[iteration == 500, density], lwd = 2, col = "#3366ff")

cc = fread("../cpp/summary/Acer_saccharum/su_54.txt")

aa = aa[(iteration == iterationBirth) & (iteration > 0)]
plot(aa$iterationBirth, aa$height)
aa = aa[density > 1/plotArea]


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
