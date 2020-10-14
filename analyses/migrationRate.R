
#### Aim of prog: Plot and analyse the ouptuts of Cpp prog
## Plot densities (initial and last states)
# sumTrunkArea is in m^2

#! I SHOULD CREATE A LOOP FOR SPECIES OR CALL THIS PROG IN A FOR LOOP

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)
library(raster)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

options(max.print = 500)

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
plotArea_ha = plotArea/1e4

## Time
t0 = as.numeric(simulationParameters[parameters == "t0", values])
tmax = as.numeric(simulationParameters[parameters == "tmax", values])
nIter = as.integer(simulationParameters[parameters == "nIter", values])
delta_t = (tmax - t0)/(nIter - 1)

#### Load results c++
## Initial condition
# List files
ls_init = list.files(initPath)
# ls_init = c("ic_5049.txt", "ic_5050.txt") # In the case 100 x 100

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

#! TEMPORARY
pathSummary = "../cpp/summary/Acer_saccharum_100x100/"
pathSummary = "../cpp/summary/Acer_saccharum_100x5_TW/"
#! END TEMPORARY

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
		transect_ns[ind_start:ind_end, signedDistance := (currentRow - row)*deltaLat] # 'negative distance' means it is south to the origin
		transect_ns[ind_start:ind_end, transectOrigin := currentOrigin]
		ind_start = ind_start + nIter
		ind_end = ind_end + nIter
	}
}

## Compute basal area
transect_ns[, basalArea := sumTrunkArea/plotArea_ha]

#### Plot density at different time
## Initiate plot
pdf("travellingWave.pdf", width = 10, height = 2*10/3)
plot(transect_ns[(iteration == 1000) & (transectOrigin == init_index[1]), distance], transect_ns[(iteration == 1000) & (transectOrigin == init_index[1]), basalArea],
	type = "l", ylim = c(0, transect_ns[transectOrigin == init_index[1], max(basalArea)]),
	xlab = "Distance", ylab = "Basal area", lwd = 2)

## For loop on time
coloursVec = c("#071B1B", "#135255", "#637872", "#B2EF80", "#F7DFC0", "#CFA47D", "#E28431")
count = 1
for (i in (seq(1300, nIter, 300) - 1))
{
	lines(transect_ns[(iteration == i) & (transectOrigin == init_index[1]), distance], transect_ns[(iteration == i) & (transectOrigin == init_index[1]), basalArea],
		lwd = 2, col = coloursVec[count])
	count = count + 1
}
round(c(1000, (seq(1300, nIter, 300) - 1))*delta_t)
legend(x = "topright", legend = paste(round(c(1000, (seq(1300, nIter, 300) - 1))*delta_t), "years"),
	lwd = 2, col = c("#000000", coloursVec), bty = "n")
dev.off()

## Dynamic plot for transectOrigin == 54 toward north
transect_ns54 = transect_ns[(transectOrigin == init_index[1]) & (signedDistance >= 0)]

# Limits xlim and ylim
max_distance = transect_ns54[, max(distance)]
max_basalArea = transect_ns54[, max(basalArea)]

# Plot
plot(transect_ns54[iteration == 0, distance], transect_ns54[iteration == 0, basalArea],
	type = "l", xlim = c(0, max_distance), ylim = c(0, max_basalArea))
invisible(sapply(seq(2, nIter, by = 10), function(x, parameters) {
	Sys.sleep(0.01)
	plot(transect_ns54[iteration == x, distance], transect_ns54[iteration == x, basalArea], type = 'l', lwd = 1,
		xlim = c(0, max_distance), ylim = c(0, max_basalArea))
}, parameters = 0)) # unused parameter, I just put it to remember how to do it in case of I need it

#### Compute the speed
## Threshold basal area o consider a plot populated
threshold_BA = 1
speed_dt = transect_ns54[basalArea > threshold_BA, min(iteration), by = distance]
setnames(speed_dt, new = c("distance", "iteration"))
setorder(speed_dt, -distance)
speed_dt[, year := iteration*delta_t]
speed_dt[, speed := c((distance[1:(.N - 1)] - distance[2:.N])/(year[1:(.N - 1)] - year[2:.N]), NA)]

#### ! CRASH TEST ZONE 1, ON RASTER
## Load everything
nbPatches = nRow_land*nCol_land
nbData = nbPatches*nIter
summary_dt = data.table(patch_id = integer(length = nbData), iteration = numeric(length = nbData),
	localSeedProduced = numeric(length = nbData), localSeedBank = numeric(length = nbData),
	sumTrunkArea = numeric(length = nbData), totalDensity = numeric(length = nbData),
	row = numeric(length = nbData), col = integer(length = nbData))

pathSummary = "../cpp/summary/Acer_saccharum/"

for (patch_id in 0:(nbPatches - 1))
{
	ind_start = patch_id*nIter + 1
	ind_end = (patch_id + 1)*nIter
	temporary = fread(paste0(pathSummary, summaryPattern, patch_id, ".txt"))
	summary_dt[ind_start:ind_end, patch_id := ..patch_id]
	summary_dt[ind_start:ind_end, c("iteration", "localSeedProduced", "localSeedBank", "sumTrunkArea", "totalDensity") := temporary]
}

summary_dt[, basalArea := sumTrunkArea/plotArea_ha]

summary_dt[, col := patch_id %% nCol_land]
summary_dt[, row := (patch_id - col)/nCol_land]

summary_dt[, x := col*deltaLon]
summary_dt[, y := (nRow_land - 1 - row)*deltaLat]

coordinates(summary_dt) = ~ x + y
gridded(summary_dt) = TRUE

indices = seq(1, nbData, nIter) # iter = 0
l_out = 10
iterLoaded = round(seq(1, nIter - 1, length.out = l_out))
summary_list_rs = vector(mode = "list", length = l_out + 1)
summary_list_rs[[1]] = raster(summary_dt[indices, "basalArea"]) # iter = 0

for (i in 1:l_out)
{
	iter = iterLoaded[i]
	indices = seq(iter + 1, nbData, nIter)
	summary_list_rs[[i + 1]] = raster(summary_dt[indices, "basalArea"])
}

# Example of cropping
crop_extent = extent(c(500, 1500, 500, 1500))

# For not cropping
crop_extent = extent(summary_list_rs[[1]])

pdf("0.pdf")
plot(crop(summary_list_rs[[1]], crop_extent), legend = FALSE)
BA_range = c(minValue(summary_list_rs[[1]]), maxValue(summary_list_rs[[1]]))
plot(summary_list_rs[[1]], legend.only = TRUE,
	legend.width = 1, legend.shrink = 0.75,
	axis.args = list(at = seq(trunc(BA_range[1]), round(BA_range[2]), 5),
		labels = seq(trunc(BA_range[1]), round(BA_range[2]), 5),
		cex.axis = 1),
	legend.args = list(text = "Basal area", side = 4, font = 2, line = 2.5, cex = 1))
dev.off()

pdf(paste0(iterLoaded[9], ".pdf"))
plot(crop(summary_list_rs[[9]], crop_extent), legend = FALSE)
BA_range = c(minValue(summary_list_rs[[9]]), maxValue(summary_list_rs[[9]]))
plot(summary_list_rs[[9]], legend.only = TRUE,
	legend.width = 1, legend.shrink = 0.75,
	axis.args = list(at = seq(trunc(BA_range[1]), round(BA_range[2]), 5),
		labels = seq(trunc(BA_range[1]), round(BA_range[2]), 5),
		cex.axis = 1),
	legend.args = list(text = "Basal area", side = 4, font = 2, line = 2.5, cex = 1))
dev.off()

pdf(paste0(iterLoaded[10], ".pdf"))
plot(crop(summary_list_rs[[10]], crop_extent), legend = FALSE)
BA_range = c(minValue(summary_list_rs[[10]]), maxValue(summary_list_rs[[10]]))
plot(summary_list_rs[[10]], legend.only = TRUE,
	legend.width = 1, legend.shrink = 0.75,
	axis.args = list(at = seq(trunc(BA_range[1]), round(BA_range[2]), 5),
		labels = seq(trunc(BA_range[1]), round(BA_range[2]), 5),
		cex.axis = 1),
	legend.args = list(text = "Basal area", side = 4, font = 2, line = 2.5, cex = 1))
dev.off()


#### ! CRASH TEST ZONE 2
aa = fread("../cpp/popDyn/Acer_saccharum/pd_54.txt")
# aa[iteration > 125 & iteration < 128]
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
plot(cc[125:135, iteration], cc[125:135, height_star])

aa = aa[(iteration == iterationBirth) & (iteration > 0)]
plot(aa$iterationBirth, aa$height)
aa = aa[density > 1/plotArea]

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
# op = par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
# plot(compReprod$time, compReprod$basalArea, type = "l", lwd = 2,
# 	xlab = "Time", ylab = "Basal area")
# dev.off()
