
#### Aim of prog: Plot and analyse the ouptuts of Cpp prog
## Plot densities (initial and last states)
# sumTrunkArea is in m^2

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)
library(terra)
# library(raster)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

options(max.print = 500)

#### Tool functions
## Function to clean strings
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

setorigin = list(plot_id = 2975)
setorigin = list(row = 270, col = 5)
plot_tw = function(transect, formatPlot, plotInfos, maxDistance = NULL, subsetIter = NULL, setorigin = NULL)
{
	# Common variables
	speciesList = plotInfos[["speciesList"]]
	species = plotInfos[["species"]]
	fromSouth = plotInfos[["fromSouth"]]

	nIter = plotInfos[["nIter"]]
	transect_index = plotInfos[["transect_index"]]

	landscapeSize = plotInfos[["landscapeSize"]]
	nRow_land = plotInfos[["nRow_land"]]
	nCol_land = plotInfos[["nCol_land"]]
	initOption = plotInfos[["initOption"]]
	climateRegion = plotInfos[["climateRegion"]]
	delta_t = plotInfos[["delta_t"]]

	maxLat_dist = plotInfos[["maxLat_dist"]]

	compute_tw = plotInfos[["compute_tw"]]
	compute_asymSpeed = plotInfos[["compute_asymSpeed"]]
	
	kernelType = plotInfos[["kernelType"]]

	coloursVec = c("#071B1B", "#135255", "#637872", "#B2EF80", "#F7DFC0", "#CFA47D", "#E28431")
	count = 1
	iterToPlot = round(seq(0, nIter - 1, length.out = length(coloursVec) + 1)) # +1 coming from the first plot (black curve)

	if (!is.null(subsetIter))
	{
		if (length(subsetIter) > length(coloursVec) + 1) # +1 coming from the first plot (black curve)
			stop("Dimensions between subsetIter and coloursVec mismatch")
		
		if (max(subsetIter) > length(coloursVec) + 1) # +1 coming from the first plot (black curve)
			stop("subsetIter goes beyond the length of coloursVec")
		
		coloursVec = coloursVec[subsetIter]
		iterToPlot = iterToPlot[subsetIter]
	}

	origin_dist = NULL
	if (!is.null(setorigin))
	{
		if (!is.list(setorigin))
			stop("setorigin must be a list containing one case among 1: (row, col), 2: distToNorth, 3: distToSouth, 4: plot_id")

		if (length(setorigin) == 2)
		{
			keyName = names(setorigin)
			if (!all(keyName %in% c("row", "col")))
				stop("setorigin must contain row and col of origin")

			row_id = setorigin[["row"]]
			col_id = setorigin[["col"]]

			if (row_id > nRow_land)
				stop("row exceed landscape dimensions")

			if (col_id > nCol_land)
				stop("col exceed landscape dimensions")

			plot_id_orig = row_id * nCol_land + col_id
		}

		if (length(setorigin) == 1)
		{
			keyName = names(setorigin)
			if (keyName != "plot_id")
				stop("setorigin must contain plot_id")
			
			plot_id_orig = setorigin[["plot_id"]]
		}

		if (!transect[transectOrigin == ls_origin[transect_index], plot_id_orig %in% patch_id])
			stop("The plot_id does not belong to the selected transect!")

		origin_dist = transect[patch_id == plot_id_orig, unique(distance)]
		if (!fromSouth[species])
			origin_dist = maxLat_dist - origin_dist
	}
	
	yMax = transect[(iteration %in% iterToPlot) & (transectOrigin == ls_origin[transect_index]) & !is.na(basalArea) &
		if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0, max(basalArea)]

	if (!(formatPlot %in% c("tex", "pdf", "pdf-tex")))
	{
		warning("Format not coded, set to pdf by default")
		formatPlot = "pdf"
	}

	if (!is.null(maxDistance))
		transect = transect[distance <= maxDistance]

	# Plot
	if (formatPlot == "pdf")
	{
		pdf(paste0("travellingWave_", species, "_", landscapeSize,"_", initOption, "_", kernelType, "_", climateRegion, ".pdf"),
			width = 17, height = 6)
	} else if (formatPlot == "tex") {
		tikz(paste0("travellingWave_", species, "_", landscapeSize,"_", initOption, "_", kernelType, "_", climateRegion, ".tex"),
			width = plotInfos[["width"]], height = plotInfos[["height"]])
	} else {
		for (format_pdfTex in c("pdf", "tex"))
			plot_tw(transect, format_pdfTex, plotInfos, maxDistance, subsetIter, setorigin) # Recursive call
		return (0)
	}

	op = par(mar = c(4, 4, 0.8, 5.5), mgp = c(2.8, 0.8, 0), tck = -0.02, xpd = TRUE)
	if (!is.null(origin_dist))
	{
		plot(transect[(iteration == iterToPlot[1]) & (transectOrigin == ls_origin[transect_index]) &
				!is.na(basalArea) & if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0,
				if (fromSouth[species]) distance - origin_dist else maxLat_dist - distance - origin_dist],
			transect[(iteration == iterToPlot[1]) & (transectOrigin == ls_origin[transect_index]) & !is.na(basalArea) &
				if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0, basalArea],
			type = "l", ylim = c(0, 1.01*yMax), las = 1,
			xlab = "Distance from leading edge (orange point in Fig. \\ref{fig::initLandscape})", ylab = "Basal area", lwd = 2)
	} else {
		plot(transect[(iteration == iterToPlot[1]) & (transectOrigin == ls_origin[transect_index]) &
				!is.na(basalArea) & if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0,
				if (fromSouth[species]) distance else maxLat_dist - distance],
			transect[(iteration == iterToPlot[1]) & (transectOrigin == ls_origin[transect_index]) & !is.na(basalArea) &
				if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0, basalArea],
			type = "l", ylim = c(0, 1.01*yMax), las = 1,
			xlab = "Distance from leading edge (orange point in Fig. \\ref{fig::initLandscape})", ylab = "Basal area", lwd = 2)
	}

	## For loop on time
	for (i in iterToPlot[2:length(iterToPlot)])
	{
		if (count > length(coloursVec))
			print("*** Warning, curve will not be plotted because of undefined colour")

		if (!is.null(origin_dist))
		{
			lines(transect[(iteration == i) & (transectOrigin == ls_origin[transect_index]) & !is.na(basalArea) &
					if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0,
					if (fromSouth[species]) distance - origin_dist else maxLat_dist - distance - origin_dist],
				transect[(iteration == i) & (transectOrigin == ls_origin[transect_index]) & !is.na(basalArea) &
					if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0, basalArea],
				lwd = 2, col = coloursVec[count])
		count = count + 1
		} else {
			lines(transect[(iteration == i) & (transectOrigin == ls_origin[transect_index]) & !is.na(basalArea) &
					if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0,
					if (fromSouth[species]) distance else maxLat_dist - distance],
				transect[(iteration == i) & (transectOrigin == ls_origin[transect_index]) & !is.na(basalArea) &
					if (fromSouth[species]) signedDistance >= 0 else signedDistance <= 0, basalArea],
				lwd = 2, col = coloursVec[count])
			count = count + 1
		}
	}

	legend(x = "topright", inset = if (formatPlot == "tex") c(-0.2, 0) else c(-0.065, 0),
		title = "Years", legend = round(iterToPlot*delta_t), lwd = 2, col = c("#000000", coloursVec), bty = "n")

	dev.off()

	## Plot speed on a 2nd graph, #! speed_dt is a global variable
	if (compute_tw[species])
	{
		if (smootherOption)
			smo = smooth.spline(x = speed_dt[!is.na(speed), year], y = speed_dt[!is.na(speed), speed], spar = 0.5)
		
		if (formatPlot == "pdf")
		{
			pdf(paste0("travellingWave_", species, "_", landscapeSize,"_", initOption, "_", kernelType, "_", climateRegion, "_speed",
				ifelse(smootherOption,"-smo", ""), ".pdf"), width = 6.25, height = 6.25/(17/6))
		} else {
			tikz(paste0("travellingWave_", species, "_", landscapeSize,"_", initOption, "_", kernelType, "_", climateRegion, "_speed",
				ifelse(smootherOption,"-smo", ""), ".tex"), width = plotInfos[["width"]], height = plotInfos[["height"]])
		}
		op = par(mar = c(4, 4, 0.8, 5.5), mgp = c(2.8, 0.8, 0), tck = -0.02)

		if (smootherOption)
		{
			plot(predict(smo, speed_dt[!is.na(speed), year]), type = "l", lwd = 2, xlab = "Year", ylab = "speed (m/yr)", las = 1)
		} else {
			plot(speed_dt[!is.na(speed), year], speed_dt[!is.na(speed), speed], type = "l",
				lwd = 2, xlab = "Year", ylab = "speed (m/yr)", las = 1)
		}

		if (compute_asymSpeed[species])
			abline(h = asymSpeed, lwd = 1.5, lty = "dashed", col = "#135255")

		dev.off()
	}
}

######## Part I: Population dynamics
#### Load parameters c++
## Simulation parameters
simulationParameters = setDT(read.table(file = "../run/simulationParameters.txt", header = FALSE,
	sep = "=", comment.char = "#", blank.lines.skip = TRUE))
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
speciesList = stringCleaner(speciesList, " ")

## Paths to files
pathCpp = "../"

pathSummary = paste0(pathCpp, simulationParameters[parameters == "summaryFilePath", values], speciesList, "/")
pathPopDyn = paste0(pathCpp, simulationParameters[parameters == "popDynFilePath", values], speciesList, "/")
initPath = paste0(pathCpp, simulationParameters[parameters == "initPath", values], speciesList, "/")

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

#! TEMPORARY ZONE
# #? Orford region, Acer saccharum alone, fat-tailed, no max disp, no min age
# pathSummary = "../run/results/withoutMaxDisp/withoutMinAge/acsa_alone_orford_fat-tail/summary/Acer_saccharum/"
# speciesList = "Acer_saccharum"
# compute_tw = TRUE
# compute_asymSpeed = FALSE
# fromSouth = TRUE
# speciesList = "Acer_saccharum"
# names(compute_tw) = speciesList
# names(compute_asymSpeed) = speciesList
# names(fromSouth) = speciesList
# tmax = 1500
# nIter = 6751

# #? Orford region, Abies balsamea and Acer saccharum
# pathSummary = c("../run/results/tmax2000_nIter9001_abba-acsa_orford/summary/Abies_balsamea/",
# 	"../run/results/tmax2000_nIter9001_abba-acsa_orford/summary/Acer_saccharum/")
# initPath = c("../run/data/initialCondition/abba-acsa_300x11_Orford/Abies_balsamea/",
# 	"../run/data/initialCondition/abba-acsa_300x11_Orford/Acer_saccharum/")
# speciesList = c("Abies_balsamea", "Acer_saccharum")

#? New Jersey region, Abies balsamea and Acer saccharum, restarted simulation (i.e., two parts)
# pathSummary = c("../run/results/mergedResults/summary/Abies_balsamea/",
# 	"../run/results/mergedResults/summary/Acer_saccharum/")
# initPath = c("../run/data/initialCondition/abba-acsa_300x11_newJersey/Abies_balsamea/",
# 	"../run/data/initialCondition/abba-acsa_300x11_newJersey/Acer_saccharum/")
# speciesList = c("Abies_balsamea", "Acer_saccharum")

# nIter = 13501
# tmax = 3000

#! END TEMPORARY ZONE

#### Common variables
## Check directories
if (any(!dir.exists(initPath)))
	stop(paste0("*** Folder <", initPath[!dir.exists(initPath)], "> does not exist ***"))

if (any(!dir.exists(pathSummary)))
	stop(paste0("*** Folder <", pathSummary[!dir.exists(pathSummary)], "> does not exist ***"))

## Define for which species the travelling wave speed should be computed
# Travelling wave speed
# compute_tw = TRUE
# compute_tw = FALSE
compute_tw = c(FALSE, TRUE)
names(compute_tw) = speciesList

# Asymptotic speed of travelling wave (it does not make sense for accelerating travelling waves for instance)
# compute_asymSpeed = FALSE
compute_asymSpeed = c(FALSE, FALSE)
names(compute_asymSpeed) = speciesList

## Define original position of the species and select transect
# fromSouth = TRUE
# fromSouth = FALSE
fromSouth = c(FALSE, TRUE)
names(fromSouth) = speciesList
maxLat_dist = (nRow_land - 1)*deltaLat # maxLat_dist - distance = distance to the south for a northern species

transect_index = 6

if (length(compute_tw) != length(speciesList))
	stop("Dimensions mismatch between `compute_tw` and species list")

if (length(compute_asymSpeed) != length(speciesList))
	stop("Dimensions mismatch between `compute_asymSpeed` and species list")

if (length(fromSouth) != length(speciesList))
	stop("Dimensions mismatch between `fromSouth` and species list")

## Plotting options and informations
landscapeSize = paste0(nRow_land, "x", nCol_land)
climateRegion = "NewJersey" # "Orford", "NewJersey"
smootherOption = TRUE

formatPlot = "pdf-tex" # "pdf", "tex", "pdf-tex"

plotInfos = vector(mode = "list", length = 14)
names(plotInfos) = c("speciesList", "species", "fromSouth", "nIter", "transect_index", "landscapeSize", "initOption",
	"climateRegion", "delta_t", "maxLat_dist", "compute_tw", "compute_asymSpeed", "width", "height")

plotInfos[["speciesList"]] = speciesList
plotInfos[["fromSouth"]] = fromSouth

plotInfos[["nIter"]] = nIter
plotInfos[["transect_index"]] = transect_index

plotInfos[["landscapeSize"]] = landscapeSize
plotInfos[["nRow_land"]] = nRow_land
plotInfos[["nCol_land"]] = nCol_land
plotInfos[["climateRegion"]] = climateRegion
plotInfos[["delta_t"]] = delta_t

plotInfos[["maxLat_dist"]] = maxLat_dist

plotInfos[["compute_tw"]] = compute_tw
plotInfos[["compute_asymSpeed"]] = compute_asymSpeed

plotInfos[["width"]] = 6
plotInfos[["height"]] = 4

plotInfos[["kernelType"]] = "2Dt_abba-vs-acsa" # "2Dt", "gaussian", "2Dt_abba-vs-acsa"

#### Loop over species to load C++ results
for (species in speciesList)
{
	## Species-specific path
	initPath_current = initPath[stri_detect(str = initPath, regex = species)]
	pathSummary_current = pathSummary[stri_detect(str = pathSummary, regex = species)]
	
	## Load initial condition
	# List files
	ls_init = list.files(initPath_current)
	initOption = paste0(length(ls_init), "RowsInit")

	# Determine their c++ coordinates (starting from 0 to n-1 rather than 1 to n)
	init_index = sort(as.integer(stri_sub(ls_init, from = stri_locate_last(ls_init, fixed = "_")[, "end"] + 1,
		to = stri_locate_last(ls_init, fixed = ".txt")[, "start"] - 1)))

	init_col = unique(init_index %% nCol_land)
	init_row = unique((init_index - init_col)/nCol_land)

	## Load files belonging to same transect (i.e., either same row = East-West or col = North-South)
	# Prepare results data tables
	nbData_ns = nRow_land*nIter*length(init_col) # nb of row x nIter x number of cols to cover

	transect_ns = data.table(patch_id = integer(length = nbData_ns), iteration = numeric(length = nbData_ns),
		localSeedProduced = numeric(length = nbData_ns), height_star = numeric(length = nbData_ns),
		sumTrunkArea = numeric(length = nbData_ns), totalDensity = numeric(length = nbData_ns),
		distance = numeric(length = nbData_ns), transectOrigin = integer(length = nbData_ns))

	ls_origin = c()	
	
	# Loop on the North-South transects, distance is from the northeast point of the transect
	for (i in 1:length(init_col))
	{
		ind_start = (i - 1)*nRow_land*nIter + 1
		ind_end = (i - 1)*nRow_land*nIter + nIter
	
		# --- Detect the northeast point colonised in column col_id
		col_id = i - 1 # C++ starts at 0
		indicesCol = seq(col_id, (nRow_land - 1)*nCol_land, nCol_land)
		currentOrigin = min(init_index[init_index %in% indicesCol]) # min because raster starts numbering from north, lower = more north
		currentRow = (currentOrigin - init_col[i])/nCol_land

		ls_origin = c(ls_origin, currentOrigin)

		# --- Loading files
		for (row in 0:(nRow_land - 1))
		{
			patch_id = init_col[i] + row*nCol_land
			temporary = fread(paste0(pathSummary_current, summaryPattern, patch_id, ".txt"))[1:nIter]
			transect_ns[ind_start:ind_end, patch_id := ..patch_id]
			transect_ns[ind_start:ind_end, c("iteration", "localSeedProduced", "height_star", "sumTrunkArea", "totalDensity") := temporary]
			transect_ns[ind_start:ind_end, distance := abs(currentRow - row)*deltaLat]
			transect_ns[ind_start:ind_end, signedDistance := (currentRow - row)*deltaLat]
			transect_ns[ind_start:ind_end, transectOrigin := currentOrigin]
			ind_start = ind_start + nIter
			ind_end = ind_end + nIter
		}
		print(paste(round(i*100/length(init_col), 2), "% done"))
	}

	if ((transect_index < 1) | (transect_index > length(ls_origin)))
		stop(paste0("The index transect_index must be between 1 and ", length(ls_origin), ". Currently, transect_index = ", transect_index))

	## Compute basal area
	transect_ns[, basalArea := sumTrunkArea/plotArea_ha]

	## Compute speed of travelling wave
	if (compute_tw[species])
	{
		# Common variables
		threshold_BA = 0.01 # Required basal area to consider a plot is populated

		# Subset data
		speed_dt = transect_ns[(transectOrigin == ls_origin[transect_index]) & (basalArea > threshold_BA) &
			(ifelse(fromSouth[species], signedDistance >= 0, signedDistance <= 0)), min(iteration, na.rm = TRUE), by = distance]

		setnames(speed_dt, new = c("distance", "iteration"))
		speed_dt = speed_dt[iteration > 0]
		
		setorder(speed_dt, -distance)
		speed_dt[, year := iteration*delta_t]
		speed_dt[, speed := c((distance[1:(.N - 1)] - distance[2:.N])/(year[1:(.N - 1)] - year[2:.N]), NA)]

		## Print results
		paste0("The averaged speed is: ", speed_dt[, round(mean(speed, na.rm = TRUE), 2)], " m/yr")
		paste0("The minimum speed is: ", speed_dt[, round(min(speed, na.rm = TRUE), 2)], " m/yr")
		paste0("The minimum positive speed is: ", speed_dt[speed >= 0, round(min(speed, na.rm = TRUE), 2)], " m/yr")

		## Compute asymptotic speed (I assume it is the most represented speed, since it should be at equilibrium)
		if (compute_asymSpeed[species])
		{
			asymSpeed = speed_dt[, .N, by = speed][which.max(N), speed]
			if (asymSpeed/speed_dt[1:round(.N/4), mean(speed)] > 1.05)
				warning("The asymptotic speed differs from the last quarter of speed_dt by more than 5%")
		}
	}

	## Plot travelling waves emanating from same origin, at different time 
	plotInfos[["initOption"]] = initOption
	plotInfos[["species"]] = species

	plot_tw(transect = transect_ns, formatPlot = formatPlot, plotInfos = plotInfos) # subsetIter = c(1, 4:8)

	print(paste(species, "done"))
}

# nbData_ew = nCol_land*nIter*length(init_row) # nb of row x nIter x number of rows to cover

# transect_ew = data.table(patch_id = integer(length = nbData_ew), iteration = numeric(length = nbData_ew),
# 	localSeedProduced = numeric(length = nbData_ew), height_star = numeric(length = nbData_ew),
# 	sumTrunkArea = numeric(length = nbData_ew), totalDensity = numeric(length = nbData_ew),
# 	distance = numeric(length = nbData_ew), transectOrigin = integer(length = nbData_ew))


#! ------------------------------------------------------------
#* ------------------------------------------------------------
#? WHAT FOLLOW IS A CRASH TEST ZONE CONTAINING MANY CRASH TESTS
#* ------------------------------------------------------------
#! ------------------------------------------------------------

#### ! CRASH TEST ZONE 1, ON RASTER
## Load everything
nbPatches = nRow_land*nCol_land
nbData = nbPatches*nIter
summary_dt = data.table(patch_id = integer(length = nbData), iteration = numeric(length = nbData),
	localSeedProduced = numeric(length = nbData), height_star = numeric(length = nbData),
	sumTrunkArea = numeric(length = nbData), totalDensity = numeric(length = nbData),
	row = numeric(length = nbData), col = integer(length = nbData))

for (patch_id in 0:(nbPatches - 1))
{
	ind_start = patch_id*nIter + 1
	ind_end = (patch_id + 1)*nIter
	temporary = fread(paste0(pathSummary, "summary/", species, "/", summaryPattern, patch_id, ".txt"))
	summary_dt[ind_start:ind_end, patch_id := ..patch_id]
	summary_dt[ind_start:ind_end, c("iteration", "localSeedProduced", "height_star", "sumTrunkArea", "totalDensity") := temporary]
}

summary_dt[, basalArea := sumTrunkArea/plotArea_ha]

summary_dt[, col := patch_id %% nCol_land]
summary_dt[, row := (patch_id - col)/nCol_land]

summary_dt[, x := col*deltaLon]
summary_dt[, y := (nRow_land - 1 - row)*deltaLat]

indices = seq(1, nbData, nIter) # iter = 0
l_out = 10
iterLoaded = round(seq(1, nIter - 1, length.out = l_out))
summary_list_rs = vector(mode = "list", length = l_out + 1)

summary_list_rs[[1]] = rast(x = summary_dt[indices, .(x, y, basalArea, totalDensity, localSeedProduced, height_star)],
	type = "xyz", crs = "EPSG:3799") # iter = 0

for (i in 1:l_out)
{
	iter = iterLoaded[i]
	indices = seq(iter + 1, nbData, nIter)
	summary_list_rs[[i + 1]] = rast(x = summary_dt[indices, .(x, y, basalArea, totalDensity, localSeedProduced, height_star)],
		type = "xyz", crs = "EPSG:3799") # raster of 4 layers
}


crop_extent = ext(summary_list_rs[[1]])

pdf("0.pdf")
par(mfrow = c(1, 4))
plot(summary_list_rs[[1]][["basalArea"]], axes = FALSE, legend = FALSE)
plot(summary_list_rs[[2]][["basalArea"]], axes = FALSE, legend = FALSE)
plot(summary_list_rs[[5]][["basalArea"]], axes = FALSE, legend = FALSE)
plot(summary_list_rs[[8]][["basalArea"]], axes = FALSE, legend = FALSE)
dev.off()




toto = sds(summary_list_rs)





# Example of cropping
crop_extent = extent(c(500, 1500, 500, 1500))

# For not cropping
crop_extent = extent(summary_list_rs[[1]])
crop_extent = ext(summary_list_rs[[1]])

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

#### Plot density at different time
## Initiate plot


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
#### ! END CRASH TEST ZONE 2

#### ! CRASH TEST ZONE 3
kernel = function(d, p, u)
	return (p/(pi*u*(1 + d^2/u)^(p + 1)))

abba = function(x)
	return (kernel(x, 0.6372, 813.1))

acsa = function(x)
	return (kernel(x, 2.015, 24355.4))

kurtosis = function(p)
	return (2*(p - 1)/(p - 2))

ratio = function(x)
	return(abba(x)/acsa(x))

curve(ratio, 0, 10000, lwd = 2, col = "#F64905")
curve(abba, 2000, 2005, lwd = 2, col = "#E28431", ylim = c(0, 2.4e-10))
curve(acsa, 2000, 2005, lwd = 2, col = "#135255", add = TRUE)


patches = ll[(patch_id > 1599) & (patch_id < 2602), unique(patch_id)]
lim = length(patches)

ll[patch_id == 1600]

pdf("test.pdf", width = 18, height = 13)
plot(ll[patch_id == patches[1], iteration], ll[patch_id == patches[1], height_star], type = "l", lwd = 0.5, ylim = c(0, 35))
for (i in 2:lim)
	lines(ll[patch_id == patches[i], iteration], ll[patch_id == patches[i], height_star], type = "l", lwd = 0.5)
dev.off()

ll[(patch_id > 1599) & (patch_id < 2602) & (iteration > 4000), mean(height_star)]

# Compute diameter for an Abies balsamea of 31.72 m high
heightToDBH = function(h, a, b) # Give the dbh in mm from height in m
	return (10^(1/b*(b - a + log(h)/log(10))))

DBH_to_height = function(dbh, a, b)
	return (exp((a - b + b*log10(dbh))*log(10)));

heightToDBH(31.72, 0.4742, 0.5743)
R0 = 0.214645

heightToDBH(20.03259, 0.4742, 0.5743)
R0 = 3.43808

DBH_to_height(heightToDBH(31.72, 0.4742, 0.5743), 0.4742, 0.5743) # Should be identity