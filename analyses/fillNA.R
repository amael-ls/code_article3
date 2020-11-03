
#### Aim of prog: fill with NA missing data from C++ program.
## Comment:
# I suspect sometimes the C++ program compiled with optimisation level -O3
# has a problem in saving population dynamics at each time step.
# This is probably due to a data race in the file access, although I am not
# sure. In any case, I count the number missing line, and fill them with NA
#
# Actually I think there is no problem with the code. Because the computer
# had a problem, I guess that some iteration skipped. Indeed, only the last
# iterations had a problem before the computer crashed.

#### Load packages
library(data.table)
library(stringi)

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
