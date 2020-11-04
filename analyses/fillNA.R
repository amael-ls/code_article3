
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

#! TEMPORARY ZONE
pathSummary = "../cpp/summary/Acer_saccharum_100x100_10RowsInit/"
nIter = 1997
delta_t = 1000/2999
tmax = (nIter - 1)*delta_t
#! END TEMPORARY ZONE

#### Load results c++
## List files
ls_files = list.files(pathSummary)

## Record plot_id and iteration missing
plot_id_rec = c()
iteration_rec = c()

## Load files and fill NA
allIterations = 0:(nIter - 1)
for (file in ls_files)
{
	summaryFile = fread(paste0(pathSummary, file))
	if (summaryFile[, .N] != nIter)
	{
		missing = allIterations[!(allIterations %in% summaryFile[, iteration])]
		plot_id_rec = c(plot_id_rec, rep(file, length(missing)))
		iteration_rec = c(iteration_rec, missing)

		summaryFile = rbind(summaryFile, data.table(iteration = missing, localSeedProduced = NA, heigh_star = NA, sumTrunkArea = NA, totalDensity = NA))
		setorder(x = summaryFile, iteration)
		fwrite(x = summaryFile, file = paste0(pathSummary, file), sep = " ", col.names = TRUE, na = "NA")
	}
}

missingLines = data.table(plot_id = as.integer(stri_sub(plot_id_rec, from = stri_locate(str = plot_id_rec, regex = "_")[, "end"] + 1,
		to = stri_locate(str = plot_id_rec, regex = ".txt")[, "start"] - 1)),
	iteration = iteration_rec)

print(paste0("Minimum iteration: ", missingLines[, min(iteration)]))

hist(missingLines[, iteration])

## Check that those whom are missing iterations before nIter - 1 are actually missing up to nIter - 1
missingLines[, freq := .N, by = plot_id]
missingLines[, minIter := min(iteration), by = plot_id]

missingLines[minIter + freq != nIter]
