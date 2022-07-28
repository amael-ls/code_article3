
#### Aim of prog: Generate the initial condition from a previous run
# There are two things to save:
#	1. The new initial condition (1 file per plot per species)
#	2. The new populated patches file. Saved in the landscape directory, but no overwrite of the previous file

#### Load library and clear memory
library(data.table)
library(stringi)

rm(list = ls())
graphics.off()

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

## Function to print initial condition
printIC = function(dt, path, filenamePattern, sep = " ", reset = TRUE)
{
	id_plots = dt[, unique(patch_id)]
	nbPlots = length(id_plots)
	for (plot in 1:nbPlots)
	{
		outfileName = paste0(path, filenamePattern, id_plots[plot], ".txt")

		if (reset & file.exists(outfileName))
			file.remove(outfileName)

		ofstream = file(outfileName)

		line = paste0("density", sep, "dbh", sep = "\n")

		paste0(line, paste0(dt[patch_id == plot, density], sep, dt[patch_id == plot, dbh], collapse = "\n"))

		writeLines(line, ofstream)
		close(ofstream)
	}
}

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

pathPopDyn = paste0(pathCpp, simulationParameters[parameters == "popDynFilePath", values], speciesList, "/")
path_ic = paste0(pathCpp, simulationParameters[parameters == "initPath", values], speciesList, "/")

## Files' patterns
popDynPattern = simulationParameters[parameters == "popDynFilePattern", values]
initPattern = simulationParameters[parameters == "initFilenamePattern", values]

## Iteration infos
nIter = as.integer(simulationParameters[parameters == "nIter", values])
saveFrom = as.integer(simulationParameters[parameters == "saveFrom", values])
freqSave = as.integer(simulationParameters[parameters == "freqSave", values])
saveOnlyLast = as.logical(simulationParameters[parameters == "saveOnlyLast", values])

## Landscape metadata
landscape_metadata = setDT(read.table(file = paste0(pathCpp, simulationParameters[parameters == "climate_file", values]),
	header = FALSE, sep = "=", blank.lines.skip = TRUE, fill = TRUE))

if (ncol(landscape_metadata) == 3)
{
	landscape_metadata[, V3 := NULL]
	landscape_metadata[V1 == "delimiter", V2 := "="]
}

setnames(landscape_metadata, new = c("parameters", "values"))

landscape_metadata[, values := stringCleaner(values, fixed = "./", skip = "../"), by = parameters]
pathLandscape = paste0(pathCpp, landscape_metadata[parameters == "path", values])

dimLandscape = as.integer(landscape_metadata[parameters == "nRow", values]) * as.integer(landscape_metadata[parameters == "nCol", values])

## Check existence directories
if (any(!dir.exists(path_ic)))
	stop(paste0("*** Directory <", path_ic[!dir.exists(path_ic)], "> does not exist ***"))

if (any(!dir.exists(pathPopDyn)))
	stop(paste0("*** Directory <", pathPopDyn[!dir.exists(pathPopDyn)], "> does not exist ***"))

if (any(!dir.exists(pathLandscape)))
	stop(paste0("*** Directory <", pathLandscape[!dir.exists(pathLandscape)], "> does not exist ***"))

#### Create initial condition from previous run
## Output directories
outputDir = paste0(stri_replace_last(str = path_ic, replacement = "", regex = "/"), "_from_prev_run/")

## Common variables
# Select the iteration that will be the initial condition
selectedIteration = 5200
if (selectedIteration < saveFrom)
	stop(paste("You must select an iteration >=", saveFrom))

if (selectedIteration > nIter - 1)
	stop(paste("You must select an iteration <=", nIter))

if (selectedIteration %% freqSave)

if ((saveOnlyLast) & (selectedIteration != nIter - 1) & (selectedIteration != 0))
	stop(paste("Only the first and last iterations were saved. You must choose between 0 and", nIter - 1))

# Save the patches that are colonised
populatedPatches_ls = vector(mode = "list", length = length(speciesList))
names(populatedPatches_ls) = speciesList

## Loop over species
for (species in speciesList)
{
	# Species specific path data
	path_data = pathPopDyn[stri_detect(str = pathPopDyn, regex = paste0(species, "/$"))]
	out_data = outputDir[stri_detect(str = outputDir, regex = species)]
	
	# List population dynamics file
	ls_files = list.files(path = path_data, pattern = paste0("^", popDynPattern, ".*.txt$"))
	if (length(ls_files) != dimLandscape)
		stop("Dimensions between landscape size and listed files mismatch")
	
	# Load files
	popDyn_ls = vector(mode = "list", length = dimLandscape)
	for (i in 1:dimLandscape)
	{
		filename = paste0(path_data, ls_files[i])
		popDyn_ls[[i]] = fread(filename)[iteration == selectedIteration]
		if (i %% 100 == 0)
			print(paste0(round(i*100/dimLandscape, digits = 3), "% done"))
	}
	popDyn_ls = rbindlist(popDyn_ls, idcol = "patch_id")

	# List colonised patches for initial condition
	populatedPatches_ls[[species]] = popDyn_ls[, .(unique(patch_id - 1))]
	setnames(populatedPatches_ls[[species]], new = "patch_id")

	# Create output directory
	if (!dir.exists(out_data))
		dir.create(out_data)

	if (length(list.files(out_data, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)) != 0)
		unlink(paste0(out_data, "*"))

	printIC(popDyn_ls, out_data, initPattern)
}

populatedPatches = rbindlist(populatedPatches_ls, idcol = "species")

saveRDS(populatedPatches, paste0(pathLandscape, "populatedPatches_restart.rds"))
