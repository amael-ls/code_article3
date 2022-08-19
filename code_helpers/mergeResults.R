
#### Aim of prog: 

#### Load packages
library(data.table)
library(stringi)

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

## Function to print results
printRes = function(dt, path, filenamePattern, sep = " ", reset = TRUE)
{
	id_plots = dt[, unique(patch_id)]
	nbPlots = length(id_plots)
	for (plot in 1:nbPlots)
	{
		outfileName = paste0(path, filenamePattern, id_plots[plot], ".txt")

		if (reset & file.exists(outfileName))
			file.remove(outfileName)

		ofstream = file(outfileName)

		line = paste0("iteration", sep, "localSeedProduced", sep, "height_star", sep, "sumTrunkArea", sep,
			"totalDensity", sep = "\n")

		line = paste0(line, paste0(dt[patch_id == id_plots[plot], iteration], sep,
			dt[patch_id == id_plots[plot], localSeedProduced], sep,
			dt[patch_id == id_plots[plot], height_star], sep,
			dt[patch_id == id_plots[plot], sumTrunkArea], sep,
			dt[patch_id == id_plots[plot], totalDensity],
			collapse = "\n"))

		writeLines(line, ofstream)
		close(ofstream)
	}
}

#### Load results (summary only)
## Paths
first_run = "../run/results/withoutMaxDisp/withoutMinAge/abba-acsa_newJersey/summary/"
second_run = "../run/results/summary/"
output = "../run/results/mergedResults/summary/"

if (!dir.exists(first_run) | !dir.exists(second_run))
	stop("Directories do not exist!")

if (!dir.exists(output))
	dir.create(output, recursive = TRUE)

## List species
ls_species = list.dirs(path = first_run, recursive = FALSE, full.names = FALSE)
if (!all(list.dirs(path = second_run, recursive = FALSE, full.names = FALSE) %in% ls_species))
	warning("It seems that the list of species are not exactly the same between first and second run")

## Time data from second run
# Load parameters
simulationParameters = setDT(read.table(file = "../run/simulationParameters.txt", header = FALSE,
	sep = "=", comment.char = "#", blank.lines.skip = TRUE))
setnames(x = simulationParameters, new = c("parameters", "values"))

# Clean parameters
simulationParameters[, parameters := stringCleaner(parameters, " ")]
simulationParameters[, values := stringCleaner(values, " ")]
simulationParameters[, values := stringCleaner(values, fixed = "./", skip = "../"), by = parameters]

# Extract temporal informations
t0 = as.numeric(simulationParameters[parameters == "t0", values])
tmax = as.numeric(simulationParameters[parameters == "tmax", values])
nIter = as.integer(simulationParameters[parameters == "nIter", values])
delta_t = (tmax - t0)/(nIter - 1)

warning("I assume that delta_t is the same for both runs, and that the second run starts exactly where the first ended")

#### Merge and save results per patch
for (species in ls_species[2])
{
	## Load species-specific results
	current_path_1 = paste0(first_run, species, "/")
	current_path_2 = paste0(second_run, species, "/")
	current_output_path = paste0(output, species, "/")

	ls_files_1 = list.files(path = current_path_1, pattern = "^su_.*.txt$")
	ls_files_2 = list.files(path = current_path_2, pattern = "^su_.*.txt$")

	if (!all(ls_files_1 %in% ls_files_2))
		stop("Mismatch between number of files first versus second run")

	if (!dir.exists(current_output_path))
		dir.create(current_output_path)

	n_plots = length(ls_files_1)
	results_ls = vector(mode = "list", length = n_plots)
	count = 1
	for (current_file in ls_files_1)
	{
		results_1 = fread(paste0(current_path_1, current_file))
		results_1[, year := iteration*delta_t]

		results_2 = fread(paste0(current_path_2, current_file))
		results_2[, year := iteration*delta_t + results_1[.N, year]]

		if (results_2[.N, year] != tmax)
			stop("Time mismatch, check delta_t, tmax, nIter, etc...")

		if ((results_2[1, sumTrunkArea] != 0 & results_1[.N, sumTrunkArea]) &
			(results_1[.N, sumTrunkArea]/results_2[1, sumTrunkArea] < 0.98 | results_1[.N, sumTrunkArea]/results_2[1, sumTrunkArea] > 1.02))
			warning(paste("Mismatch trunk area larger than 2% for", current_file))

		if ((results_2[1, totalDensity] != 0 & results_1[.N, totalDensity]) &
			(results_1[.N, totalDensity]/results_2[1, totalDensity] < 0.95 | results_1[.N, totalDensity]/results_2[1, totalDensity] > 1.05))
			warning(paste("Mismatch total density larger than 5% for", current_file))

		results_2 = results_2[iteration != 0]
		results_2[, iteration := iteration + results_1[.N, iteration]]

		results_ls[[count]] = rbind(results_1, results_2)
		count = count + 1

		if (count %% 100 == 0)
			print(paste0(round(count*100/n_plots), "% done"))
	}
	results_ls = rbindlist(results_ls)
	results_ls[, filename := rep(ls_files_1, each = 2*nIter - 1)]
	results_ls[, patch_id := stri_sub(filename,
		from = stri_locate_first(filename, regex = "^su_")[, "end"] + 1,
		to = stri_locate_first(filename, regex = ".txt$")[, "start"] - 1), by = filename]
	results_ls[, filename := NULL]
	
	## Write species-specific results
	printRes(dt = results_ls, path = current_output_path, filenamePattern = "su_")
}
