
#### Aim of prog: Generate the initial condition from a previous run
#

#### Load library and clear memory
library(data.table)
library(stringi)

rm(list = ls())
graphics.off()

#### Tool function
printIC = function(density_vec, dbh_vec, patch_id, filenamePattern, outputPath, sep = " ", rm = FALSE)
{
	patch_id = unique(patch_id)
	nbCohorts = length(density_vec)

	outfileName = paste0(outputPath, filenamePattern, patch_id, ".txt")
	if (rm & file.exists(outfileName))
		file.remove(outfileName)

	ofstream = file(outfileName)
	line = paste0("density", sep, "dbh", sep = "\n")

	for (cohort in 1:nbCohorts)
		line = paste0(line, density_vec[cohort], sep, dbh_vec[cohort], sep = "\n")

	writeLines(line, ofstream)
	close(ofstream)
}

#### Parameters
## Folder and id plots (for names initial condition)
path_ic = "./randomInitialCondition/Acer_saccharum/"
filenamePattern = "ic_"
id_plots = readRDS("../createLandscape/populatedPatches.rds")

## Get results previous run
# Path
path_data = "../cpp/popDyn/Acer_saccharum/"

# Pattern population dnamics
loadPattern = "pd_"

# List population dynamics file
files_ls = list.files(path = path_data, pattern = paste0("^", loadPattern, ".*.txt"))
nbFiles = length(files_ls)

# Load files
popDyn_ls = vector(mode = "list", length = nbFiles)

for (i in 1:nbFiles)
{
	filename = paste0(path_data, files_ls[i])
	popDyn_ls[[i]] = fread(filename)[iteration == max(iteration)]
}

popDyn_ls = rbindlist(popDyn_ls, idcol = "patch_id")

# Filter low densities because std::stod in C++ cannot handle smaller numbers than std::numeric_limits<double>::min()
popDyn_ls = popDyn_ls[density > 1e-50]

#### Few plots
plot(popDyn_ls[patch_id == 1, dbh], popDyn_ls[patch_id == 1, density], type = "l", lwd = 2,
	xlab = "dbh", ylab = "density")
lines(popDyn_ls[patch_id == 200, dbh], popDyn_ls[patch_id == 200, density], lwd = 2, col = "blue")

#### Write files initial condition
## Create folder
if (!dir.exists(path_ic))
	dir.create(path_ic)

## Write files
popDyn_ls[, printIC(density, dbh, patch_id, filenamePattern, path_ic, sep = " ", rm = TRUE), by = patch_id]
