
#### Aim of prog: Create initial condition for analytical cases:
# Case 1: No dispersal, no fecundity, no competition. Path = ic_1/
# Case 2: No dispersal, no fecundity. Path = ic_2/
# Case 3: No dispersal, no competition. Path = ic_3/
# Case 4: No dispersal. Path = ic_4/

#### Load library and clear memory
library(data.table)

rm(list = ls())
graphics.off()

#### Tool function
printIC = function(densities, dbh, path, filenamePattern, id_plots = 1:nrow(densities), sep = " ", reset = TRUE)
{
	nbCohorts = nrow(densities)
	nbPlots = ncol(densities)
	if (length(dbh) != nbCohorts)
		stop(paste0("Number of cohorts (", nbCohorts, ") and dbh (", length(dbh), ") mismatch"))
	
	if (length(id_plots) != nbPlots)
		stop(paste0("Number of plots (", nbPlots, ") and id_plots (", length(id_plots), ") mismatch"))

	for (plot in 1:nbPlots)
	{
		outfileName = paste0(path, filenamePattern, id_plots[plot], ".txt")

		if (reset & file.exists(outfileName))
			file.remove(outfileName)

		ofstream = file(outfileName)

		line = paste0("density", sep, "dbh", sep = "\n")

		for (cohort in 1:nbCohorts)
			line = paste0(line, densities[cohort, plot], sep, dbh[cohort], sep = "\n")

		writeLines(line, ofstream)
		close(ofstream)
	}
}

#### Parameters
## Folder and id plots (for names initial condition)
filenamePattern = "ic_"

## Cohorts 
nbCohorts = 150
maxDiameter = 51
minDiameter = 2

#### Case number 1
## Parameters
path = "./ic_1/Acer_saccharum/"
id_plots = readRDS("../createLandscape/clim_5x5/populatedPatches.rds")
nbPlots = length(id_plots)

if (!dir.exists(path))
	dir.create(path)

## Generate cohorts
dbh = seq(minDiameter, maxDiameter, length.out = nbCohorts)
densities = matrix(data = 0, nrow = nbCohorts, ncol = nbPlots)

for (col in 1:nbPlots)
	densities[, col] = 1 + dbh

## Write files
printIC(densities, dbh, path, filenamePattern, id_plots, sep = " ", reset = TRUE)
