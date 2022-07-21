
#### Aim of prog: Create initial condition for analytical case

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
	
	if (id_plots[, .N] != nbPlots)
		stop(paste0("Number of plots (", nbPlots, ") and id_plots (", id_plots[, .N], ") mismatch"))

	for (plot in 1:nbPlots)
	{
		outfileName = paste0(path, filenamePattern, id_plots[plot, patch_id], ".txt")
	
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
path = "../run/data/initialCondition/Acer_saccharum/"
id_plots = readRDS("../run/data/landscape_5x5_acsa/populatedPatches.rds")
nbPlots = id_plots[, .N]

if (!dir.exists(path))
	dir.create(path)

## Generate cohorts
dbh = seq(minDiameter, maxDiameter, length.out = nbCohorts)
densities = matrix(data = 0, nrow = nbCohorts, ncol = nbPlots)

for (col in 1:nbPlots)
	densities[, col] = 1 + dbh

## Write files
printIC(densities, dbh, path, filenamePattern, id_plots, sep = " ", reset = TRUE)
