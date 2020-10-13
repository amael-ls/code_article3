
#### Aim of prog: Generate n random cohorts such that basal areas equal alpha
## Comments
# This program generates n random densities of cohort and their associated dbh, such that the basal area is
# alpha m^2/ha
# It is using a dirichlet distribution which generates n numbers that sum to 1.
# Multiply the output by alpha to get a sum to the desired BA
# For instance, to generate a basal area 25 m^2/ha = 2.5e-3 (dimensionless number)
# with 150 individuals we would do:
#	- Generate 150 dbh (if in cm, convert to meters)
#	- Generate 150 numbers β using dirichlet distribution: β_i = N_i * π * dbh_i^2/(4 . 10^4) such that the sum is alpha
#	- From β_i, deduce N_i
#	- Check that Σ N_i * π * dbh_i^2/(4 . 10^4) = alpha
#
## Output:
#	- The output is a matrix of size nbPlots x nbCohorts
#	- Each line is for one plot
#	- Each column is associated to a dbh (i.e., same size-structure among plots)
#
# See https://reference.wolfram.com/language/ref/DirichletDistribution.html for more explanation on this distribtion
#
## Remarks
# 1/. If "id_plots" is requiered, then it is necessary to first create the landscape (see ../createLandscape/)
# 2/. To give you an idea at different scales:
# 		- 1 km^2 might contain more than 60,000 trees (https://www.reddit.com/r/MapPorn/comments/47f9s5/trees_per_square_km_in_europe_1200x1000/)
#		- 1 ha might contain up to 2,500 trees (plantation), between 620-1,500 (regeneration) or 15-150 (old-growth forest)
# 			source: https://naldc.nal.usda.gov/download/34466/PDF
#

#### Load library and clear memory
library(data.table)
library(rBeta2009) # If you cannot download it, MCMCpack also have a dirichlet random generator

rm(list = ls())
graphics.off()

#### Tool function
printIC = function(densities, dbh, path, filenamePattern, id_plots = 1:nrow(densities), sep = " ", reset = TRUE)
{
	nbCohorts = ncol(densities)
	nbPlots = nrow(densities)
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
			line = paste0(line, densities[plot, cohort], sep, dbh[cohort], sep = "\n")

		writeLines(line, ofstream)
		close(ofstream)
	}
}

#### Parameters
## Folder and id plots (for names initial condition)
path = "./randomInitialCondition/Acer_saccharum/"
filenamePattern = "ic_"
id_plots = readRDS("../createLandscape/populatedPatches.rds")

## Cohorts 
nbCohorts = 150
maxDiameter = 960
minDiameter = 2
nbPlots = length(id_plots)

#### Generate cohorts
set.seed(1969-08-18) # Woodstock seed

alphaShape = abs(rnorm(nbCohorts, 0.5, 0.5))
BA = 25 # m^2/ha
areaPlot = 20*20 # 1ha = 100 x 100 m^2
targetedArea = areaPlot*BA/10000

## Generate beta numbers (read comments at the beginning)
betas = targetedArea*rdirichlet(nbPlots, alphaShape)

## Generate dbh (in mm)
dbh = sample(x = minDiameter:maxDiameter, size = nbCohorts, replace = TRUE)

dbh_meters = dbh/1000

## Deduce densities N
densities = matrix(data = 0, nrow = nrow(betas), ncol = ncol(betas))
for (i in 1:nbPlots)
	for(j in 1:nbCohorts)
		densities[i, j] = betas[i, j]*4/(pi*dbh_meters[j]^2)

## Check cohorts basal area
for (i in 1:nbPlots)
	if (!all.equal((10^4*dbh_meters^2 %*% t(densities)[,i]*pi/(4*areaPlot))[1,1], BA))
		print(paste0("Check row ", i, ". The basal area of this plot might not be the value expected"))

## Number of trees per patch
data.table(id_plots, density_patch = rowSums(densities))

## Remove all values smaller than 1e-50
# I do this because std::stod in C++ cannot handle smaller numbers than std::numeric_limits<double>::min() which is around 1e-308 in my computer
densities[densities < 1e-50] = 0

## Check cohorts basal area again
for (i in 1:nbPlots)
	if (!all.equal((10^4*dbh_meters^2 %*% t(densities)[,i]*pi/(4*areaPlot))[1,1], BA))
		print(paste0("Check row ", i, ". The basal area of this plot might not be the value expected"))

#### Write files initial condition
## Create folder
if (!dir.exists(path))
	dir.create(path)

## Write files
printIC(densities, dbh, path, filenamePattern, id_plots, sep = " ", reset = TRUE)
