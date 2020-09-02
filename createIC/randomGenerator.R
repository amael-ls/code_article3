
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
## Remark
# If "id_plots" is requiered, then it is necessary to first create the landscape (see ../createLandscape/)
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
		file = paste0(path, filenamePattern, id_plots[plot], ".txt")

		if (reset & file.exists(file))
			file.remove(file)

		sink(file = file, append = TRUE)
		cat(paste0("density", sep, "dbh", sep = "\n"))

		for (cohort in 1:nbCohorts)
			cat(paste0(densities[plot, cohort], sep, dbh[cohort], sep = "\n"))

		sink(file = NULL)
	}
}

#### Parameters
## Folder and id plots (for names initial condition)
path = "./randomInitialCondition/"
filenamePattern = "ic_"
id_plots = readRDS("../createLandscape/populatedPatches.rds")

## Cohorts 
nbCohorts = 5
maxDiameter = 960
minDiameter = 2
nbPlots = length(id_plots)

#### Generate cohorts
set.seed(1969-08-18) # Woodstock seed

alphaShape = abs(rnorm(nbCohorts, 0.5, 0.5))
BA = 25 # m^2/ha

## Generate beta numbers (read comments at the beginning)
betas = BA*rdirichlet(nbPlots, alphaShape)

## Generate dbh (in cm)
dbh = sample(x = minDiameter:maxDiameter, size = nbCohorts, replace = TRUE)

dbh_meters = dbh/100

## Deduce densities N
densities = matrix(data = 0, nrow = nrow(betas), ncol = ncol(betas))
for (i in 1:nbPlots)
	for(j in 1:nbCohorts)
		densities[i, j] = betas[i, j]*4*10^4/(pi*dbh_meters[j]^2)

## Check cohorts basal area
for (i in 1:nbPlots)
	if (!all.equal((dbh_meters^2 %*% t(densities)[,i]*pi/(4*10^4))[1,1], BA))
		print(paste0("Check row ", i, ". The basal area of this plot might not be the value expected"))

## Number of trees per patch
data.table(id_plots, density_ha = rowSums(densities))

#### Write files initial condition
## Create folder
if (!dir.exists(path))
	dir.create(path)

## Write files
printIC(densities, dbh, path, filenamePattern, id_plots, sep = " ", reset = TRUE)
