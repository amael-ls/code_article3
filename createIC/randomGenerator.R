
#### Aim of prog: Generate n random cohorts such that basal areas equal alpha
## Comments
# This program generates n random densities of cohort and their associated dbh, such that the basal area is
# alpha m^2/ha
# It is using a dirichlet distribution which generates n numbers that sum to 1.
# Multiply the output by alpha to get a sum to the desired BA
# For instance, to generate a basal area 25 m^2/ha = 2.5e-3 (dimensionless number)
# with 150 individuals we would do:
#	- Generate 150 dbh (in cm)
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

#### Load library and clear memory
library(rBeta2009) # If you cannot download it, MCMCpack also have a dirichlet random generator

rm(list = ls())
graphics.off()

#### Tool function
printIC = function(densities, dbh, path, filenamePattern, sep = " ", reset = TRUE)
{
	nbCohorts = ncol(densities)
	nbPlots = nrow(densities)
	if (length(dbh) != nbCohorts)
		stop(paste0("The number of cohorts (", nbCohorts, ") should be the same number of dbh provided (", length(dbh), ")"))

	for (plot in 1:nbPlots)
	{
		file = paste0(path, filenamePattern, plot, ".txt")

		if (reset & file.exists(file))
			file.remove(file)

		sink(file = file, append = TRUE)
		cat(paste0("density", sep, "dbh", sep = "\n"))

		for (cohort in 1:nbCohorts)
			cat(paste0(densities[plot, cohort], sep, dbh[cohort], sep = "\n"))

		sink(file = NULL)
	}
}

#### Generate cohorts
## Parameters
nbCohorts = 5
maxDiameter = 960
minDiameter = 2
nbPlots = 3

set.seed(1969-08-18) # Woodstock seed

alphaShape = abs(rnorm(nbCohorts, 0.5, 0.5))
BA = 25 # m^2/ha

## Generate beta numbers (read comments at the beginning)
betas = BA*rdirichlet(nbPlots, alphaShape)

## Generate dbh
dbh = sample(x = minDiameter:maxDiameter, size = nbCohorts, replace = TRUE)

## Deduce densities N
densities = matrix(data = 0, nrow = nrow(betas), ncol = ncol(betas))
for (i in 1:nbPlots)
	for(j in 1:nbCohorts)
		densities[i, j] = betas[i, j]*4*10^4/(pi*dbh[j]^2)

## Check cohorts basal area
for (i in 1:nbPlots)
	if (!all.equal((dbh^2 %*% t(densities)[,i]*pi/(4*10^4))[1,1], BA))
		print(paste0("Check row ", i, ". The basal area of this plot might not be the value expected"))

#### Write files initial condition
## Parameters
path = "./randomInitialCondition/"
filenamePattern = "ic_"

## Create folder
if (!dir.exists(path))
	dir.create(path)

## Write files
printIC(densities, dbh, path, filenamePattern, sep = " ", reset = TRUE)
