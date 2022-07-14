
#### Aim of prog: compute yearly average for the four variables tmax, tmean, tmin, and prec
## *Variables and units:
# tmax			mean daily maximum air temperature					K/10
# tmean			mean daily air temperature							K/10
# tmin			mean daily minimum air temperature					K/10
# prec			precipitation amount								kg.m-2
#
# tas			near-surface (usually, 2 meter) air temperature		K/10
# tasmin		idem, but minimum									K/10
# tasmax		idem, but maximum									K/10
# pr			precipitation flux									kg.m-2/100
#
## *Comments
# Many important informations can be found on the official website of the data (especially the units!):
#	https://chelsa-climate.org
#
# Here is the pdf containing the informations for the version 2 (version read: 24 May 2021)
#	https://chelsa-climate.org/wp-admin/download-page/CHELSA_tech_specification_V2.pdf
#
# The three temperature variables are the daily maximum/average/minimum air temperatures at
#	2 metres since previous post-processing from ERA interim averaged over 1 month
#
# For the precipiations, "Amount" means mass per unit area while "Precipitation" in the earth's atmosphere
#	means precipitation of water in all phases.
#
## *Explanation parallel
# This script is parallelised on years. Therefore, each call will process all the variables for a specific year
#
# Due to the parallel programming extra care must be done on the creation of folders. If more than one core tries to create a same folder,
#	dir.create will throw a warning (but might overwrite some files). I added a security in the bash script (more explanation below).
#
# As you may noticed, the parallelisation is done outside of this program! Indeed, I use a for loop in a bash script. The variable array_id
# 	is inspired from SLURM, a workload manager on super computers. More informations at: https://slurm.schedmd.com/job_array.html. By
#	security,the first array is started one minute before the others, to make sure that all the folders are created (see problem that
#	could occur above).
#
## *Shapefile
#	The shapefile was created with the script create_usa-canada_shapefile.R

#### Clear memory and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(stringi)
library(terra)
