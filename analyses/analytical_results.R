
#### Aim of prog:  Create initial condition for analytical cases:
# Case 1: No dispersal, no fecundity, no competition. Path = ../code/popDyn/case1/
# Case 2: No dispersal, no fecundity. Path = ../code/popDyn/case2/
# Case 3: No dispersal, no competition. Path = ../code/popDyn/case3/
# Case 4: No dispersal. Path = ../code/popDyn/case4/

#### Load library and clear memory
library(data.table)

rm(list = ls())
graphics.off()

#### Case 1
## Parameters
# pathPopDyn = "../cpp/popDyn/case1/"
pathPopDyn = "../cpp/popDyn/Acer_saccharum/"
nIter = 500
tmax = 100

delta_t = tmax/(nIter - 1)

## Mortality
mu = 1/100

## Load results
results_case1 = fread(paste0(pathPopDyn, "pd_12.txt"))
results_case1[, time := iteration*delta_t]

results_case1[, range(iterationBirth)]

## Analytical solution
ode_II = function(N0, dbh0, mu, t)
	return(list(density = N0*exp(-mu*t), dbh = sqrt(4*t + (dbh0 + 1)^2) - 1))

N0_initCond = results_case1[iteration == 0, density]
dbh0_initCond = results_case1[iteration == 0, dbh]
# N0_birth_145 = results_case1[iterationBirth == 145, density]

## Comparison
# Initial condition (there should be small difference due to rounding error when writing ouptut)
max(abs(N0_initCond - results_case1[iteration == 0, density]))
max(abs(ode_II(N0_initCond, dbh0_initCond, mu, 0)[["density"]] - results_case1[iteration == 0, density]))

# Check analytical solution of the initial cohort
density_diff_init = -1
dbh_diff_init = -1

for (i in 1:499)
{
	delta_density = max(abs(ode_II(N0_initCond, dbh0_initCond, mu, i*delta_t)[["density"]] - results_case1[iteration == i & iterationBirth == 0, density]))
	if (density_diff_init < delta_density)
		density_diff_init = delta_density
	
	delta_dbh = max(abs(ode_II(N0_initCond, dbh0_initCond, mu, i*delta_t)[["dbh"]] - results_case1[iteration == i & iterationBirth == 0, dbh]))
	if (dbh_diff_init < delta_dbh)
		dbh_diff_init = delta_dbh
}

# Check for cohorts born later
ls_iterationBirth = results_case1[iterationBirth > 0, unique(iterationBirth)]
density_diff = data.table(iterationBirth = ls_iterationBirth, delta_density = -1)
dbh_diff = data.table(iterationBirth = ls_iterationBirth, delta_density = -1)

for (iterBirth in ls_iterationBirth)
{
	initDensity = results_case1[iteration == iterBirth & iterationBirth == iterBirth, density]
	initDBH = results_case1[iteration == iterBirth & iterationBirth == iterBirth, dbh]
	density_maxDiff = -1
	dbh_maxDiff = -1
	
	for (i in iterBirth:499)
	{	
		delta_density = max(abs(ode_II(initDensity, initDBH, mu, (i - iterBirth)*delta_t)[["density"]] - results_case1[iteration == i & iterationBirth == iterBirth, density]))
		if (density_maxDiff < delta_density)
			density_maxDiff = delta_density
		
		delta_dbh = max(abs(ode_II(initDensity, initDBH, mu, (i - iterBirth)*delta_t)[["dbh"]] - results_case1[iteration == i & iterationBirth == iterBirth, dbh]))
		if (dbh_maxDiff < delta_dbh)
			dbh_maxDiff = delta_dbh
	}
	density_diff[iterationBirth == iterBirth, delta_density := ..density_maxDiff]
	dbh_diff[iterationBirth == iterBirth, delta_density := ..dbh_maxDiff]
}

print(paste0("Max difference density: ", density_diff))
print(paste0("Max difference dbh: ", dbh_diff))

print(paste0("Max difference density: ", density_diff_init))
print(paste0("Max difference dbh: ", dbh_diff_init))
