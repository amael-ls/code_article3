
#### Aim of prog: Compare results with analytical solution
## Reference:
# For the analytical solution, see the development in the paper
# Development of a spatial Escalator Boxcar Train algorithm for sessile species: application to the boreal-temperate forest ecotone
# Link: 

#### Load library and clear memory
library(data.table)
library(stringi)

rm(list = ls())
graphics.off()

#### Tool function
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
species = simulationParameters[parameters == "species_filenames", values]
species = stringCleaner(species, ".txt")

## Paths to files
pathCpp = "../"

pathSummary = paste0(pathCpp, simulationParameters[parameters == "summaryFilePath", values], species, "/")
pathPopDyn = paste0(pathCpp, simulationParameters[parameters == "popDynFilePath", values], species, "/")
initPath = paste0(pathCpp, simulationParameters[parameters == "initPath", values], species, "/")

## Files' patterns
summaryPattern = simulationParameters[parameters == "summaryFilePattern", values]
initPattern = simulationParameters[parameters == "initFilenamePattern", values]

## Time
t0 = as.numeric(simulationParameters[parameters == "t0", values])
tmax = as.numeric(simulationParameters[parameters == "tmax", values])
(nIter = as.integer(simulationParameters[parameters == "nIter", values]))
delta_t = (tmax - t0)/(nIter - 1)
print(log10(delta_t))

## Mortality (see reference at the beginning of this script)
mu = 1/100

## Load results
results = fread(paste0(pathPopDyn, "pd_22.txt"))
results[, time := iteration*delta_t]

results[, unique(iterationBirth)]

## Analytical solution (see reference at the beginning of this script)
ode_II = function(N0, dbh0, mu, t)
	return(list(density = N0*exp(-mu*t), dbh = sqrt(4*t + (dbh0 + 1)^2) - 1))

N0_initCond = results[iteration == 0, density]
dbh0_initCond = results[iteration == 0, dbh]

## Comparison
# Initial condition (there might be small difference due to rounding error when writing output)
max(abs(N0_initCond - results[iteration == 0, density]))
max(abs(ode_II(N0_initCond, dbh0_initCond, mu, 0)[["density"]] - results[iteration == 0, density]))

# Check analytical solution of the initial cohort
density_diff_init = -Inf
dbh_diff_init = -Inf

for (i in 1:(nIter - 1))
{
	delta_density = max(abs(ode_II(N0_initCond, dbh0_initCond, mu, i*delta_t)[["density"]] - results[iteration == i & iterationBirth == 0, density]))
	if (density_diff_init < delta_density)
		density_diff_init = delta_density
	
	delta_dbh = max(abs(ode_II(N0_initCond, dbh0_initCond, mu, i*delta_t)[["dbh"]] - results[iteration == i & iterationBirth == 0, dbh]))
	if (dbh_diff_init < delta_dbh)
		dbh_diff_init = delta_dbh
}

# Check for cohorts born later
ls_iterationBirth = results[iterationBirth > 0, unique(iterationBirth)]
density_diff = data.table(iterationBirth = ls_iterationBirth, delta_density = -Inf)
dbh_diff = data.table(iterationBirth = ls_iterationBirth, delta_dbh = -Inf)

for (iterBirth in ls_iterationBirth)
{
	initDensity = results[iteration == iterBirth & iterationBirth == iterBirth, density]
	initDBH = results[iteration == iterBirth & iterationBirth == iterBirth, dbh]
	density_maxDiff = -1
	dbh_maxDiff = -1
	
	for (i in iterBirth:(nIter - 1))
	{	
		delta_density = max(abs(ode_II(initDensity, initDBH, mu, (i - iterBirth)*delta_t)[["density"]] - results[iteration == i & iterationBirth == iterBirth, density]))
		if (density_maxDiff < delta_density)
			density_maxDiff = delta_density
		
		delta_dbh = max(abs(ode_II(initDensity, initDBH, mu, (i - iterBirth)*delta_t)[["dbh"]] - results[iteration == i & iterationBirth == iterBirth, dbh]))
		if (dbh_maxDiff < delta_dbh)
			dbh_maxDiff = delta_dbh
	}
	density_diff[iterationBirth == iterBirth, delta_density := ..density_maxDiff]
	dbh_diff[iterationBirth == iterBirth, delta_dbh := ..dbh_maxDiff]
}

print(density_diff)
print(dbh_diff)

print(paste0("Max difference density init: ", density_diff_init))
print(paste0("Max difference dbh init: ", dbh_diff_init))

if (!dir.exists("./diffResults"))
	dir.create("./diffResults")

saveRDS(density_diff, paste0("./diffResults/density_", nIter, ".rds"))
saveRDS(dbh_diff, paste0("./diffResults/dbh_", nIter, ".rds"))
saveRDS(density_diff_init, paste0("./diffResults/density_init_", nIter, ".rds"))
saveRDS(dbh_diff_init, paste0("./diffResults/dbh_init_", nIter, ".rds"))
