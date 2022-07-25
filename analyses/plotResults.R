
#### Aim of prog: Plot results (initial and last densities) for the appendix ???
## Reference:
# See appendix ??? in the paper:
# Development of a spatial Escalator Boxcar Train algorithm for sessile species: application to the boreal-temperate forest ecotone
# Link: 

#### Load library and clear memory
library(data.table)
library(tikzDevice)
library(stringi)

rm(list = ls())
graphics.off()

#### Options tikzDevice
options(tikzDefaultEngine = "luatex")

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
pathPopDyn = paste0(pathCpp, simulationParameters[parameters == "popDynFilePath", values], species, "/")

## Time
t0 = as.numeric(simulationParameters[parameters == "t0", values])
tmax = as.numeric(simulationParameters[parameters == "tmax", values])
(nIter = as.integer(simulationParameters[parameters == "nIter", values]))
delta_t = (tmax - t0)/(nIter - 1)

#### Results
## Loading
results = fread(paste0(pathPopDyn, "pd_0.txt"))

## Add column and keep only row/columns of interest
results[, time := iteration*delta_t]
results = results[iteration %in% c(0, nIter - 1)]
results[, c("iterationBirth", "height") := NULL]

## Initial condition
N0_initCond = results[iteration == 0, density]
dbh0_initCond = results[iteration == 0, dbh]

## Last state
N0_end = results[iteration == nIter - 1, density]
dbh0_end = results[iteration == nIter - 1, dbh]

#### Plot
tikz("noDisp.tex", width = 3.1, height = 3.1,
	packages = c(getOption("tikzLatexPackages"), "\\usetikzlibrary{arrows}"))
op = par(mar = c(4, 4, 3, 0.8), mgp = c(2, 0.9, 0), tck = -0.015, xpd = TRUE)
plot(0, type = "n", axes = FALSE, xlim = c(0, max(dbh0_end)), ylim = c(0, max(N0_end)), ann = FALSE)
axis(side = 1, at = 0:max(dbh0_end))
axis(side = 2, at = 0:max(N0_end), las = 1)
title(xlab = "Size", ylab = "Density")
for (i in 1:length(N0_end))
{
	# Densities
	segments(x0 = dbh0_initCond[i], y0 = 0, x1 = dbh0_initCond[i], y1 = N0_initCond[i], lty = "solid", lwd = 1)
	segments(x0 = dbh0_end[i], y0 = 0, x1 = dbh0_end[i], y1 = N0_end[i], lty = "dashed", lwd = 1)
	
	# Arrows from initial condition to last state
	tikzCoord(dbh0_initCond[i], N0_initCond[i], paste0("arrow_start_", i))
	tikzCoord(dbh0_end[i], N0_end[i], paste0("arrow_end_", i))
	tikzAnnotate(paste0("\\draw[arrow, draw = orange] (arrow_start_", i,") -- (arrow_end_", i,");"))
}
legend(x = "topleft", inset = c(0, -0.2), legend = c("t = 0", "t = 3"), lty = c("solid", "dashed"), bty = "n")
dev.off()
