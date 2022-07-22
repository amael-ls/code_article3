#### Aim of prog: Plot results analytical test
# Case 1: No dispersal, no competition.

#### Load packages and clear memory
library(data.table)
library(tikzDevice)
library(stringi)

rm(list = ls())
graphics.off()
options(max.print = 500)

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

#### Load results
pathResults = "./diffResults/"
if (!dir.exists(pathResults))
	stop("Directory is missing")

ls_files = list.files(pathResults)
nbFiles = length(ls_files)
nbFiles_density = sum(stri_detect(str = ls_files, regex = "density"))
nbFiles_dbh = sum(stri_detect(str = ls_files, regex = "dbh"))

#### Load parameters c++
## Simulation parameters
simulationParameters = setDT(read.table(file = "../run/simulationParameters.txt", header = FALSE,
	sep = "=", comment.char = "#", blank.lines.skip = TRUE))
setnames(x = simulationParameters, new = c("parameters", "values"))

## Clean parameters
simulationParameters[, parameters := stringCleaner(parameters, " ")]
simulationParameters[, values := stringCleaner(values, " ")]
simulationParameters[, values := stringCleaner(values, fixed = "./", skip = "../"), by = parameters]

tmax = as.numeric(simulationParameters[parameters == "tmax", values])

ls_results_density = vector(mode = "list", length = nbFiles_density)
ls_results_dbh = vector(mode = "list", length = nbFiles_dbh)
density_counter = 0
dbh_counter = 0

for (i in 1:nbFiles)
{
	file = ls_files[i]
	if (stri_detect(str = file, regex = "density"))
	{
		density_counter = density_counter + 1
		temp = readRDS(paste0(pathResults, file))
		if (stri_detect(str = file, regex = "init"))
			temp = data.table(iterationBirth = 0, delta_density = temp)

		steps = as.integer(stri_sub(str = file, from = stri_locate_last(str = file, regex = "_")[, "end"] + 1,
			to = stri_locate(str = file, regex = ".rds")[, "start"] - 1))

		temp[, nbSteps := steps]
		ls_results_density[[density_counter]] = temp
	}
		
	if (stri_detect(str = file, regex = "dbh"))
	{
		dbh_counter = dbh_counter + 1
		temp = readRDS(paste0(pathResults, file))
		if (stri_detect(str = file, regex = "init"))
			temp = data.table(iterationBirth = 0, delta_dbh = temp)

		steps = as.integer(stri_sub(str = file, from = stri_locate_last(str = file, regex = "_")[, "end"] + 1,
			to = stri_locate(str = file, regex = ".rds")[, "start"] - 1))

		temp[, nbSteps := steps]
		ls_results_dbh[[dbh_counter]] = temp
	}
}

ls_results_density = setorder(rbindlist(ls_results_density), nbSteps)
ls_results_dbh = setorder(rbindlist(ls_results_dbh), nbSteps)

ls_results_density[, deltaT := tmax/(nbSteps - 1)]
ls_results_density[, nbPoints := .N, by = nbSteps]

ls_results_dbh[, deltaT := tmax/(nbSteps - 1)]
ls_results_dbh[, nbPoints := .N, by = nbSteps]

#### Compute slopes
slopesZero = ls_results_density[iterationBirth == 0]
slopesZero[, log10_ddensity := log10(delta_density)]
slopesZero[, log10_deltaT := log10(deltaT)]

slopesZero[, slopes := (shift(log10_ddensity, 1, NA, "lead") - log10_ddensity)/(shift(log10_deltaT, 1, NA, "lead") - log10_deltaT)]

sl0 = mean(slopesZero[, slopes], na.rm = TRUE)

coeffs = coef(lm(slopesZero[, log10_ddensity] ~ slopesZero[, log10_deltaT]))
names(coeffs) = c("intercept", "slope")

#### Plot results
## Plot for density
tikz("error_density.tex", width = 3.1, height = 3.1)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(0, pch = "", xlim = slopesZero[, range(log10_deltaT)], ylim = slopesZero[, range(log10_ddensity)],
	xlab = "$\\log_{10}(\\Delta t) $", ylab = "$ \\log_{10}(Error density) $")
abline(coef = coeffs, col = "orange", lwd = 2)
points(slopesZero[, log10(deltaT)], slopesZero[, log10(delta_density)], pch = 19)
dev.off()

## Plot for dbh
tikz("error_dbh.tex", width = 3.1, height = 3.1)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(0, pch = "", xlim = ls_results_dbh[iterationBirth == 0, range(log10(deltaT))],
	ylim = ls_results_dbh[iterationBirth == 0, range(log10(delta_dbh))],
	xlab = "$\\log_{10}(\\Delta t) $", ylab = "$ \\log_{10}(Error dbh) $")
abline(
	coef = coef(lm(ls_results_dbh[iterationBirth == 0, log10(delta_dbh)] ~ ls_results_dbh[iterationBirth == 0, log10(deltaT)])),
	col = "orange", lwd = 2)
points(ls_results_dbh[iterationBirth == 0, log10(deltaT)], ls_results_dbh[iterationBirth == 0, log10(delta_dbh)], pch = 19)
dev.off()

print(paste0("The slope is: ", round(coeffs["slope"], 5)))
