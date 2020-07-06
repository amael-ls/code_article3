
#### Aim of prog: Plot and analyse the ouptuts of Cpp prog
## Plot densities (initial and last states)
#

#### Load packages
library(data.table)
library(tikzDevice)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

######## Part I: No dispersion effect (goudriann, 1986)
#### Load results c++
## Path
pathCpp = "../cpp/"

## Initial state
init = fread(paste0(pathCpp, "initNoDisp.txt"))

## Last state
end = fread(paste0(pathCpp, "endNoDisp.txt"))

## Competition, reproduction, basal area and total density
compReprod = fread(paste0(pathCpp, "compReprodNoDisp.txt"))

## Dynamics
# Load the density and dbh dynamics of cohorts
dyn = fread(paste0(pathCpp, "popDynNoDisp.txt"))

# Determine the number of cohorts and time steps
nbTimeSteps = dim(compReprod)[1]
nbCohorts = dim(dyn)[1]/nbTimeSteps - ifelse(end[.N, density] == 0, 1, 0)

#### Tikz plots
tikz("./noDisp.tex", width = 3.1, height = 3.1) #, standAlone = TRUE)
op <- par(mar = c(2.5, 2.5, 2.5, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
xmax = max(end[, dbh])
ymax = max(end[, density])
plot(x = NULL, y = NULL, xlim = c(0, xmax + 0.1),
	ylim = c(0, ymax + 0.1), axes = FALSE, xlab = "Size",
	ylab = "Density")
for (i in 1:nbCohorts)
{
	# Initial state
	tikzCoord(init[i, dbh], 0, paste0("init_start_", i))
	tikzCoord(init[i, dbh], init[i, density], paste0("init_end_", i))
	tikzAnnotate(paste0("\\draw (init_start_", i, ") -- (init_end_", i, ");"))

	# End state
	tikzCoord(end[i, dbh], 0, paste0("end_start_", i))
	tikzCoord(end[i, dbh], end[i, density], paste0("end_end_", i))
	tikzAnnotate(paste0("\\draw[dashed] (end_end_", i, ") -- (end_start_", i, ");"))

	# Arrow movement
	tikzAnnotate(paste0("\\draw[arrow, draw = orange] (init_end_", i, ") -- (end_end_", i, ");"))
}
# Axes
axis(side = 1, at = 0:xmax, labels = 0:xmax)
axis(side = 2, at = 0:ymax, labels = 0:ymax)

# Legend
legend(x = "topleft", legend = c("t = 0", "t = 3"), xpd = TRUE,
	lty = c("solid", "dashed"), lwd = 2, bty = "n", inset = c(0, -0.15))

dev.off()

######## Part II: Population dynamics
#### Load results c++
## Path
pathCpp = "../cpp/"

	## Initial state
init = fread(paste0(pathCpp, "init.txt"))

## Last state
end = fread(paste0(pathCpp, "end.txt"))

## Competition, reproduction, basal area and total density
compReprod = fread(paste0(pathCpp, "compReprod.txt"))

## Dynamics
# Load the density and dbh dynamics of cohorts
dyn = fread(paste0(pathCpp, "popDyn.txt"))

# Determine the number of cohorts and time steps
nbTimeSteps = dim(compReprod)[1]
nbCohorts = dim(dyn)[1]/nbTimeSteps

#### Plots
## Initial state
plot(init$dbh, init$density)

## Last state
plot(end$dbh, end$density)

## Dynamic plot
# Limits xlim and ylim
max_dbh = dyn[, max(dbh)]
max_density = dyn[, max(density)]

# Plot
plot(dyn[1:nbCohorts, dbh], dyn[1:nbCohorts, density], pch = '', # type = 'l',
	xlim = c(0, max_dbh), ylim = c(0, max_density))
invisible(sapply(2:nbTimeSteps, function(x, nbCohorts) {
	Sys.sleep(0.005)
	ind_start = (x - 1)*nbCohorts + 1
	ind_end = x*nbCohorts
	plot(dyn[ind_start:ind_end, dbh], dyn[ind_start:ind_end, density], type = 'l', lwd = 1,
		xlim = c(0, max_dbh), ylim = c(0, max_density))
}, nbCohorts = nbCohorts))

## Reproduction and competition
par(mfrow = c(4,1))
plot(compReprod$time, compReprod$reproduction, type = "l", lwd = 2,
	xlab = "Time", ylab = "Reproduction")
plot(compReprod$time, compReprod$competition, type = "l", lwd = 2,
	xlab = "Time", ylab = "Competition")
plot(compReprod$time, compReprod$basalArea, type = "l", lwd = 2,
	xlab = "Time", ylab = "Basal area")
plot(compReprod$time, compReprod$totalDensity, type = "l", lwd = 2,
	xlab = "Time", ylab = "Total density")

#### Tikz version
tikz(paste0("./ba_plot.tex"), width = 3.1, height = 3.1) #, standAlone = TRUE)
op <- par(mar = c(2.5, 2.5, 0.8, 0.8), mgp = c(1.5, 0.3, 0), tck = -0.015)
plot(compReprod$time, compReprod$basalArea, type = "l", lwd = 2,
	xlab = "Time", ylab = "Basal area")
dev.off()
