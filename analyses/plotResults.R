
#### Aim of prog: Plot and analyse the ouptuts of Cpp prog
## Plot densities (initial and last states)
#

#### Load packages
library(data.table)
library(tikzDevice)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

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
