
#### Aim of prog: Plot growth and mortality versus climate

#### Load packages and clear memory
rm(list = ls())
graphics.off()

library(data.table)
library(tikzDevice)
library(stringi)

options(max.print = 500)

#### Tool functions
## Demography
# Function to compute growth in function of a climate value x (either annual mean temperature or precipitation)
growth = function(x, otherClim, dbh, coeff_glm, is_x_temp = TRUE, scaling = NULL)
{
	scaling_G_mu = 0
	scaling_G_sd = 1
	if (is_x_temp)
	{
		params_G = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "growth", temp = x,
			precip = otherClim, isScaled = TRUE)
	} else {
		params_G = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "growth", temp = otherClim,
			precip = x, isScaled = TRUE)
	}

	if (!is.null(scaling))
	{
		if (!all(colnames(scaling) %in% c("parameters", "values")))
			stop("Scaling must contain 'parameters' and 'values'")
		
		setkey(scaling, parameters)
		scaling_G_mu = scaling["scaling_G_mu", values]
		scaling_G_sd = scaling["scaling_G_sd", values]
		scaling_dbh_mu_G = scaling["scaling_dbh_mu_G", values]
		scaling_dbh_sd_G = scaling["scaling_dbh_sd_G", values]

		dbh = (dbh - scaling_dbh_mu_G)/scaling_dbh_sd_G
		
		if (is_x_temp)
		{
			params_G = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "growth", temp = x,
				precip = otherClim, isScaled = FALSE, scaling = scaling)
		} else {
			params_G = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "growth", temp = otherClim,
				precip = x, isScaled = FALSE, scaling = scaling)
		}
	}

	c0 = params_G[stri_detect(names(params_G), regex = "^a0")]
	c1 = params_G[stri_detect(names(params_G), regex = "^a1")]
	c2 = params_G[stri_detect(names(params_G), regex = "^a2")]

	return (unname(exp(scaling_G_mu + scaling_G_sd*(c0 + c1*dbh + c2*dbh^2))))
}

# Function to compute mortality rate in function of a climate value x (either min annual temperature or precipitation driest quarter)
mortality = function(x, otherClim, dbh, coeff_glm, is_x_temp = TRUE, scaling = NULL)
{
	if (is_x_temp)
	{
		params_M = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "mortality", temp = x,
			precip = otherClim, isScaled = TRUE)
	} else {
		params_M = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "mortality", temp = otherClim,
			precip = x, isScaled = TRUE)
	}

	if (!is.null(scaling))
	{
		if (!all(colnames(scaling) %in% c("parameters", "values")))
			stop("Scaling must contain 'parameters' and 'values'")
		
		setkey(scaling, parameters)
		scaling_dbh_mu_M = scaling["scaling_dbh_mu_M", values]
		scaling_dbh_sd_M = scaling["scaling_dbh_sd_M", values]

		dbh = (dbh - scaling_dbh_mu_M)/scaling_dbh_sd_M

		if (is_x_temp)
		{
			params_M = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "mortality", temp = x,
				precip = otherClim, isScaled = FALSE, scaling = scaling)
		} else {
			params_M = prepareConstants(dt = coeff_glm, layer = "overstorey", demography = "mortality", temp = otherClim,
				precip = x, isScaled = FALSE, scaling = scaling)
		}
	}

	d0 = params_M[stri_detect(names(params_M), regex = "^a0")]
	d1 = rep(params_M["a1"], length(d0))
	d2 = rep(params_M["a2"], length(d0))

	return (unname(exp(d0 + d1*dbh + d2*dbh^2)))
}

## Prepare constants
prepareConstants = function(dt, layer, demography, temp, precip, isScaled = FALSE, ...)
{
	providedArgs = list(...)
	ls_names = names(providedArgs)

	if (!(demography %in% c("growth", "mortality")))
		stop("Demography must be either 'growth' or 'mortality'")

	if (!all(colnames(dt) %in% c("parameters", "values")))
		stop("dt must contains the columns 'parameters' and 'values'")
	
	if (!(layer %in% c("understorey", "overstorey")))
		stop("Layer must be either 'understorey' or 'overstorey'")

	if (!isScaled)
	{
		if (!("scaling" %in% ls_names))
			stop("You must provide the scaling data")
		scaling = providedArgs[["scaling"]]
		setkey(scaling, parameters)
		
		if (demography == "growth")
		{
			temp_mu = scaling["scaling_temp_mu_G", values]
			temp_sd = scaling["scaling_temp_sd_G", values]
			precip_mu = scaling["scaling_precip_mu_G", values]
			precip_sd = scaling["scaling_precip_sd_G", values]
		} else {
			temp_mu = scaling["scaling_temp_mu_M", values]
			temp_sd = scaling["scaling_temp_sd_M", values]
			precip_mu = scaling["scaling_precip_mu_M", values]
			precip_sd = scaling["scaling_precip_sd_M", values]
		}

		temp = (temp - temp_mu)/temp_sd
		precip = (precip - precip_mu)/precip_sd
	}

	setkey(dt, parameters)
	isOverstorey = layer == "overstorey"

	canopyFactor = 0
	if (isOverstorey)
		canopyFactor = 1
	
	if (demography == "growth")
	{
		a0 = dt["intercept", values] + canopyFactor*dt["cs", values] +
			(dt["T", values] + canopyFactor*dt["cs_T", values])*temp + (dt["T_sq", values] + canopyFactor*dt["cs_T_sq", values])*temp^2 +
			(dt["P", values] + canopyFactor*dt["cs_P", values])*precip + (dt["P_sq", values] + canopyFactor*dt["cs_P_sq", values])*precip^2

		a1 = dt["dbh", values] + dt["dbh_T", values]*temp + dt["dbh_T_sq", values]*temp^2 +
		dt["dbh_P", values]*precip + dt["dbh_P_sq", values]*precip^2

		a2 = dt["dbh_sq", values] + dt["dbh_sq_T", values]*temp + dt["dbh_sq_T_sq", values]*temp^2 +
		dt["dbh_sq_P", values]*precip + dt["dbh_sq_P_sq", values]*precip^2
	}

	if (demography == "mortality")
	{
		a0 = dt["intercept", values] + canopyFactor*(dt["cs", values] + dt["cs_T", values]*temp + dt["cs_T_sq", values]*temp^2 +
			dt["cs_P", values]*precip + dt["cs_P_sq", values]*precip^2) +
			dt["T", values]*temp + dt["T_sq", values]*temp^2 +dt["P", values]*precip + dt["P_sq", values]*precip^2

		a1 = dt["dbh", values]

		a2 = dt["dbh_sq", values]
	}

	if (any(is.na(c(a0, a1, a2))))
		warning("NA detected, check that dt contains all the required parameters")

	return(c(a0 = a0, a1 = a1, a2 = a2))
}

#### Read parameters
## Abies balsamea, growth (G) and mortality (M)
abba_G = fread("../run/data/speciesParameters/Abies_balsamea_G.txt", sep = "=", header = FALSE)
abba_M = fread("../run/data/speciesParameters/Abies_balsamea_M.txt", sep = "=", header = FALSE)

scaling_abba_G = fread("../run/data/speciesParameters/Abies_balsamea_scaling_G.txt", sep = "=", header = FALSE)
scaling_abba_M = fread("../run/data/speciesParameters/Abies_balsamea_scaling_M.txt", sep = "=", header = FALSE)

setnames(abba_G, new = c("parameters", "values"))
setnames(abba_M, new = c("parameters", "values"))
setnames(scaling_abba_G, new = c("parameters", "values"))
setnames(scaling_abba_M, new = c("parameters", "values"))

## Acer saccharrum, growth (G) and mortality (M)
acsa_G = fread("../run/data/speciesParameters/Acer_saccharum_G.txt", sep = "=", header = FALSE)
acsa_M = fread("../run/data/speciesParameters/Acer_saccharum_M.txt", sep = "=", header = FALSE)

scaling_acsa_G = fread("../run/data/speciesParameters/Acer_saccharum_scaling_G.txt", sep = "=", header = FALSE)
scaling_acsa_M = fread("../run/data/speciesParameters/Acer_saccharum_scaling_M.txt", sep = "=", header = FALSE)

setnames(acsa_G, new = c("parameters", "values"))
setnames(acsa_M, new = c("parameters", "values"))
setnames(scaling_acsa_G, new = c("parameters", "values"))
setnames(scaling_acsa_M, new = c("parameters", "values"))

#### Read climate values
## Climate in Orford
climate = fread("../run/data/landscape_300x11_abba-acsa_orford/climate_0.txt", sep = "=", header = FALSE) # All climate files are the same
setnames(climate, new = c("variables", "values"))
setkey(climate, variables)
climate = climate[stri_detect(str = variables, regex = c("precipitation")) | stri_detect(str = variables, regex = c("temperature")),]
climate[, values := as.numeric(values)]

temperature_G = climate["annual_mean_temperature", values]
precipitations_G = climate["annual_precipitation", values]

temperature_M = climate["min_temperature_of_coldest_month", values]
precipitations_M = climate["precipitation_of_driest_quarter", values]

#### Plots
## Growth vs ...
# ... temperature
tikz("growth_vs_temperature_newJersey.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.2, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(growth(x, otherClim = precipitations_G, dbh = 500, coeff_glm = abba_G, is_x_temp = TRUE, scaling = scaling_abba_G), from = -5,
	to = 20, ylim = c(0, 2.6), xlab = "Temperature (in Celsius)", ylab = "Growth (in mm/yr)", lwd = 3, col = "#1D4965")
curve(growth(x, otherClim = precipitations_G, dbh = 500, coeff_glm = acsa_G, is_x_temp = TRUE, scaling = scaling_acsa_G), add = TRUE,
	lwd = 3, col = "#2F90A8")

axis(side = 1, at = temperature_G, labels = "", las = 1, col.ticks = "#E5494B", lwd.ticks = 2)

points(x = round(temperature_G, digits = 1), y = growth(round(temperature_G, digits = 1), otherClim = precipitations_G, dbh = 500,
	coeff_glm = abba_G, is_x_temp = TRUE, scaling = scaling_abba_G), col = "#E5494B", pch = 20, cex = 1.4)
points(x = round(temperature_G, digits = 1), y = growth(round(temperature_G, digits = 1), otherClim = precipitations_G, dbh = 500,
	coeff_glm = acsa_G, is_x_temp = TRUE, scaling = scaling_acsa_G), col = "#E5494B", pch = 20, cex = 1.4)

dev.off()

# ... precipitations
tikz("growth_vs_precipitations_newJersey.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.2, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(growth(x, otherClim = temperature_G, dbh = 500, coeff_glm = abba_G, is_x_temp = FALSE, scaling = scaling_abba_G), from = 500,
	to = 1800, ylim = c(0, 3.2), xlab = "Precipitations (in mm)", ylab = "Growth (in mm/yr)", lwd = 3, col = "#1D4965")
curve(growth(x, otherClim = temperature_G, dbh = 500, coeff_glm = acsa_G, is_x_temp = FALSE, scaling = scaling_acsa_G), add = TRUE,
	lwd = 3, col = "#2F90A8")

axis(side = 1, at = precipitations_G, labels = "", las = 1, col.ticks = "#E5494B", lwd.ticks = 2)

points(x = round(precipitations_G, digits = 1), y = growth(round(precipitations_G, digits = 1), otherClim = temperature_G, dbh = 500,
	coeff_glm = abba_G, is_x_temp = FALSE, scaling = scaling_abba_G), col = "#E5494B", pch = 20, cex = 1.4)
points(x = round(precipitations_G, digits = 1), y = growth(round(precipitations_G, digits = 1), otherClim = temperature_G, dbh = 500,
	coeff_glm = acsa_G, is_x_temp = FALSE, scaling = scaling_acsa_G), col = "#E5494B", pch = 20, cex = 1.4)

dev.off()

## Mortality vs ...
# ... temperature
tikz("mortality_vs_temperature_newJersey.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.2, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(mortality(x, otherClim = precipitations_M, dbh = 500, coeff_glm = abba_M, is_x_temp = TRUE, scaling = scaling_abba_M), from = -10,
	to = 10, ylim = c(0, 0.012), xlab = "Temperature (in Celsius)", ylab = "Mortality", lwd = 3, col = "#1D4965")
curve(mortality(x, otherClim = precipitations_M, dbh = 500, coeff_glm = acsa_M, is_x_temp = TRUE, scaling = scaling_acsa_M), add = TRUE,
	lwd = 3, col = "#2F90A8")

axis(side = 1, at = temperature_M, labels = "", las = 1, col.ticks = "#E5494B", lwd.ticks = 2)

points(x = round(temperature_M, digits = 1), y = mortality(round(temperature_M, digits = 1), otherClim = precipitations_M, dbh = 500,
	coeff_glm = abba_M, is_x_temp = TRUE, scaling = scaling_abba_M), col = "#E5494B", pch = 20, cex = 1.4)
points(x = round(temperature_M, digits = 1), y = mortality(round(temperature_M, digits = 1), otherClim = precipitations_M, dbh = 500,
	coeff_glm = acsa_M, is_x_temp = TRUE, scaling = scaling_acsa_M), col = "#E5494B", pch = 20, cex = 1.4)

dev.off()

# ... precipitations
tikz("mortality_vs_precipitations_newJersey.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.2, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(mortality(x, otherClim = temperature_M, dbh = 500, coeff_glm = abba_M, is_x_temp = FALSE, scaling = scaling_abba_M), from = 150,
	to = 300, ylim = c(0, 0.017), xlab = "Precipitations (in mm)", ylab = "Mortality", lwd = 3, col = "#1D4965")
curve(mortality(x, otherClim = temperature_M, dbh = 500, coeff_glm = acsa_M, is_x_temp = FALSE, scaling = scaling_acsa_M), add = TRUE,
	lwd = 3, col = "#2F90A8")

axis(side = 1, at = precipitations_M, labels = "", las = 1, col.ticks = "#E5494B", lwd.ticks = 2)

points(x = round(precipitations_M, digits = 1), y = mortality(round(precipitations_M, digits = 1), otherClim = temperature_M, dbh = 500,
	coeff_glm = abba_M, is_x_temp = FALSE, scaling = scaling_abba_M), col = "#E5494B", pch = 20, cex = 1.4)
points(x = round(precipitations_M, digits = 1), y = mortality(round(precipitations_M, digits = 1), otherClim = temperature_M, dbh = 500,
	coeff_glm = acsa_M, is_x_temp = FALSE, scaling = scaling_acsa_M), col = "#E5494B", pch = 20, cex = 1.4)

dev.off()
