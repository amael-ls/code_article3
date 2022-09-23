
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
# Function to compute growth in function of dbh
growth = function(dbh, temp, precip, coeff_glm, layer = "overstorey", scaling = NULL)
{
	scaling_G_mu = 0
	scaling_G_sd = 1

	params_G = prepareConstants(dt = coeff_glm, layer = layer, demography = "growth", temp = temp, precip = precip, isScaled = TRUE)

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

		params_G = prepareConstants(dt = coeff_glm, layer = layer, demography = "growth", temp = temp, precip = precip,
			isScaled = FALSE, scaling = scaling)
	}

	c0 = params_G[stri_detect(names(params_G), regex = "^a0")]
	c1 = params_G[stri_detect(names(params_G), regex = "^a1")]
	c2 = params_G[stri_detect(names(params_G), regex = "^a2")]

	return (unname(exp(scaling_G_mu + scaling_G_sd*(c0 + c1*dbh + c2*dbh^2))))
}

# Function to compute mortality rate in function of dbh
mortality = function(dbh, temp, precip, coeff_glm, layer = "overstorey", scaling = NULL)
{
	params_M = prepareConstants(dt = coeff_glm, layer = layer, demography = "mortality", temp = temp, precip = precip, isScaled = TRUE)

	if (!is.null(scaling))
	{
		if (!all(colnames(scaling) %in% c("parameters", "values")))
			stop("Scaling must contain 'parameters' and 'values'")
		
		setkey(scaling, parameters)
		scaling_dbh_mu_M = scaling["scaling_dbh_mu_M", values]
		scaling_dbh_sd_M = scaling["scaling_dbh_sd_M", values]

		dbh = (dbh - scaling_dbh_mu_M)/scaling_dbh_sd_M

		params_M = prepareConstants(dt = coeff_glm, layer = layer, demography = "mortality", temp = temp, precip = precip,
			isScaled = FALSE, scaling = scaling)
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
orford = fread("../run/data/landscape_300x11_abba-acsa_orford/climate_0.txt", sep = "=", header = FALSE) # All climate files are the same
setnames(orford, new = c("variables", "values"))
setkey(orford, variables)
orford = orford[stri_detect(str = variables, regex = c("precipitation")) | stri_detect(str = variables, regex = c("temperature")),]
orford[, values := as.numeric(values)]

temperature_G_orford = orford["annual_mean_temperature", values]
precipitations_G_orford = orford["annual_precipitation", values]

temperature_M_orford = orford["min_temperature_of_coldest_month", values]
precipitations_M_orford = orford["precipitation_of_driest_quarter", values]

## Climate in New-Jersey
newJersey = fread("../run/data/landscape_300x11_abba-acsa_newJersey/climate_0.txt", sep = "=", header = FALSE)
setnames(newJersey, new = c("variables", "values"))
setkey(newJersey, variables)
newJersey = newJersey[stri_detect(str = variables, regex = c("precipitation")) | stri_detect(str = variables, regex = c("temperature")),]
newJersey[, values := as.numeric(values)]

temperature_G_newJersey = newJersey["annual_mean_temperature", values]
precipitations_G_newJersey = newJersey["annual_precipitation", values]

temperature_M_newJersey = newJersey["min_temperature_of_coldest_month", values]
precipitations_M_newJersey = newJersey["precipitation_of_driest_quarter", values]

#### Plots
## Growth vs ...
# ... dbh in Orford
tikz("growth_vs_dbh_orford.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.8, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(growth(x, temp = temperature_G_orford, precip = precipitations_G_orford, coeff_glm = abba_G, scaling = scaling_abba_G),
	from = 50, to = 800, ylim = c(0, 3.5), xlab = "Diameter (in mm)", ylab = "Growth (in mm/yr)", lwd = 3, col = "#1D4965")

curve(growth(x, temp = temperature_G_orford, precip = precipitations_G_orford, coeff_glm = acsa_G, scaling = scaling_acsa_G), add = TRUE,
	lwd = 3, col = "#2F90A8")
dev.off()

# ... dbh in New-Jersey
tikz("growth_vs_dbh_newJersey.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.8, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(growth(x, temp = temperature_G_newJersey, precip = precipitations_G_newJersey, coeff_glm = abba_G, scaling = scaling_abba_G),
	from = 50, to = 800, ylim = c(0, 4.3), xlab = "Diameter (in mm)", ylab = "Growth (in mm/yr)",
	lwd = 3, col = "#1D4965")

curve(growth(x, temp = temperature_G_newJersey, precip = precipitations_G_newJersey, coeff_glm = acsa_G, scaling = scaling_acsa_G),
	add = TRUE, lwd = 3, col = "#2F90A8")
dev.off()



## Mortality vs ...
# ... dbh in Orford
tikz("mortality_vs_dbh_orford.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.8, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(mortality(x, temp = temperature_M_orford, precip = precipitations_M_orford, coeff_glm = abba_M, scaling = scaling_abba_M),
	from = 50, to = 800, ylim = c(0, 0.05), xlab = "Diameter (in mm)", ylab = "Mortality (in $ yr^{-1} $)", lwd = 3, col = "#1D4965")

curve(mortality(x, temp = temperature_M_orford, precip = precipitations_M_orford, coeff_glm = acsa_M, scaling = scaling_acsa_M), add = TRUE,
	lwd = 3, col = "#2F90A8")
dev.off()

# ... dbh in Orford
tikz("mortality_vs_dbh_newJersey.tex", width = 2.6, height = 2.6)
par(mar = c(4, 4, 0.8, 0.8), mgp = c(2.8, 0.8, 0), tck = -0.02)
curve(mortality(x, temp = temperature_M_newJersey, precip = precipitations_M_newJersey, coeff_glm = abba_M, scaling = scaling_abba_M),
	from = 50, to = 800, ylim = c(0, 0.05), xlab = "Diameter (in mm)", ylab = "Mortality (in $ yr^{-1} $)",
	lwd = 3, col = "#1D4965")

curve(mortality(x, temp = temperature_M_newJersey, precip = precipitations_M_newJersey, coeff_glm = acsa_M, scaling = scaling_acsa_M),
	add = TRUE, lwd = 3, col = "#2F90A8")
dev.off()
