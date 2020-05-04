
#### Aim of prog: Create the simulation parameters for the C++ prog

#### Load packages
library(data.table)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

#### Parameters
## File's name
outfileName = "./simulationParameters.txt"

## List of the simulation params
simulation_ls = list(
	maxCohorts = 3000,
	n_t = 5000,
	t0 = 0,
	t_max = 500,
	climate_file = "../createParams/clim_1.txt",
	species_file = "../createParams/Populus_tremuloides.txt",
	init_file = "../createIC/ic_1.txt")

varNames = names(simulation_ls)

#### Write the simultation parameters
## Remove file if exists
if (file.exists(outfileName))
	file.remove(outfileName)

## Change the std::out
sink(file = outfileName, append = TRUE)

## Fill
for (i in 1:length(simulation_ls))
{
	cat(paste0(varNames[i], ": ", simulation_ls[[varNames[i]]]), eol = "\n")
}

## Come back to classic std::out
sink(file = NULL)
