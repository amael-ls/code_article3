
#### Aim of prog: Create the growth and mortality parameters for C++ program
## Growth:
#		- 20 parameters (intercep, slopes of T, P, dbh, dbh², and cs + clim interactions with dbh, dbh², and cs)
#		- The parameters are stored in the article 1 folder on Beluga (Compute Canada)
#
## Mortality
#		- x parameters
#		- The parameters are stored in the article 1 folder on Beluga (Compute Canada)
#

#### Load packages
library(data.table)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

#### Tool functions
## Write parameters for C++ program, with cppNames for the C++ names
writeCppParams = function(params, cppNames, species, suffix, rm = FALSE)
{
	if (length(params) != length(cppNames))
		print("*** ERROR (from writeCppParams) ***: dimensions mismatch")

	outfileName = paste0("./", species, suffix, ".txt")
	if (rm & file.exists(outfileName))
		file.remove(outfileName)

	sink(file = outfileName, append = TRUE)
	for (i in 1:length(params))
		cat(paste0(cppNames[i], " = ", params[i]), sep = "\n")
	sink(file = NULL)
}

## Coerce C0_C1 to vector
makeC0_C1 = function(C0_C1, crownNames)
{
	C0_C1_names = paste0(rep(C0_C1[, parameter], each = 2), rep(c("_C0", "_C1"), 4))
	if (any(crownNames != C0_C1_names))
		print("*** ERROR (from makeC0_C1) ***: names mismatch")

	C0_C1_vec = numeric(2*C0_C1[, .N])
	for (i in 1:C0_C1[, .N])
		C0_C1_vec[(2*i-1):(2*i)] = unlist(C0_C1[i, .(C0, C1)])

	return (C0_C1_vec)
}

#### Common variables
## Folders
loadPath_G = "~/projects/def-dgravel/amael/article1/progToSendToReview/growth/"
loadPath_M = "~/projects/def-dgravel/amael/article1/progToSendToReview/mortality/"
ls_folders = list.files(loadPath_G, pattern = "array_[0-9]")

if (any(ls_folders != list.files(loadPath_M, pattern = "array_[0-9]")))
	print("*** ERROR (from parametersForCpp.R) ***: growth and mortality folders mismatch")

matlabPath = "~/projects/def-dgravel/amael/article1/progToSendToReview/createMatlabData/"

## Scalings
scaling_G = fread(paste0(matlabPath, "growthScaling.csv"))
scaling_dbh_G = fread(paste0(matlabPath, "growthDbhScaling.csv"))
scaling_temp_G = fread(paste0(matlabPath, "growthTempScaling.csv"))
scaling_precip_G = fread(paste0(matlabPath, "growthPrecipScaling.csv"))

scaling_dbh_M = fread(paste0(matlabPath, "mortalityDbhScaling.csv"))
scaling_temp_M = fread(paste0(matlabPath, "mortalityTempScaling.csv"))
scaling_precip_M = fread(paste0(matlabPath, "mortalityPrecipScaling.csv"))

## Allometries
abT_params = fread(paste0(matlabPath, "purves2007_allometries.csv"))
C0_C1 = fread(paste0(matlabPath, "C0_C1.csv"))
C0_C1 = C0_C1[!(parameter %in% c("Rus", "Vus"))]

## Fecundity (from Purves et al. 2008: http://www.pnas.org/content/105/44/17018.abstract)
fec = 0.0071

## Max diameter (it correspond to the infinite dbh, which is the dbh at 45 meters)
# cf article 1, to check the time it takes to grow up to 45m (folder article1/code/time)
dbh_bounds = fread(paste0(matlabPath, "dbh_params.csv"));

## C++ names
# Growth
cppNames_G = c("intercept", "cs", "dbh", "dbh_sq",
	"T", "T_sq", "P", "P_sq",
	"cs_T", "cs_T_sq", "cs_P", "cs_P_sq",
	"dbh_T", "dbh_T_sq", "dbh_P", "dbh_P_sq",
	"dbh_sq_T", "dbh_sq_T_sq", "dbh_sq_P", "dbh_sq_P_sq")

cppNames_scaling_growth_G = c("scaling_G_mu", "scaling_G_sd")
cppNames_scaling_dbh_G = c("scaling_dbh_mu_G", "scaling_dbh_sd_G")
cppNames_scaling_temp_G = c("scaling_temp_mu_G", "scaling_temp_sd_G")
cppNames_scaling_precip_G = c("scaling_precip_mu_G", "scaling_precip_sd_G")

cppNames_scaling_G = c(cppNames_scaling_growth_G, cppNames_scaling_dbh_G,
	cppNames_scaling_temp_G, cppNames_scaling_precip_G)

# Mortality
cppNames_M = c("intercept", "cs",
	"T", "T_sq", "P", "P_sq",
	 "dbh", "dbh_sq",
	"cs_T", "cs_T_sq", "cs_P", "cs_P_sq")

cppNames_scaling_dbh_M = c("scaling_dbh_mu_M", "scaling_dbh_sd_M")
cppNames_scaling_temp_M = c("scaling_temp_mu_M", "scaling_temp_sd_M")
cppNames_scaling_precip_M = c("scaling_precip_mu_M", "scaling_precip_sd_M")

cppNames_scaling_M = c(cppNames_scaling_dbh_M, cppNames_scaling_temp_M,
	cppNames_scaling_precip_M)

# Crown and height allometries
cppNames_height = c("a", "b", "T_param")
cppNames_crown = c("R0_C0", "R0_C1", "R40_C0", "R40_C1",
	"B_C0", "B_C1", "M_C0", "M_C1")

cppNames_allometries = c(cppNames_height, cppNames_crown)

## Coerce C0_C1 to vector
C0_C1 = makeC0_C1(C0_C1, cppNames_crown)

## Species
species_dt = readRDS("~/projects/def-dgravel/amael/article1/progToSendToReview/createData/speciesTable.rds")

#### Write files for Cpp
for (folder in ls_folders)
{
	species_txt = list.files(paste0(loadPath_G, folder), ".txt")
	species_txt = stri_sub(str = species_txt, to = stri_locate_last(str = species_txt, regex = ".txt")[1] - 1)
	species_scientific = species_dt[species_id == species_txt, latin]
	species_scientific = stri_replace_all(species_scientific, regex = " ", "_")

	# Growth parameters (growth function)
	params = readRDS(paste0(loadPath_G, folder, "/fixef_growth.rds"))
	writeCppParams(params, cppNames_G, species_scientific, "_G", rm = TRUE)

	# Mortality parameters
	params = readRDS(paste0(loadPath_M, folder, "/fixef.rds"))
	writeCppParams(params, cppNames_M, species_scientific, "_M", rm = TRUE)

	# Allometries
	params = c(unlist(abT_params[species == species_txt, .(a, b, T)]), C0_C1)
	writeCppParams(params, cppNames_allometries, species_scientific, "_allometries", rm = TRUE)

	# Scaling (growth)
	params = unlist(scaling_G[species_id == species_txt, .(mu, sd)])
	params = c(params, unlist(scaling_dbh_G[species_id == species_txt, .(mu, sd)]))
	params = c(params, unlist(scaling_temp_G[species_id == species_txt, .(mu, sd)]))
	params = c(params, unlist(scaling_precip_G[species_id == species_txt, .(mu, sd)]))
	writeCppParams(params, cppNames_scaling_G, species_scientific, "_scaling_G", rm = TRUE)

	# Scaling (mortality)
	params = unlist(scaling_dbh_M[species_id == species_txt, .(mu, sd)])
	params = c(params, unlist(scaling_temp_M[species_id == species_txt, .(mu, sd)]))
	params = c(params, unlist(scaling_precip_M[species_id == species_txt, .(mu, sd)]))
	writeCppParams(params, cppNames_scaling_M, species_scientific, "_scaling_M", rm = TRUE)

	# Species name and fecundity parameter
	outfileName = paste0("./", species_scientific, ".txt")
	if (file.exists(outfileName))
		file.remove(outfileName)

	sink(file = outfileName, append = TRUE)
	cat(paste0("species = ", species_scientific), sep = "\n")
	cat(paste0("fecundity = ", fec), sep = "\n")
	cat(paste0("maxDiameter = ", dbh_bounds[species_id == species_txt, dbh_infinity]), sep = "\n")
	sink(file = NULL)

	# Species dispersal parameters
	outfileName = paste0("./", species_scientific, "_dispersal.txt")
	if (file.exists(outfileName))
		file.remove(outfileName)

	keysToRead = c("dispersalThreshold", "refKernel_doi", "propLDD", "relLDDtoSDD")
	sink(file = outfileName, append = TRUE)
	cat(paste0("keysToRead = ", paste0(keysToRead, collapse = ", ")), sep = "\n")
	cat(paste0("dispersalProbaThresold = ", 0.01), sep = "\n")
	cat(paste0("refKernel_doi = ", "10.1016/j.jtbi.2005.12.019"), sep = "\n")
	cat(paste0("propLDD = ", 0.1), sep = "\n")
	cat(paste0("relLDDtoSDD = ", 0.1), sep = "\n")
	cat(paste0("dispersalDistThresold = ", 5), sep = "\n") # Maximal distance dispersal in km
	sink(file = NULL)

	print(paste0("Species: ", species_scientific, " done"))
}
