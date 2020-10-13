
#ifndef SPECIES_C
#define SPECIES_C

// My headers
#include "Species.h++"
#include "Params.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Species::Species(std::string const& species_filename, std::string const& species_path, const std::string& delim)
{
	// Read species' main file from the parameters got at the execution
	std::string file_sp = species_path + species_filename;
	par::Params speciesParams(file_sp.c_str(), delim);

	// Species name
	m_speciesName = speciesParams.get_val<std::string>("species");

	// Create the files' name for a given species
	std::string file_G = species_path + m_speciesName;
	file_G.append("_G.txt");

	std::string file_M = species_path + m_speciesName;
	file_M.append("_M.txt");

	std::string file_allometries = species_path + m_speciesName;
	file_allometries.append("_allometries.txt");

	std::string file_scaling_G = species_path + m_speciesName;
	file_scaling_G.append("_scaling_G.txt");

	std::string file_scaling_M = species_path + m_speciesName;
	file_scaling_M.append("_scaling_M.txt");

	std::string file_scaling_dispersal = species_path + m_speciesName;
	file_scaling_dispersal.append("_dispersal.txt");

	// Load parameters from files
	par::Params speciesParams_G(file_G.c_str(), delim);
	par::Params speciesParams_M(file_M.c_str(), delim);
	par::Params speciesParams_allometries(file_allometries.c_str(), delim);
	par::Params speciesParams_scaling_G(file_scaling_G.c_str(), delim);
	par::Params speciesParams_scaling_M(file_scaling_M.c_str(), delim);
	par::Params speciesParams_dispersal(file_scaling_dispersal.c_str(), delim);

	// Allometry parameters (height from dbh)
	a = speciesParams_allometries.get_val<double>("a");
	b = speciesParams_allometries.get_val<double>("b");

	// Allometry parameters (crown area from dbh)
	T_param = speciesParams_allometries.get_val<double>("T_param");
	R0_C0 = speciesParams_allometries.get_val<double>("R0_C0");
	R0_C1 = speciesParams_allometries.get_val<double>("R0_C1");
	R40_C0 = speciesParams_allometries.get_val<double>("R40_C0");
	R40_C1 = speciesParams_allometries.get_val<double>("R40_C1");
	M_C0 = speciesParams_allometries.get_val<double>("M_C0");
	M_C1 = speciesParams_allometries.get_val<double>("M_C1");
	B_C0 = speciesParams_allometries.get_val<double>("B_C0");
	B_C1 = speciesParams_allometries.get_val<double>("B_C1");

	// Growth parameters (20)
	intercept_G = speciesParams_G.get_val<double>("intercept");
	beta_dbh = speciesParams_G.get_val<double>("dbh");
	beta_dbh_T = speciesParams_G.get_val<double>("dbh_T");
	beta_dbh_T_sq = speciesParams_G.get_val<double>("dbh_T_sq");
	beta_dbh_P = speciesParams_G.get_val<double>("dbh_P");
	beta_dbh_P_sq = speciesParams_G.get_val<double>("dbh_P_sq");
	beta_dbh_sq = speciesParams_G.get_val<double>("dbh_sq");
	beta_dbh_sq_T = speciesParams_G.get_val<double>("dbh_sq_T");
	beta_dbh_sq_T_sq = speciesParams_G.get_val<double>("dbh_sq_T_sq");
	beta_dbh_sq_P = speciesParams_G.get_val<double>("dbh_sq_P");
	beta_dbh_sq_P_sq = speciesParams_G.get_val<double>("dbh_sq_P_sq");
	beta_cs = speciesParams_G.get_val<double>("cs");
	beta_cs_T = speciesParams_G.get_val<double>("cs_T");
	beta_cs_T_sq = speciesParams_G.get_val<double>("cs_T_sq");
	beta_cs_P = speciesParams_G.get_val<double>("cs_P");
	beta_cs_P_sq = speciesParams_G.get_val<double>("cs_P_sq");
	beta_T = speciesParams_G.get_val<double>("T");
	beta_T_sq = speciesParams_G.get_val<double>("T_sq");
	beta_P = speciesParams_G.get_val<double>("P");
	beta_P_sq = speciesParams_G.get_val<double>("P_sq");

	// Mortality parameters (12 + 1)
	intercept_M = speciesParams_M.get_val<double>("intercept");
	beta_cs_M = speciesParams_M.get_val<double>("cs");
	beta_T_M = speciesParams_M.get_val<double>("T");
	beta_T_sq_M = speciesParams_M.get_val<double>("T_sq");
	beta_P_M = speciesParams_M.get_val<double>("P");
	beta_P_sq_M = speciesParams_M.get_val<double>("P_sq");
	beta_dbh_M = speciesParams_M.get_val<double>("dbh");
	beta_dbh_sq_M = speciesParams_M.get_val<double>("dbh_sq");
	beta_cs_T_M = speciesParams_M.get_val<double>("cs_T");
	beta_cs_T_sq_M = speciesParams_M.get_val<double>("cs_T_sq");
	beta_cs_P_M = speciesParams_M.get_val<double>("cs_P");
	beta_cs_P_sq_M = speciesParams_M.get_val<double>("cs_P_sq");

	// Scaling (growth)
	scaling_G_mu = speciesParams_scaling_G.get_val<double>("scaling_G_mu");
	scaling_G_sd = speciesParams_scaling_G.get_val<double>("scaling_G_sd");
	scaling_dbh_mu_G = speciesParams_scaling_G.get_val<double>("scaling_dbh_mu_G");
	scaling_dbh_sd_G = speciesParams_scaling_G.get_val<double>("scaling_dbh_sd_G");
	scaling_temp_mu_G = speciesParams_scaling_G.get_val<double>("scaling_temp_mu_G");
	scaling_temp_sd_G = speciesParams_scaling_G.get_val<double>("scaling_temp_sd_G");
	scaling_precip_mu_G = speciesParams_scaling_G.get_val<double>("scaling_precip_mu_G");
	scaling_precip_sd_G = speciesParams_scaling_G.get_val<double>("scaling_precip_sd_G");

	// Scaling (mortality)
	scaling_dbh_mu_M = speciesParams_scaling_M.get_val<double>("scaling_dbh_mu_M");
	scaling_dbh_sd_M = speciesParams_scaling_M.get_val<double>("scaling_dbh_sd_M");
	scaling_temp_mu_M = speciesParams_scaling_M.get_val<double>("scaling_temp_mu_M");
	scaling_temp_sd_M = speciesParams_scaling_M.get_val<double>("scaling_temp_sd_M");
	scaling_precip_mu_M = speciesParams_scaling_M.get_val<double>("scaling_precip_mu_M");
	scaling_precip_sd_M = speciesParams_scaling_M.get_val<double>("scaling_precip_sd_M");

	// Fecundity
	// --- Get provided parameters
	fecundity = speciesParams.get_val<double>("fecundity");
	minAgeReproduction = speciesParams.get_val<double>("minAgeReproduction");

	// --- Compute minimal size to reproduce
	// ...... Define void *ptr for parameters
	double allParameters[24] = {scaling_G_mu, scaling_G_sd,
		scaling_dbh_mu_G, scaling_dbh_sd_G,
		scaling_temp_mu_G, scaling_temp_sd_G,
		scaling_precip_mu_G, scaling_precip_sd_G,
		intercept_G, beta_dbh, beta_dbh_sq,
		beta_T, beta_T_sq, beta_dbh_T, beta_dbh_T_sq, beta_dbh_sq_T, beta_dbh_sq_T_sq,
		beta_P, beta_P_sq, beta_dbh_P, beta_dbh_P_sq, beta_dbh_sq_P, beta_dbh_sq_P_sq};
	void *ptr = &allParameters;

	// ...... Define solver parameters
	std::string timeSpan = "[0, " + std::to_string(minAgeReproduction) + "]";
	alglib::real_1d_array timeSpanArray(timeSpan.c_str());
	alglib::real_1d_array y = "[0]"; // Initial condition
	double eps = 0.00001;
	double initialTimeStep = 0;
	alglib::odesolverstate s;
	alglib::ae_int_t m;
	alglib::real_1d_array xtbl;
	alglib::real_2d_array ytbl;
	alglib::odesolverreport rep;
	alglib::odesolverrkck(y, timeSpanArray, eps, initialTimeStep, s);

	// ...... Solve ds/dt = G(s, t)
	alglib::odesolversolve(s, Species::growth_callback, ptr);
	odesolverresults(s, m, xtbl, ytbl, rep);

	// ...... Set minimal size to reproduce
	if (*(ytbl[1]) < 0 || m != 2)
		throw Except_Species(m, *(ytbl[1]));
	
	minDbhReproduction = *(ytbl[1]);
	minHeightReproduction = std::exp((a - b + b*std::log10(minDbhReproduction))*std::log(10));

	// Dispersal parameters, not necessarily all provided by the user. Note: I could have used a map to optimise the coding...
	keysToRead = speciesParams_dispersal.get_val<std::string>("keysToRead");
	min_dispersalProba = false;
	max_dispersalDist = false;
	if (keysToRead.find("dispersalProbaThreshold") != std::string::npos)
	{
		dispersalProbaThreshold = speciesParams_dispersal.get_val<double>("dispersalProbaThreshold");
		min_dispersalProba = true;
	}

	if (keysToRead.find("refKernel_doi") != std::string::npos)
		refKernel_doi = speciesParams_dispersal.get_val<std::string>("refKernel_doi");

	if (keysToRead.find("propLDD") != std::string::npos)
		propLDD = speciesParams_dispersal.get_val<double>("propLDD");

	if (keysToRead.find("relLDDtoSDD") != std::string::npos)
		relLDDtoSDD = speciesParams_dispersal.get_val<double>("relLDDtoSDD");

	if (keysToRead.find("dispersalDistThreshold") != std::string::npos)
	{
		dispersalDistThreshold = speciesParams_dispersal.get_val<double>("dispersalDistThreshold");
		max_dispersalDist = true;
	}

	if (keysToRead.find("twoDt_a") != std::string::npos)
		twoDt_a = speciesParams_dispersal.get_val<double>("twoDt_a");

	if (keysToRead.find("twoDt_b") != std::string::npos)
		twoDt_b = speciesParams_dispersal.get_val<double>("twoDt_b");

	// Others
	maxDiameter = speciesParams.get_val<double>("maxDiameter");
}

/**************************************/
/******        Demography        ******/
/**************************************/
// Individual growth rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::v(double s, double const s_star, double temp, double precip) const
{
	/*
		It is assumed that none of the variables is scaled.
	*/
	bool cs = s_star < s; // ? false : true;
	double beta_0, beta_1, beta_2;

	// Scaling all the variables (size, temperature and precipitation)
	s = (s - scaling_dbh_mu_G)/scaling_dbh_sd_G;
	temp = (temp - scaling_temp_mu_G)/scaling_temp_sd_G;
	precip = (precip - scaling_precip_mu_G)/scaling_precip_sd_G;

	// Define the coefficient of the polynomial of s
	beta_0 = intercept_G +
	(beta_cs + beta_cs_T*temp + beta_cs_T_sq*temp*temp +
		beta_cs_P*precip + beta_cs_P_sq*precip*precip)*cs +
	beta_T*temp + beta_T_sq*temp*temp + beta_P*precip + beta_P_sq*precip*precip;

	beta_1 = beta_dbh + beta_dbh_T*temp + beta_dbh_T_sq*temp*temp +
		beta_dbh_P*precip + beta_dbh_P_sq*precip*precip;

	beta_2 = beta_dbh_sq + beta_dbh_sq_T*temp + beta_dbh_sq_T_sq*temp*temp +
		beta_dbh_sq_P*precip + beta_dbh_sq_P_sq*precip*precip;

	// Growth function
	return std::exp(scaling_G_mu + scaling_G_sd * (beta_0 + beta_1*s + beta_2*s*s));
}

// Individual death rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::d(double s, double const s_star, double temp, double precip) const
{
	/*
		It is assumed that none of the variables is scaled.
	*/

	bool cs = s_star < s; // ? false : true;
	double beta_0;

	// Scaling all the variables (size, temperature and precipitation)
	s = (s - scaling_dbh_mu_M)/scaling_dbh_sd_M;
	temp = (temp - scaling_temp_mu_M)/scaling_temp_sd_M;
	precip = (precip - scaling_precip_mu_M)/scaling_precip_sd_M;

	// Define the coefficient of the polynomial of s
	beta_0 = intercept_M +
	(beta_cs_M + beta_cs_T_M*temp + beta_cs_T_sq_M*temp*temp
		+ beta_cs_P_M*precip + beta_cs_P_sq_M*precip*precip)*cs +
	beta_T_M*temp + beta_T_sq_M*temp*temp +beta_P_M*precip + beta_P_sq_M*precip*precip;

	return std::exp(beta_0 + beta_dbh_M*s + beta_dbh_sq_M*s*s);
	// return 1 - std::exp(-std::exp(beta_0 + beta_dbh_M*s + beta_dbh_sq_M*s*s));
}

// Differentiate individual growth rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::dv_ds(double s, double const s_star, double temp, double precip) const
{
	/*
		It is assumed that none of the variables is scaled.
		This function is only weakly differentiable, due to the jump at s_star (hereafter s*)
		Its weak derivative is (δ_0(s - s*) + β_1 + 2 β_2 s) Exp[β_0 + β_1 s + β_2 s^2]
		where the βs are defined in v.

		However, I guess I can also work in [0, s*[ and ]s*, +Inf[, and set a continuity condition at s*.
		This was done in the paper of Adams et al. 2007:
		Understanding height-structured competition in forests: is there an R* for light?
		the total rate at which cohorts 'arrive' at size s* as canopy trees must equal the rate at which they
		'leave' s* as understory trees.
	*/

	bool cs = s_star < s; // ? false : true;
	double beta_0, beta_1, beta_2;

	// Scaling all the variables (size, temperature and precipitation)
	s = (s - scaling_dbh_mu_G)/scaling_dbh_sd_G;
	temp = (temp - scaling_temp_mu_G)/scaling_temp_sd_G;
	precip = (precip - scaling_precip_mu_G)/scaling_precip_sd_G;

	// Define the coefficient of the polynomial of s
	beta_0 = intercept_G +
	(beta_cs + beta_cs_T*temp + beta_cs_T_sq*temp*temp +
		beta_cs_P*precip + beta_cs_P_sq*precip*precip)*cs +
	beta_T*temp + beta_T_sq*temp*temp + beta_P*precip + beta_P_sq*precip*precip;

	beta_1 = beta_dbh + beta_dbh_T*temp + beta_dbh_T_sq*temp*temp +
		beta_dbh_P*precip + beta_dbh_P_sq*precip*precip;

	beta_2 = beta_dbh_sq + beta_dbh_sq_T*temp + beta_dbh_sq_T_sq*temp*temp +
		beta_dbh_sq_P*precip + beta_dbh_sq_P_sq*precip*precip;

	// Polynomial of s (order 2)
	double dbh_polynomial = beta_0 + beta_1*s + beta_2*s*s;

	return scaling_G_sd*(beta_1 + 2*beta_2*s) * std::exp(scaling_G_mu + scaling_G_sd * dbh_polynomial);
}

// Differentiate individual mortality rate, s is the diameter, s_star is dbh_star (obtained from height_star)
double Species::dd_ds(double s, double const s_star, double temp, double precip) const
{
	/*
		It is assumed that none of the variables is scaled.
		This function is only weakly differentiable, due to the jump at s_star (hereafter s*)
		Its weak derivative is + (δ_0(s - s*) + β_1 + 2 β_2 s) Exp[-u]/(1 + Exp[-u])^2
		where the βs are defined in d, and u = β_0 + β_1 s + β_2 s^2.
		an equivalent expression is: u'/(2 + Exp[-u] + Exp[u]), which tends to 0 when u tends to +Inf

		However, I guess I can also work in [0, s*[ and ]s*, +Inf[, and set a continuity condition at s*.
		This was done in the paper of Adams et al. 2007:
		Understanding height-structured competition in forests: is there an R* for light?
		the total rate at which cohorts 'arrive' at size s* as canopy trees must equal the rate at which they
		'leave' s* as understory trees.
	*/

	bool cs = s_star < s; // ? false : true;
	double beta_0;

	// Scaling all the variables (size, temperature and precipitation)
	s = (s - scaling_dbh_mu_M)/scaling_dbh_sd_M;
	temp = (temp - scaling_temp_mu_M)/scaling_temp_sd_M;
	precip = (precip - scaling_precip_mu_M)/scaling_precip_sd_M;

	// Define the coefficient of the polynomial of s
	beta_0 = intercept_M +
	(beta_cs_M + beta_cs_T_M*temp + beta_cs_T_sq_M*temp*temp
		+ beta_cs_P_M*precip + beta_cs_P_sq_M*precip*precip)*cs +
	beta_T_M*temp + beta_T_sq_M*temp*temp +beta_P_M*precip + beta_P_sq_M*precip*precip;

	// Polynomial of s (order 2)
	double dbh_polynomial = beta_0 + beta_dbh_M*s + beta_dbh_sq_M*s*s;

	return (beta_dbh_M + 2*beta_dbh_sq_M*s)*std::exp(dbh_polynomial);
	// return -(beta_dbh_M + 2*beta_dbh_sq_M*s)*std::exp(dbh_polynomial - std::exp(dbh_polynomial));
}

/*************************************/
/******        Dispersal        ******/
/*************************************/
// From Moorcroft & Lewis 2006: Potential role of natural enemies during tree range expansions following climate change.
// Rk: There kernel has small mistakes that have been corrected here
double Species::K(double const distance) const
{
	double proba = 0;
	if (refKernel_doi == "10.1016/j.jtbi.2005.12.019") // Moorcroft2006
		proba = (1.0 - propLDD)/2.0*exp(-std::abs(distance)) + propLDD*relLDDtoSDD/2.0*exp(-relLDDtoSDD*distance);

	if (refKernel_doi == "laplacian") // Laplacian with parameter = 100 (for testing)
		proba = 1.0/(2*M_PI*100*100) * std::exp(- distance/100);

	if (refKernel_doi == "10.2307/176541") // Clark1999, Boisvert-Marsh2020 (might be 2021)
		proba = twoDt_a/(M_PI*twoDt_b) * std::exp((-twoDt_a - 1)*std::log(1 + distance*distance/twoDt_b));
	
	return proba;
}

double Species::K(double delta_lon, double delta_lat) const
{
	double proba = 0;
	double distance = sqrt(delta_lon*delta_lon + delta_lat*delta_lat);
	if (refKernel_doi == "10.1016/j.jtbi.2005.12.019") // Moorcroft2006
		proba = (1.0 - propLDD)/2.0*exp(-std::abs(distance)) + propLDD*relLDDtoSDD/2.0*exp(-relLDDtoSDD*distance);
	
	if (refKernel_doi == "laplacian") // Laplacian with parameter = 100 (for testing)
		proba = 1.0/(2*M_PI*100*100) * std::exp(- distance/100);

	if (refKernel_doi == "10.2307/176541") // Clark1999, Boisvert-Marsh2020 (might be 2021)
		proba = twoDt_a/(M_PI*twoDt_b) * std::exp((-twoDt_a - 1)*std::log(1 + distance*distance/twoDt_b));

	return proba;
}

double Species::K(double const longitude1, double const latitude1, double const longitude2, double const latitude2) const
{
	double proba = 0;
	double distance = sqrt((longitude1 - longitude2)*(longitude1 - longitude2) + (latitude1 - latitude2)*(latitude1 - latitude2));
	if (refKernel_doi == "10.1016/j.jtbi.2005.12.019") // Moorcroft2006
		proba = (1.0 - propLDD)/2.0*exp(-std::abs(distance)) + propLDD*relLDDtoSDD/2.0*exp(-relLDDtoSDD*distance);
	
	return proba;
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream &operator<<(std::ostream &os, Species const& species)
{
	os << std::setprecision(5);
	os << species.m_speciesName << std::endl;
	os << std::endl;

	os << "Allometry parameters height = f(dbh):" << std::endl;
	os << species.a << "\t" << species.b << std::endl;
	os << std::endl;

	os << "Allometry parameters crownArea = f(dbh):" << std::endl;
	os << species.T_param << "\t" << species.R0_C0 << "\t" << species.R0_C1 << "\t" << species.R40_C0
	<< "\t" << species.R40_C1 << "\t" << species.M_C0 << "\t" << species.M_C1
	<< "\t" << species.B_C0 << "\t" << species.B_C1 << std::endl;
	os << std::endl;

	os << "Growth parameters (intercept):" << std::endl << species.intercept_G << std::endl;
	os << std::endl;

	os << "Growth parameters (dbh and climate interactions):" << std::endl;
	os << species.beta_dbh << "\t" << species.beta_dbh_T << "\t" << species.beta_dbh_T_sq
		<< "\t" << species.beta_dbh_P << "\t" << species.beta_dbh_P_sq << std::endl;
	os << std::endl;

	os << "Growth parameters (dbh² and climate interactions):" << std::endl;
	os << species.beta_dbh_sq << "\t" << species.beta_dbh_sq_T << "\t" << species.beta_dbh_sq_T_sq
		<< "\t" << species.beta_dbh_sq_P << "\t" << species.beta_dbh_sq_P_sq << std::endl;
	os << std::endl;

	os << "Growth parameters (canopy status and climate interactions):" << std::endl;
	os << species.beta_cs << "\t" << species.beta_cs_T << "\t" << species.beta_cs_T_sq
		<< "\t" << species.beta_cs_P << "\t" << species.beta_cs_P_sq << std::endl;
	os << std::endl;

	os << "Growth parameters (climate):" << std::endl;
	os << species.beta_T << "\t" << species.beta_T_sq
		<< "\t" << species.beta_P << "\t" << species.beta_P_sq << std::endl;
	os << std::endl;

	os << "Mortality parameters (intercept): " << std::endl << species.intercept_M << std::endl;
	os << std::endl;

	os << "Mortality parameters (dbh and dbh²):" << std::endl;
	os << species.beta_dbh_M << "\t" << species.beta_dbh_sq_M << std::endl;
	os << std::endl;

	os << "Mortality parameters (canopy status and climate interactions):" << std::endl;
	os << species.beta_cs_M << "\t" << species.beta_cs_T_M << "\t" << species.beta_cs_T_sq_M << "\t"
		<< species.beta_cs_P_M << "\t" << species.beta_cs_P_sq_M << std::endl;
	os << std::endl;

	os << "Mortality parameters (climate):" << std::endl;
	os << species.beta_T_M << "\t" << species.beta_T_sq_M << "\t" << species.beta_P_M << "\t"
		<< species.beta_P_sq_M << "\t" << std::endl;
	os << std::endl;

	os << "Fecundity parameters:" << std::endl;
	os << species.fecundity << "\t" << species.minHeightReproduction << std::endl;
	os << std::endl;

	os << "Scaling (growth):" << std::endl;
	os << species.scaling_G_mu << "\t" << species.scaling_G_sd << std::endl;
	os << std::endl;

	os << "Scaling (growth dbh):" << std::endl;
	os << species.scaling_dbh_mu_G << "\t" << species.scaling_dbh_sd_G << std::endl;
	os << std::endl;

	os << "Scaling (growth temp):" << std::endl;
	os << species.scaling_temp_mu_G << "\t" << species.scaling_temp_sd_G << std::endl;
	os << std::endl;

	os << "Scaling (growth precip):" << std::endl;
	os << species.scaling_precip_mu_G << "\t" << species.scaling_precip_sd_G << std::endl;
	os << std::endl;

	os << "Scaling (mortality dbh):" << std::endl;
	os << species.scaling_dbh_mu_M << "\t" << species.scaling_dbh_sd_M << std::endl;
	os << std::endl;

	os << "Scaling (mortality temp):" << std::endl;
	os << species.scaling_temp_mu_M << "\t" << species.scaling_temp_sd_M << std::endl;
	os << std::endl;

	os << "Scaling (mortality precip):" << std::endl;
	os << species.scaling_precip_mu_M << "\t" << species.scaling_precip_sd_M << std::endl;
	os << std::endl;

	os << "Dispersal parameters:" << std::endl;
	if (species.keysToRead.find("dispersalProbaThreshold") != std::string::npos)
		os << "dispersalProbaThreshold: " << species.dispersalProbaThreshold << std::endl;
	if (species.keysToRead.find("refKernel_doi") != std::string::npos)
		os << "refKernel_doi: " << species.refKernel_doi << std::endl;
	if (species.keysToRead.find("propLDD") != std::string::npos)
		os << "propLDD: " << species.propLDD << std::endl;
	if (species.keysToRead.find("relLDDtoSDD") != std::string::npos)
		os << "relLDDtoSDD: " << species.relLDDtoSDD << std::endl;
	if (species.keysToRead.find("dispersalDistThreshold") != std::string::npos)
		os << "dispersalDistThreshold: " << species.dispersalDistThreshold << std::endl;

	os << "Maximal diameter (correspond to a 45m height):" << std::endl;
	os << species.maxDiameter << std::endl;
	os << std::endl;

	return os;
}

// Sort by alphabetical order
bool operator<(Species const& species1, Species const& species2)
{
	return (species1.m_speciesName < species2.m_speciesName);
}

void Species::printName(std::ostream& os) const
{
	os << m_speciesName;
}

// Getter
std::string Species::getName() const
{
	return m_speciesName;
}

void Species::growth_callback(const alglib::real_1d_array &y, double x, alglib::real_1d_array &dy, void *ptr) 
{
	/*
		x = time
		y = function to find dy/dx = growth_callback(x, y(x))

		It is assumed that size is unscaled.
		It is assumed that the tree is in the understorey, i.e., canopy status = false
		It is assumed that temperature and precipitation are at average, i.e. = 0 when scaled
	*/
	// Access parameters from void *ptr
	// --- Casting
	double *parameters = (double*) ptr;

	// --- Growth scaling
	double scaling_G_mu = parameters[0];
	double scaling_G_sd = parameters[1];

	// --- Size scaling
	double scaling_dbh_mu_G = parameters[2];
	double scaling_dbh_sd_G = parameters[3];

	// --- Temperature scaling
	double scaling_temp_mu_G = parameters[5];
	double scaling_temp_sd_G = parameters[6];

	// --- Precipitation scaling
	double scaling_precip_mu_G = parameters[7];
	double scaling_precip_sd_G = parameters[8];

	// --- Coefficients (size y)
	double intercept_G = parameters[9];
	double beta_dbh = parameters[10];
	double beta_dbh_sq = parameters[11];

	// --- Coefficients (temperature)
	double beta_T = parameters[12];
	double beta_T_sq = parameters[13];
	double beta_dbh_T = parameters[14];
	double beta_dbh_T_sq = parameters[15];
	double beta_dbh_sq_T = parameters[16];
	double beta_dbh_sq_T_sq = parameters[17];

	// --- Coefficients (precipitation)
	double beta_P = parameters[18];
	double beta_P_sq = parameters[19];
	double beta_dbh_P = parameters[20];
	double beta_dbh_P_sq = parameters[21];
	double beta_dbh_sq_P = parameters[22];
	double beta_dbh_sq_P_sq = parameters[23];

	// Scaling all the variables (size y, temperature temp, and precipitation precip)
	double temp = 0; // (temp - scaling_temp_mu_G)/scaling_temp_sd_G;
	double precip = 0; // (precip - scaling_precip_mu_G)/scaling_precip_sd_G;

	// Define the coefficient of the polynomial of s
	double beta_0 = intercept_G +
	beta_T*temp + beta_T_sq*temp*temp + beta_P*precip + beta_P_sq*precip*precip;

	double beta_1 = beta_dbh + beta_dbh_T*temp + beta_dbh_T_sq*temp*temp +
		beta_dbh_P*precip + beta_dbh_P_sq*precip*precip;

	double beta_2 = beta_dbh_sq + beta_dbh_sq_T*temp + beta_dbh_sq_T_sq*temp*temp +
		beta_dbh_sq_P*precip + beta_dbh_sq_P_sq*precip*precip;

	// Growth function
	dy[0] = std::exp(scaling_G_mu + scaling_G_sd * (beta_0 + beta_1 * (y[0] - scaling_dbh_mu_G)/scaling_dbh_sd_G
		+ beta_2 * (y[0] - scaling_dbh_mu_G)/scaling_dbh_sd_G * (y[0] - scaling_dbh_mu_G)/scaling_dbh_sd_G));
}

#endif
