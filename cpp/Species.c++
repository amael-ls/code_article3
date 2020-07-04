
#ifndef SPECIES_C
#define SPECIES_C

// Official headers
#include <iomanip> // std::setw, std::left, std::setprecision
#include <cmath> // for log

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

	// Load parameters from files
	par::Params speciesParams_G(file_G.c_str(), delim);
	par::Params speciesParams_M(file_M.c_str(), delim);
	par::Params speciesParams_allometries(file_allometries.c_str(), delim);
	par::Params speciesParams_scaling_G(file_scaling_G.c_str(), delim);
	par::Params speciesParams_scaling_M(file_scaling_M.c_str(), delim);

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

	// Mortality parameters
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

	// Fecundity
	fecundity = speciesParams.get_val<double>("fecundity");

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

	// Others
	maxDiameter = speciesParams.get_val<double>("maxDiameter");
}

/**************************************/
/******        Demography        ******/
/**************************************/
// Individual growth rate (Checked, it equals the lme4 prediction)
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

	// Define the coefficient of the polynom of s
	beta_0 = intercept_G +
	(beta_cs + beta_cs_T*temp + beta_cs_T_sq*temp*temp +
		beta_cs_P*precip + beta_cs_P_sq*precip*precip)*cs +
	beta_T*temp + beta_T_sq*temp*temp + beta_P*precip + beta_P_sq*precip*precip;

	beta_1 = beta_dbh + beta_dbh_T*temp + beta_dbh_T_sq*temp*temp +
		beta_dbh_P*precip + beta_dbh_P_sq*precip*precip;

	beta_2 = beta_dbh_sq + beta_dbh_sq_T*temp + beta_dbh_sq_T_sq*temp*temp +
		beta_dbh_sq_P*precip + beta_dbh_sq_P_sq*precip*precip;

	// Polynom of s (order 2)
	double dbh_polynom = beta_0 + beta_1*s + beta_2*s*s;

	// Growth function
	return 1; // std::exp(scaling_G_mu + scaling_G_sd * dbh_polynom);

	// double results = ((2 + s_star)*std::exp(-s_star))/(1 + s); // with feedback loop
	// double results = 2/(1 + s); // without feedback loop
	// return (results);
}

// Individual death rate // mean(aa) = 0.7318187
double Species::d(double s, double const s_star, double temp, double precip) const
{
	/*
		It is assumed that none of the variables is scaled.
	*/

	bool cs = s_star < s; // ? false : true;
	double beta_0, beta_1, beta_2;

	// Scaling all the variables (size, temperature and precipitation)
	s = (s - scaling_dbh_mu_M)/scaling_dbh_sd_M;
	temp = (temp - scaling_temp_mu_M)/scaling_temp_sd_M;
	precip = (precip - scaling_precip_mu_M)/scaling_precip_sd_M;

	// Define the coefficient of the polynom of s
	beta_0 = intercept_M +
	(beta_cs_M + beta_cs_T_M*temp + beta_cs_T_sq_M*temp*temp
		+ beta_cs_P_M*precip + beta_cs_P_sq_M*precip*precip)*cs +
	beta_T_M*temp + beta_T_sq_M*temp*temp +beta_P_M*precip + beta_P_sq_M*precip*precip;

	beta_1 = beta_dbh_M;

	beta_2 = beta_dbh_sq_M;

	// Polynom of s (order 2)
	double dbh_polynom = beta_0 + beta_1*s + beta_2*s*s;
	return 0; // 1 / (1 + std::exp(-dbh_polynom));

	// return (-0.3);
}

// Differentiate individual growth rate
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
		the total rate at which cohorts ‘arrive’ at size s* as canopy trees must equal the rate at which they
		‘leave’ s* as understory trees.
	*/

	bool cs = s_star < s; // ? false : true;
	double beta_0, beta_1, beta_2;

	// Scaling all the variables (size, temperature and precipitation)
	s = (s - scaling_dbh_mu_G)/scaling_dbh_sd_G;
	temp = (temp - scaling_temp_mu_G)/scaling_temp_sd_G;
	precip = (precip - scaling_precip_mu_G)/scaling_precip_sd_G;

	// Define the coefficient of the polynom of s
	beta_0 = intercept_G +
	(beta_cs + beta_cs_T*temp + beta_cs_T_sq*temp*temp +
		beta_cs_P*precip + beta_cs_P_sq*precip*precip)*cs +
	beta_T*temp + beta_T_sq*temp*temp + beta_P*precip + beta_P_sq*precip*precip;

	beta_1 = beta_dbh + beta_dbh_T*temp + beta_dbh_T_sq*temp*temp +
		beta_dbh_P*precip + beta_dbh_P_sq*precip*precip;

	beta_2 = beta_dbh_sq + beta_dbh_sq_T*temp + beta_dbh_sq_T_sq*temp*temp +
		beta_dbh_sq_P*precip + beta_dbh_sq_P_sq*precip*precip;

	// Polynom of s (order 2)
	double dbh_polynom = beta_0 + beta_1*s + beta_2*s*s;

	return 0; // scaling_G_sd*(beta_1 + 2*beta_2*s) * std::exp(scaling_G_mu + scaling_G_sd * dbh_polynom);

	// double results = -((2 + s_star)*std::exp(-s_star))/((1 + s)*(1 + s)); // with feedback loop
	// double results = -2/((1 + s)*(1 + s)); // without feedback loop
	// return (results);
}

// Differentiate individual mortality rate
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
		the total rate at which cohorts ‘arrive’ at size s* as canopy trees must equal the rate at which they
		‘leave’ s* as understory trees.
	*/

	bool cs = s_star < s; // ? false : true;
	double beta_0, beta_1, beta_2;

	// Scaling all the variables (size, temperature and precipitation)
	s = (s - scaling_dbh_mu_M)/scaling_dbh_sd_M;
	temp = (temp - scaling_temp_mu_M)/scaling_temp_sd_M;
	precip = (precip - scaling_precip_mu_M)/scaling_precip_sd_M;

	// Define the coefficient of the polynom of s
	beta_0 = intercept_M +
	(beta_cs_M + beta_cs_T_M*temp + beta_cs_T_sq_M*temp*temp
		+ beta_cs_P_M*precip + beta_cs_P_sq_M*precip*precip)*cs +
	beta_T_M*temp + beta_T_sq_M*temp*temp +beta_P_M*precip + beta_P_sq_M*precip*precip;

	beta_1 = beta_dbh_M;

	beta_2 = beta_dbh_sq_M;

	// Polynom of s (order 2)
	double dbh_polynom = beta_0 + beta_1*s + beta_2*s*s;
	double dbh_polynom_differentiate = beta_1 + 2*beta_2*s;

	return 0; // dbh_polynom_differentiate / (2 + std::exp(-dbh_polynom) + std::exp(dbh_polynom));
	// return (0);
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

	os << "fecundity:" << std::endl << species.fecundity << std::endl;
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

#endif
