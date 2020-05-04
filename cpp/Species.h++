
/* Description of the class
The class species defines a tree species by its parameters and vital rates.

It contains:
	- all the species specific parameters required to define a species (including allometries)
	- vital rates (functions)

I list the functions here, but describe them in the associated c++ file:
	- Constructor
	- Overloading for ease
	- Vital rates and their differentiate:
		* v (growth speed)
		* d (deat rate)
		* dv_ds differentiate of v with respect to s (size variable)
		* dd_ds differentiate of d with respect to s (size variable)
*/

#ifndef SPECIES_H
#define SPECIES_H

// Official headers
#include <string>
#include <cmath>

/*************************************/
/******      CLASS SPECIES      ******/
/*************************************/

class Species
{
	friend class Cohort;
	friend class Population;
	public :
		// constructors
		// Species();
		Species(const char* species_name, const std::string& delim);

		// friend functions and overload
		// Species& operator=(const Species& species);
		friend std::ostream& operator<<(std::ostream& os, const Species& species);
		friend bool operator<(Species const& species1, Species const& species2);
		// void printSpecies(std::ostream &os) const;

		// Demographic functions
		double v(double s, double const s_star, double temp, double precip) const;
		double d(double s, double const s_star, double temp, double precip) const;
		double dv_ds(double s, double const s_star, double temp, double precip) const;
		double dd_ds(double s, double const s_star, double temp, double precip) const;

		// others
		// bool compareName(Species* species_1, Species* species_2) const;

	private :
		// Species' name
		std::string m_speciesName;

		// Allometry parameters (height from dbh)
		double a, b;

		// Allometry parameters (crown area from dbh)
		double T_param, R0_C0, R0_C1, R40_C0, R40_C1, M_C0, M_C1, B_C0, B_C1;

		// Growth parameters (20)
		double intercept_G; // intercept
		double beta_dbh, beta_dbh_T, beta_dbh_T_sq, beta_dbh_P, beta_dbh_P_sq; // dbh (and climate interactions)
		double beta_dbh_sq, beta_dbh_sq_T, beta_dbh_sq_T_sq, beta_dbh_sq_P, beta_dbh_sq_P_sq; // dbh² (and climate interactions)
		double beta_cs, beta_cs_T, beta_cs_T_sq, beta_cs_P, beta_cs_P_sq; // competition (and climate interactions)
		double beta_T, beta_T_sq, beta_P, beta_P_sq; // climate

		// Mortality parameters (12)
		double intercept_M;
		double beta_dbh_M, beta_dbh_sq_M; // dbh and dbh²
		double beta_cs_M, beta_cs_T_M, beta_cs_T_sq_M, beta_cs_P_M, beta_cs_P_sq_M; // competition (and climate interactions)
		double beta_T_M, beta_T_sq_M, beta_P_M, beta_P_sq_M; // climate

		// Fecundity
		double fecundity;

		// Scaling (growth)
		double scaling_G_mu, scaling_G_sd;
		double scaling_dbh_mu_G, scaling_dbh_sd_G;
		double scaling_temp_mu_G, scaling_temp_sd_G, scaling_precip_mu_G, scaling_precip_sd_G;

		// Scaling (mortality)
		double scaling_dbh_mu_M, scaling_dbh_sd_M;
		double scaling_temp_mu_M, scaling_temp_sd_M, scaling_precip_mu_M, scaling_precip_sd_M;

		// Others
		double maxDiameter;
};

#endif
