
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

Allometries:
1. Calculate the height of a tree given its dbh:
	- Parameters for 106 species, source table S3, Purves2007, a and b
	- Function define at the end of the program, source appendix S2, end 1st page

2. Calculate the radius of the crown at a distance y from the top, given dbh of a tree:
	- Parameters for 106 species, source table S3, Purves2007, T
	- Function define ate the end of the program, source appendix S1, eq S1.6, S1.7
	
3. Data C0-C1:
	- Source = table S2 Purves2007

4. Bibliography:
		- Crown plasticity and competition for canopy space: a new spatially implicit model parameterized for 250 North American tree species
*/

#ifndef SPECIES_H
#define SPECIES_H

// Official headers
#include <iomanip> // std::setw, std::left, std::setprecision
#include <string>
#include <cmath> // for log

// Alglib header
#include "diffequations.h" // To compute minSizeReproduction

// My headers
#include "Error_classes.h++"
#include "Environment.h++"

/*************************************/
/******      CLASS SPECIES      ******/
/*************************************/

class Species
{
	// Friendship
	friend class Population;
	friend class Dispersal;
	friend class Forest;
	friend class Cohort;
	friend class Patch;

	public :
		// constructors
		Species(std::string const& species_filename, std::string const& species_path, const std::string& delim);

		// friend functions and overload
		friend std::ostream& operator<<(std::ostream& os, const Species& species);
		friend bool operator<(Species const& species1, Species const& species2);

		// Demographic functions
		double v(double s, double const s_star, double temp, double precip) const;
		double d(double s, double const s_star, double temp, double precip) const;
		double dv_ds(double s, double const s_star, double temp, double precip) const;
		double dd_ds(double s, double const s_star, double temp, double precip) const;

		// Dispersal
		double K(double const distance) const;
		double K(double delta_lon, double delta_lat) const;
		double K(double const longitude1, double const latitude1, double const longitude2, double const latitude2) const;
		bool isDirac() const;

		// others
		static void growth_callback(const alglib::real_1d_array &y, double x, alglib::real_1d_array &dy, void *ptr);
		void printName(std::ostream& os) const;
		std::string getName() const;

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

		// Mortality parameters (12 + 1)
		double intercept_M;
		double beta_dbh_M, beta_dbh_sq_M; // dbh and dbh²
		double beta_cs_M, beta_cs_T_M, beta_cs_T_sq_M, beta_cs_P_M, beta_cs_P_sq_M; // competition (and climate interactions)
		double beta_T_M, beta_T_sq_M, beta_P_M, beta_P_sq_M; // climate

		// Fecundity
		double fecundity;
		double minAgeReproduction;
		double minDbhReproduction;
		double minHeightReproduction;

		// Scaling (growth)
		double scaling_G_mu, scaling_G_sd;
		double scaling_dbh_mu_G, scaling_dbh_sd_G;
		double scaling_temp_mu_G, scaling_temp_sd_G, scaling_precip_mu_G, scaling_precip_sd_G;

		// Scaling (mortality)
		double scaling_dbh_mu_M, scaling_dbh_sd_M;
		double scaling_temp_mu_M, scaling_temp_sd_M, scaling_precip_mu_M, scaling_precip_sd_M;

		// Dispersal parameters, not necessarily all provided by the user
		// --- General parameters
		std::string keysToRead;
		double dispersalProbaThreshold; // Value between 0 and 1. Threshold from which probability K(x, y) is considered 0
		bool min_dispersalProba; // True if dispersalProbaThreshold is defined
		std::string refKernel_doi; // The DOI from which the kernel is from
		double dispersalDistThreshold; // Threshold distance from which probability K(x, y) is considered 0 (i.e., max dist)
		bool max_dispersalDist; // True if dispersalDistThreshold is defined

		// --- Moorcroft 2006
		double relLDDtoSDD; // Average long-distance dispersal relative to short-distance dispersal ---> 1/beta_l = 10 => beta_l = 0.1
		double propLDD; // Fraction of fecundity going to long-distance dispersal ---> p in Moorcroft, 0.1

		// Clark 1999, cf Boisvert-Marsh 2020 for some values
		double twoDt_a;
		double twoDt_b;

		// Gaussian Kernel, Cousens 2008 - Dispersal in plants, p. 82
		double gaussParam;

		// Laplace Kernel, Cousens 2008 - Dispersal in plants, p. 82
		double laplaceParam;

		// Others
		double maxDiameter;
};

#endif
