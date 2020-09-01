
/* Description of the class
The class population represents an ensemble of cohorts of the same species at a
given place. I organise it:
	- A vector of cohorts sorted by decreasing diameter (or height)
	- A species constant pointer
	- A height threshold s* (s_star, not to be mistaken with a pointer of type s)

It contains:
	- maximal number of cohorts
	- Number of non empty cohorts
	- Last cohort able to reproduce (given sorted by decreasing order)
	- The population of cohorts
	- The species
	- s*

I list the functions here, but describe them in the associated c++ file:
	- Two constructors, 1st with the composant of the cohorts, 2nd with a vector of cohorts
	- Euler, reproduction, and competition to solve the PDEs
	- Overloading for ease
	- Sorting functions to keep the cohorts organised:
		* Sort them by decreasing order
		* Find the last tree able to reproduce (only canopy trees, understorey trees do not reproduce)
		* Merge nearly equal cohorts when required (lack of space)
		* Reset cohorts
		* Print non zero cohorts
*/

#ifndef POPULATION_H
#define POPULATION_H

// Official headers
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// My headers
#include "Error_classes.h++"
#include "Environment.h++"
#include "Species.h++"
#include "Cohort.h++"

class Population
{
	friend class Forest;

	public :
		// Constructors
		Population(unsigned int const maxCohorts, double const s_inf, double const delta_s, Species const * const species);
		Population(unsigned int const maxCohorts, double const s_inf, double const delta_s, Species const * const species, std::string const& filename);

		// *********
		// Population(unsigned int const maxCohorts, Species* const sp, std::vector<double> const & lambda, std::vector<double> const & mu,
		// 	Environment* const env, unsigned int currentIter, std::string const compReprodFilename, std::string const popDynFilename);
		// Population(unsigned int const maxCohorts, Species* const sp, std::vector<Cohort> const & cohorts, Environment* const env,
		// 	unsigned int currentIter, std::string const compReprodFilename, std::string const popDynFilename);
		// Population(unsigned int const maxCohorts, Species* const sp, double const lambda, Environment* const env,
		// 	unsigned int currentIter, std::string const compReprodFilename, std::string const popDynFilename);
		// Population(unsigned int const maxCohorts, Species* const sp, std::string const& fileName, Environment* const env,
		// 	unsigned int currentIter, std::string const compReprodFilename, std::string const popDynFilename);

		// Dynamics

		// *********
		void euler(double const t, double const delta_t, double const s_star, Environment const & env);
		void recruitment(double const t, double const delta_t, double const s_star, Environment const & env);
		double reproduction(double const s_star);
		// void competition();
		// void competition(double const t);
		void totalDensity_basalArea();

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Population const &pop);

		// Sorting and organising
		void sort(bool const decreasingOrder);
		bool mergeCohorts(double const thresholdSimilarity, double const thresholdDensity);
		void resetCohorts(std::vector<Cohort>::iterator const it);
		void printNonZero() const;

	private :
	// Length vector cohort and discretisation
		unsigned int const m_maxCohorts; // Maximal number of cohorts
		double const m_s_inf; // Maximal possible size
		double const m_delta_s; // Size step for integration
	
	// Cohorts of species
		std::vector<Cohort> m_cohortsVec; // The population of cohorts
		Species const * const m_species; // The pointer is constant, as a population shall not change species

	// Dynamics
	// --- Reproduction
		unsigned int m_nonZeroCohort; // Number of non empty cohorts
		double m_localSeedBank; // Seeds received from local source and from external source (dispersal)
	
	// --- Total population
		double m_basalArea; // Basal area, an output variable
		double m_totalDensity; // Total density, an output variable: Integral[N(s, t) ds, from = 0, to = +Inf]

	// --- Others
		unsigned int m_currentIter;
		
		// double m_localProducedSeeds; // Seeds produced in local patch (part of it are propagated)
		// double m_s_star; // Competition, both used in dynamics and output variable
	
		
	// *********
		// Environment* const m_env; // The pointer is constant, as a population shall not change of location. However, the environment properties can change (only ptr is constant)
	// Saving files
		// std::ofstream m_compReprod_ofs; // ofstream that will be open at creation of Population. The non destruction of Pop can be a problem
		// std::ofstream m_popDyn_ofs; // ofstream that will be open at creation of Population. The non destruction of Pop can be a problem
		
};

#endif
