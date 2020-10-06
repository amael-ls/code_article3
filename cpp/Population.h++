
/* Description of the class
The class population represents an ensemble of cohorts of the same species at a
given place. I organise it:
	- A vector of cohorts sorted by decreasing diameter (or height)
	- A species constant pointer
	- A height threshold height_star (and dbh_star for demography parameterised for dbh only)

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
#include <filesystem>
#include <execution>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
// #include <omp.h>

// My headers
#include "Error_classes.h++"
#include "Environment.h++"
#include "Species.h++"
#include "Cohort.h++"

class Population
{
	friend class Patch;

	public :
		// Constructors
		Population(unsigned int const maxCohorts, Species const * const species, std::string const summaryFilename, std::string const popDynFilename);
		Population(unsigned int const maxCohorts, Species const * const species, std::string const& initFilename,
			std::string const summaryFilename, std::string const popDynFilename);

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Population const &pop);
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
	// --- Functions
		void cohortDynamics(double const t, double const delta_t, double const height_star, Environment const & env); // Call euler and seedProduction
		void euler(double const t, double const delta_t, double const dbh_star, Environment const & env);
		void seedProduction(double const height_star); // Compute the local seed production (private function)
		void recruitment(double const t, double const delta_t, double const dbh_star, Environment const & env, bool& isPopulated); // Compute recruitment (dispersal + euler)
		void totalDensity_totalTrunkArea();
	
	// --- Sorting and organising
		void sort(bool const decreasingOrder);
		bool mergeCohorts(double const thresholdSimilarity, double const thresholdDensity);
		void resetCohorts(std::vector<Cohort>::iterator const it);

	// --- Reproduction
		unsigned int m_nonZeroCohort; // Number of non empty cohorts
		double m_localProducedSeeds; // Seeds produced in local patch (part of it are propagated)
		double m_localSeedBank; // Seeds received from local source and from external source (dispersal)
	
	// --- Total population
		double m_sumTrunkArea; // Basal area, an output variable
		double m_totalDensity; // Total density, an output variable: Integral[N(s, t) ds, from = 0, to = +Inf]

	// --- Others
		unsigned int m_currentIter;
		
	// Saving files
	// --- Ofstreams names
		std::string m_summary_ofs; // Summary file containing basal area, total density, and competition
		std::string m_popDyn_ofs; // Population dynamics file containing all the non zero cohorts

	// --- Output writing function
		void saveResults(double const height_star, bool const writeSummary);
		void summary(double const height_star);
		void closeOutputFiles();
		
};

#endif
