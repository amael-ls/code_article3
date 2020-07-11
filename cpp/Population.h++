
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
#include <vector>
#include <string>

// My headers
#include "Environment.h++"
#include "Species.h++"
#include "Cohort.h++"

class Population
{
	public :
		// Constructors
		Population(unsigned int const maxCohorts, Species* const sp, std::vector<double> const & lambda,
			std::vector<double> const & mu, Environment* const env, unsigned int currentIter);
		Population(unsigned int const maxCohorts, Species* const sp, std::vector<Cohort> const & cohorts,
			Environment* const env, unsigned int currentIter);
		Population(unsigned int const maxCohorts, Species* const sp, std::string const& fileName,
			Environment* const env, unsigned int currentIter);

		// Dynamics
		void euler(unsigned int n_t, double t0, double t_max,
			std::string const& outCompReprod = "compReprod.txt", std::string const& popTimeFile = "popDyn.txt");
		void rk4(unsigned int n_t, double t0, double t_max,
			std::string const& outCompReprod = "compReprod.txt", std::string const& popTimeFile = "popDyn.txt");
		double reproduction();
		void competition();
		void competition(double const t);
		void totalDensity_basalArea();

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Population const &pop);

		// Sorting and organising
		void sort(bool const decreasingOrder);
		// void lastReproducer();
		bool mergeCohorts(double const thresholdSimilarity, double const thresholdDensity);
		void resetCohorts(std::vector<Cohort>::iterator const it);
		void printNonZero() const;

	private :
		unsigned int m_maxCohorts; // maximal number of cohorts
		unsigned int m_nonZeroCohort; // Number of non empty cohorts
		// std::vector<Cohort>::iterator m_lastReproducer; // Last cohort able to reproduce (given sorted by decreasing order)
		unsigned int m_currentIter;
		double const m_s_inf, m_delta_s; // max possible size and size step for integration
		std::vector<Cohort> m_cohortsVec; // The population of cohorts
		Species* const m_species; // The pointer is constant, as a population shall not change species
		Environment* const m_env; // The pointer is constant, as a population shall not change environment (but Temp/precip in Env could be dynamic)
		double m_s_star; // Competition, both used in dynamics and output variable
		double m_basalArea; // Basal area, an output variable
		double m_totalDensity; // Total density, an output variable: Integral[N(s, t) ds, from = 0, to = +Inf]
		double m_localProducedSeeds; // Seeds produced in local patch (part of it are propagated)
		double m_localSeedBank; // Seeds received from local source and from external source (dispersal)
};

#endif
