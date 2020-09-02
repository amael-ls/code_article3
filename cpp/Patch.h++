
/*
	This class handles ...
*/

#include <vector>
#include <string>
#include <map>

#include "Population.h++"
#include "Species.h++"

#ifndef PATCH_H
#define PATCH_H

class Patch
{
	friend class Forest;
	public:
	// Constructors
		Patch(Environment const& env, std::vector<Species*> speciesList, unsigned int const maxCohorts);
		Patch(Environment const& env, Species* species, unsigned int const maxCohorts);
		Patch(Environment const& env, Species* species, std::string const initFilename, unsigned int const maxCohorts);

	// Dynamics
		std::vector<Cohort*> getNonZeroCohorts() const;
		void competition(double const tolHeight);

	// Overloading
		friend std::ostream& operator<<(std::ostream& os, Patch const &patch);
		friend bool operator<(Patch const& patch1, Patch const& patch2);
		friend bool operator>(Patch const& patch1, Patch const& patch2);

	private:
	// Utilities
		std::map <Species*, std::string> m_filenamePattern_map;
		unsigned int m_maxCohorts;

	// Populations
		std::map <Species*, Population> m_pop_map; // One population (vector of cohorts) of a species present in Environment
	
	// Localisation
		Environment m_env; // Contains climate, coordinates and patch id

	// Dynamics
	// --- Competition
		double m_s_star;

	// --- Reproduction
		std::map <Species*, double> m_localSeedBank_map;
		std::map <Species*, double> m_localSeedProduced_map;
};

#endif
