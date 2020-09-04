
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
	// Constructor
		Patch(Environment const& env, std::vector<Species*> const speciesList, std::string const initPath,
			std::string const initFilenamePattern, unsigned int const maxCohorts);

	// Add population
		void addPopulation();

	// Dynamics
		void populationDynamics(double const t, double const delta_t);
		void competition(double const tolHeight);

	// Overloading
		friend std::ostream& operator<<(std::ostream& os, Patch const &patch);
		friend bool operator<(Patch const& patch1, Patch const& patch2);
		friend bool operator>(Patch const& patch1, Patch const& patch2);

	private:
	// Utilities
		std::map <Species*, std::string> m_filenamePattern_map;
		unsigned int m_maxCohorts;
		double m_minDelta_s;
		bool m_isPopulated;

	// Populations
		std::map <Species*, Population> m_pop_map; // One population (vector of cohorts) of a species present in Environment
	
	// Localisation
		Environment m_env; // Contains climate, coordinates and patch id

	// Dynamics
	// --- Competition
		double m_height_star; //The height at which the canopy is closed

	// --- Reproduction
		std::map <Species*, double> m_localSeedBank_map; // Might be useless, it is stored in each population anyway
		std::map <Species*, double> m_localSeedProduced_map; // Might be useless, it is stored in each population anyway

	// Private functions
	void getAllNonZeroCohorts(std::vector<Cohort *> nonZeroCohorts) const;
};

#endif
