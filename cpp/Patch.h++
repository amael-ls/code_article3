
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
		Patch(Environment const& env, std::vector<Species*> const speciesList,
			std::string const initPath, std::string const initFilenamePattern,
			std::string const summaryPath, std::string const summaryFilenamePattern, 
			std::string const popDynPath, std::string const popDynFilenamePattern,
			unsigned int const maxCohorts);

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

	// --- Functions
		void populationDynamics(double const t, double const delta_t);
		void dispersal(std::vector<Patch>::iterator targetPatch, Patch* sourcePatch, Species* species);
		void recruitment(std::vector<Patch>::iterator targetPatch, Species* species, double const t, double const delta_t);
		void competition(double const tolHeight);

	// Output writing function
		void saveResults();
		void closeOutputFiles();

	// Private functions
	void getAllNonZeroCohorts(std::vector<Cohort *> nonZeroCohorts) const;
};

#endif
