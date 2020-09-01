
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
	public:
		Patch(Environment const env, std::vector<Species*> speciesList);

	private:
		std::map <Species*, Population> m_pop_map; // One population (vector of cohorts) of a species present in Environment
		std::map <Species*, std::string> m_filenamePattern_map;
	
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
