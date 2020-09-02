
#ifndef PATCH_C
#define PATCH_C

// My headers
#include "Patch.h++"

// Define typedef shortcuts
typedef std::map<Species*, Population>::iterator population_it;
typedef std::map<Species*, Population>::const_iterator c_population_it;

typedef std::map<Species*, std::string>::iterator filename_it;

Patch::Patch(Environment const& env, std::vector<Species*> speciesList, unsigned int const maxCohorts):
	m_maxCohorts(maxCohorts), m_env(env)
{
	std::vector<Species*>::const_iterator species_it = speciesList.cbegin();
	for (; species_it != speciesList.cend(); ++species_it)
	{
		m_filenamePattern_map[*species_it] = "not initialised.txt";
		m_pop_map.emplace(*species_it, Population(maxCohorts, *species_it));
	}
}

Patch::Patch(Environment const& env, Species* species, unsigned int const maxCohorts):
	m_maxCohorts(maxCohorts), m_env(env)
{
	m_pop_map.emplace(species, Population(maxCohorts, species));
}

Patch::Patch(Environment const& env, Species* species, std::string const initFilename, unsigned int const maxCohorts):
	m_maxCohorts(maxCohorts), m_env(env)
{
	m_pop_map.emplace(species, Population(maxCohorts, species, initFilename));
}


// population_it pop_it = m_pop_map.begin();
// filename_it file_it = m_filenamePattern_map.begin();





// /* Competition calculation:
// To calculate the competition, I start from the biggest tree, and calculate the
// its crown area. As long as I have not found s* such that the sum of the crown
// areas of trees above s* equals the ground area, I decrement s* of delta_s. If
// such an s* cannot be found, I set s* = 0. It means there is a gap in the canopy.

// Remarks:
// 1. This calculation is costly, and has to be done at each time step!
// 2. In the case of the flat-top model, there is no need to integrate over s, we
// can directly jump from cohorts to cohorts as the crown area is a step function
// of size s in this particular case. The flat-top case is also discontinuous, and
// there is no guarantee of having a solution for summedArea == 1.
// */
// void Population::competition()
// {
// 	double summedArea = 0;
// 	c_cohort_it it_cohort = m_cohortsVec.cbegin();
// 	m_s_star = it_cohort->m_mu;

// 	while ((m_s_star > 0) && (summedArea < 1)) // needs a plot_area later? Is it a density and thus 1?
// 	{
// 		// Reset the values
// 		summedArea = 0;
// 		it_cohort = m_cohortsVec.cbegin();
// 		while ((it_cohort->m_mu >= m_s_star) && (it_cohort != m_cohortsVec.cend()))
// 		{
// 			// crown area evaluated at s* given cohort of size m_mu, multiplied by density m_lambda of the cohort
// 			summedArea = it_cohort->crownArea(m_s_star)*it_cohort->m_lambda;
// 			++it_cohort;
// 		}
// 		m_s_star = m_s_star - m_delta_s;
// 	}
// 	if (m_s_star < 0)
// 		m_s_star = 0;
// }

// void Population::competition(double const t)
// {
// 	m_s_star = std::sqrt(1 - 2*m_species->d(0, 0, 0, 0)*t) - 1; // death rate is constant, so no problem
// }


#endif
