
#ifndef PATCH_C
#define PATCH_C

// My headers
#include "Patch.h++"

// Define typedef shortcuts
typedef std::map<Species*, Population>::iterator population_it;
typedef std::map<Species*, Population>::const_iterator c_population_it;

typedef std::map<Species*, std::string>::iterator filename_it;

Patch::Patch(Environment const& env, std::vector<Species*> speciesList, unsigned int const maxCohorts):
	m_maxCohorts(maxCohorts), m_env(env), m_height_star(0)
{
	std::vector<Species*>::const_iterator species_it = speciesList.cbegin();
	for (; species_it != speciesList.cend(); ++species_it)
	{
		m_filenamePattern_map[*species_it] = "not initialised.txt";
		m_pop_map.emplace(*species_it, Population(maxCohorts, *species_it));
	}
}

Patch::Patch(Environment const& env, Species* species, unsigned int const maxCohorts):
	m_maxCohorts(maxCohorts), m_env(env), m_height_star(0)
{
	m_pop_map.emplace(species, Population(maxCohorts, species));
}

Patch::Patch(Environment const& env, Species* species, std::string const initFilename, unsigned int const maxCohorts):
	m_maxCohorts(maxCohorts), m_env(env), m_height_star(0)
{
	m_pop_map.emplace(species, Population(maxCohorts, species, initFilename));
}


// population_it pop_it = m_pop_map.begin();
// filename_it file_it = m_filenamePattern_map.begin();





/* Competition calculation:
To calculate the competition, I start from the biggest tree, and calculate its
crown area. As long as I have not found s* such that the sum of the crown areas
of trees above s* equals the ground area, I decrement s* of delta_s. If such an
s* cannot be found, I set s* = 0. It means there is a gap in the canopy.

Remarks:
1. This calculation is costly, and has to be done at each time step!
2. In the case of the flat-top model, there is no need to integrate over s, we
can directly jump from cohorts to cohorts as the crown area is a step function
of size s in this particular case. The flat-top case is also discontinuous, and
there is no guarantee of having a solution for summedArea == 1.
*/
void Patch::getAllNonZeroCohorts(std::vector<Cohort *> nonZeroCohorts) const
{
	c_population_it it_map = m_pop_map.cbegin();
	std::vector<Cohort>::const_iterator it_cohort, lim_it;
	Cohort* cohort_ptr;

	for (; it_map != m_pop_map.cend(); ++it_map)
	{
		it_cohort = ((it_map->second).m_cohortsVec).cbegin();
		/*
			Problems on lim_it:
				- If there is no cohort?
				- Am I going one step too far? Segmentation fault, random cohort assignation value
		*/
		lim_it = it_cohort + (it_map->second).m_nonZeroCohort; // To avoid taking the zeros
		for (; it_cohort != lim_it; ++it_cohort)
		{
			cohort_ptr = new Cohort(*it_cohort);
			nonZeroCohorts.emplace_back(cohort_ptr);
		}
	}

	delete cohort_ptr;

	// Sort vector
	std::sort(nonZeroCohorts.begin(), nonZeroCohorts.end(), greaterCohortPtr);

	std::vector<Cohort *>::iterator it = nonZeroCohorts.begin();
	for (; it != nonZeroCohorts.end(); ++it)
		std::cout << **it << std::endl;
}

void Patch::competition(double const tolHeight) // The best tolHeight is delta_s
{
	// Get all non zero cohorts
	std::vector<Cohort *> nonZeroCohorts;

	this->getAllNonZeroCohorts(nonZeroCohorts);

	if (nonZeroCohorts.size() == 0)
		m_height_star = 0;
	else
	{
		double supBound = nonZeroCohorts[0]->m_height; // Vector sorted by decreasing height, so the max is at 0
		double infBound = 0;
		double sumArea = 0;

		unsigned int index = 0;

		while (supBound - infBound > tolHeight)
		{
			m_height_star = (supBound + infBound)/2;

			while ((index < nonZeroCohorts.size()) && nonZeroCohorts[index]->m_height >= m_height_star)
			{
				sumArea += nonZeroCohorts[index]->crownArea(m_height_star) * nonZeroCohorts[index]->m_lambda;
				++index;
			}

			if (sumArea < m_env.plotArea)
				supBound = m_height_star;
			if (sumArea >= m_env.plotArea)
				infBound = m_height_star;

			sumArea = 0;
			index = 0;
		}
		if (index >= nonZeroCohorts.size())
			m_height_star= 0;
	}
}

// /************************************/
// /******        Overload        ******/
// /************************************/
std::ostream& operator<<(std::ostream& os, Patch const &patch)
{
	c_population_it pop_it = patch.m_pop_map.cbegin();
	pop_it->first->printName(os); os << std::endl;
	for (; pop_it != patch.m_pop_map.cend(); ++pop_it)
	{
		(patch.m_env).printId(os); os << std::endl;
		os << pop_it->second << std::endl;
	}

	return os;
}

bool operator<(Patch const& patch1, Patch const& patch2)
{
	return (patch1.m_env < patch2.m_env);
}

bool operator>(Patch const& patch1, Patch const& patch2)
{
	return (patch1.m_env > patch2.m_env);
}

#endif
