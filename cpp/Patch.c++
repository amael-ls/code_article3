
#ifndef PATCH_C
#define PATCH_C

// My headers
#include "Patch.h++"

// Define typedef shortcuts
typedef std::map<Species*, Population>::iterator population_it;
typedef std::map<Species*, Population>::const_iterator c_population_it;

typedef std::map<Species*, std::string>::iterator filename_it;

Patch::Patch(Environment const& env, std::vector<Species*> const speciesList,
	std::string const initPath, std::string const initFilenamePattern,
	std::string const summaryPath, std::string const summaryFilenamePattern, 
	std::string const popDynPath, std::string const popDynFilenamePattern,
	unsigned int const maxCohorts):
		m_maxCohorts(maxCohorts), m_isPopulated(false), m_env(env), m_height_star(0)
{
	// Initialisation data
	std::string initFile;
	bool foundAnInitFile = false;

	// Species
	std::vector<std::string> speciesNames;
	std::vector<Species*>::const_iterator species_it = speciesList.cbegin();

	// Output data
	std::string summaryFile;
	std::string popDynFile;

	// Others
	m_minDelta_s = std::numeric_limits<double>::infinity();

	for (; species_it != speciesList.cend(); ++species_it)
	{
		speciesNames.emplace_back((*species_it)->m_speciesName);
		initFile = initPath + (*species_it)->m_speciesName + "/" + initFilenamePattern +
		std::to_string(env.m_patchId) + ".txt";

		summaryFile = summaryPath + (*species_it)->m_speciesName + "/" + summaryFilenamePattern +
		std::to_string(env.m_patchId) + ".txt";

		popDynFile = popDynPath + (*species_it)->m_speciesName + "/" + popDynFilenamePattern +
		std::to_string(env.m_patchId) + ".txt";
		
		if (std::filesystem::exists(initFile))
		{
			m_pop_map.emplace(*species_it, Population(maxCohorts, *species_it, initFile, summaryFile, popDynFile));
			m_filenamePattern_map[*species_it] = initFile;
			foundAnInitFile = true;
			m_isPopulated = true;
		}
		else
		{
			m_pop_map.emplace(*species_it, Population(maxCohorts, *species_it, summaryFile, popDynFile));
			m_filenamePattern_map[*species_it] = "not initialised";
		}

// Would be smarter to do it in forest rather than here, and to transmit it as an argument rather than having it as a member
		if (((m_pop_map.find(*species_it))->second).m_delta_s < m_minDelta_s)  
			m_minDelta_s = ((m_pop_map.find(*species_it))->second).m_delta_s;
	}

	if (m_env.m_initPopulated && !foundAnInitFile)
		throw Except_Patch(m_env.m_patchId, speciesNames);

	this->competition(m_minDelta_s);
}

/************************************/
/******        dynamics        ******/
/************************************/
// The following function just call competition and Euler for all species (i.e. populations of one patch)
void Patch::populationDynamics(double const t, double const delta_t)
{
	this->competition(m_minDelta_s); // Update m_height_star
	population_it pop_it = m_pop_map.begin();
	for (; pop_it != m_pop_map.end(); ++pop_it) // Go over all the species
		(pop_it->second).cohortDynamics(t, delta_t, m_height_star, m_env);
}

void Patch::dispersal(Patch* sourcePatch, Species* species,
	std::map<Distance, double> const& distToIntegral, double const deltaLat, double const deltaLon)
{
	Distance dist(m_env, sourcePatch->m_env, deltaLat, deltaLon);

	// withdrawn seeds from source = local production of source times amount dispersed by K
	double withdrawnSeeds = ((sourcePatch->m_pop_map).at(species)).m_localProducedSeeds * distToIntegral.at(dist);
	(m_pop_map.at(species)).m_localSeedBank += withdrawnSeeds;
	((sourcePatch->m_pop_map).at(species)).m_localProducedSeeds -= withdrawnSeeds;
}

void Patch::recruitment(Species* species, double const t, double const delta_t)
{
	double const dbh_star = std::exp(1/species->b*(std::log10(m_height_star) - species->a)*std::log(10));
	(m_pop_map.at(species)).recruitment(t, delta_t, dbh_star, m_env, m_isPopulated);
}

/* Competition calculation:
To compute competition within a patch, it is necessary to first get all
the cohorts (that have a non zeroheight), regardless species. Once the
cohorts are all gathered, I use dichotomy:
	- First set s* = maxHeight/2
	- Compute the sum of the crown areas for such s*
	- Change the boundaries of the interval in which looking for s*
		(if sumArea is to high, then the lower boundary)
	- Set s* halfway to the new boundaries

If the interval becomes to narrow (i.e., the difference of its upper and
lower boundaries is below heightTol), then I set s* = 0. It means there is
a gap in the canopy.

Remarks:
1. This calculation is costly, and has to be done at each time step!
2. In the case of the flat-top model, there is no need to integrate over s, we
can directly jump from cohorts to cohorts as the crown area is a step function
of size s in this particular case. The flat-top case is also discontinuous, and
there is no guarantee of having a solution for summedArea == 1.
*/
void Patch::getAllNonZeroCohorts(std::vector<Cohort *>& nonZeroCohorts) const
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

	// Sort vector
	std::sort(std::execution::par_unseq, nonZeroCohorts.begin(), nonZeroCohorts.end(), greaterCohortPtr);
}

void Patch::competition(double const heightTol) // The best heightTol is delta_s
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

		while (supBound - infBound > heightTol)
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

/**************************************/
/******        Organising        ******/
/**************************************/
void Patch::savePatch() const
{
	c_population_it pop = m_pop_map.cbegin();
	for (; pop != m_pop_map.cend(); ++ pop)
		(pop->second).savePopulation(m_height_star);
}

void Patch::summary() const
{
	c_population_it pop = m_pop_map.cbegin();
	for (; pop != m_pop_map.cend(); ++ pop)
		(pop->second).summary(m_height_star);
}

/************************************/
/******        Overload        ******/
/************************************/
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
