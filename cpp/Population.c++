
#ifndef POPULATION_C
#define POPULATION_C

// Official headers

// My headers
#include "Population.h++"

// Define typedef shortcuts
typedef std::vector<Cohort>::iterator cohort_it;
typedef std::vector<Cohort>::const_iterator c_cohort_it;

/****************************************/
/******        Constructors        ******/
/****************************************/
Population::Population(unsigned int const maxCohorts, Species const * const species, std::string const summaryFilename, std::string const popDynFilename):
	m_maxCohorts(maxCohorts), m_s_inf(species->maxDiameter), m_delta_s(m_s_inf/maxCohorts), m_species(species),
	m_nonZeroCohort(0), m_localProducedSeeds(0), m_localSeedBank(0), m_currentIter(0),
	m_summary_ofs(summaryFilename), m_popDyn_ofs(popDynFilename)
{
	// Assign vector of cohorts
	for (unsigned int i = 0; i < m_maxCohorts; ++i)
		m_cohortsVec.emplace_back(Cohort(species, 0));

	// Sort and compute total trunk area and total density
	this->sort(true); // true to sort by decreasing size
	this->totalDensity_totalTrunkArea();

	/***** Open ofstreams to save initial condition and close them (too many to be kept open) *****/
	// Open ofstreams m_summary_ofs
	// if (std::filesystem::exists(summaryFilename)) // Remove file if already exists
	// 	std::filesystem::remove(summaryFilename);
	
	std::ofstream summary(m_summary_ofs);

	if(!summary.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR (from constructor Population): cannot open output file <" << summaryFilename << ">"
			<< "> for species <" << species->m_speciesName << ">";
		throw (std::runtime_error (ss.str()));
	}
	
	summary << "iteration localSeedProduced localSeedBank sumTrunkArea totalDensity" << std::endl;
	summary << m_currentIter << " " << m_localProducedSeeds << " " << m_localSeedBank << " " <<
		m_sumTrunkArea << " " << m_totalDensity << std::endl;

	// Open ofstream m_popDyn_ofs
	if (std::filesystem::exists(popDynFilename)) // Remove file if already exists
		std::filesystem::remove(popDynFilename);
	
	std::ofstream popDyn(m_popDyn_ofs);

	if(!popDyn.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR (from constructor Population): cannot open output file <" << popDynFilename << ">"
			<< "> for species <" << species->m_speciesName << ">";
		throw (std::runtime_error (ss.str()));
	}

	popDyn << "iteration iterationBirth density dbh height" << std::endl;
	popDyn << *this;

	// Close ofstreams
	summary.close();
	popDyn.close();
}

Population::Population(unsigned int const maxCohorts, Species const * const species, std::string const& initFilename,
	std::string const summaryFilename, std::string const popDynFilename):
	m_maxCohorts(maxCohorts), m_s_inf(species->maxDiameter), m_delta_s(m_s_inf/maxCohorts), m_species(species),
	m_nonZeroCohort(0), m_localProducedSeeds(0), m_localSeedBank(0), m_currentIter(0),
	m_summary_ofs(summaryFilename), m_popDyn_ofs(popDynFilename)
{
	// Assign vector of cohorts
	// --- Open input file
	std::ifstream inputFile(initFilename);
	if(!inputFile.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR: cannot open file <" << initFilename << ">";
		throw (std::runtime_error (ss.str()));
	}

	std::string line;
	getline(inputFile, line);
	if (line.find("density", 0) == std::string::npos || line.find("dbh", 0) == std::string::npos)
	{
		std::stringstream ss;
		ss << "*** ERROR: file <" << initFilename << "> must have density and dbh headers";
		throw (std::runtime_error (ss.str()));
	}
	
	double density(0), dbh(0);
	size_t afterNumVal; // This value is set by std::stod to position of the next character in str after the numerical value

	while(inputFile.good())
	{
		getline(inputFile, line);
		if (line.empty())
			continue;
		density = std::stod(line, &afterNumVal);
		dbh = std::stod(line.substr(afterNumVal));

		if (density != 0)
		{
			m_cohortsVec.emplace_back(Cohort(density, dbh, species, 0));
			++m_nonZeroCohort;
		}

		if (m_maxCohorts < m_nonZeroCohort)
			throw(Except_Population(m_maxCohorts, initFilename));
	}

	inputFile.close();

	// Fill with zero cohorts up to m_maxCohorts. No problem if m_cohortsVec full
	for (int count = m_nonZeroCohort; count < m_maxCohorts; ++count)
		m_cohortsVec.emplace_back(Cohort(species, 0));
	
	double tallest_tree = std::max_element(m_cohortsVec.cbegin(), m_cohortsVec.cend())->m_mu;
	if (m_s_inf < tallest_tree)
		throw(Except_Population(m_s_inf, tallest_tree, initFilename));

	// Sort and compute total trunk area and total density
	this->sort(true); // true to sort by decreasing size
	this->totalDensity_totalTrunkArea();

/***** Open ofstreams to save initial condition and close them (too many to be kept open) *****/
	// Open ofstreams m_summary_ofs
	// if (std::filesystem::exists(summaryFilename)) // Remove file if already exists
	// 	std::filesystem::remove(summaryFilename);
	
	std::ofstream summary(m_summary_ofs);

	if(!summary.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR (from constructor Population): cannot open output file <" << summaryFilename
			<< "> for species <" << species->m_speciesName << ">";
		throw (std::runtime_error (ss.str()));
	}
	
	summary << "iteration localSeedProduced localSeedBank sumTrunkArea totalDensity" << std::endl;
	summary << m_currentIter << " " << m_localProducedSeeds << " " << m_localSeedBank << " " <<
		m_sumTrunkArea << " " << m_totalDensity << std::endl;

	// Open ofstream m_popDyn_ofs
	if (std::filesystem::exists(popDynFilename)) // Remove file if already exists
		std::filesystem::remove(popDynFilename);
	
	std::ofstream popDyn(m_popDyn_ofs);

	if(!popDyn.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR (from constructor Population): cannot open output file <" << popDynFilename << ">"
			<< "> for species <" << species->m_speciesName << ">";
		throw (std::runtime_error (ss.str()));
	}

	popDyn << "iteration iterationBirth density dbh height" << std::endl;
	popDyn << *this;

	// Close ofstreams
	summary.close();
	popDyn.close();
}

/********************************************/
/******        Euler & dynamics        ******/
/********************************************/
/* Euler (explicit) method:
To solve equation y'(t) = f(t, y). Iterative method with a time step delta_t
	y_{n + 1} = y_n + delta_t f(t_n, y_n)
For stability reason, I use Runge-Kutta 4 method later, but I left Euler methods

In this method, I had to manage the negligeable cohorts. Indeed, there is only a
finite number of cohorts, although along the characteristics, the population de-
crease exponentially (and therefore never reach 0). Hence, cohorts that are si-
milar are merged with their characteristics (density and diameter) averaged.

Once this is done, I integrate within the size space Omega (Ω), i.e., everywhere
but the boundary condition. This use the ODE II system of equations.

Lastly, I compute the dynamics at the boundary condition. This use the ODE V
system of equations.

Remarks:
1. m_cohortsVec.begin() + m_nonZeroCohort means that I will iterate on m_nonZeroCohort
values. Indeed, c++ starts iterations from 0. So when writing
	xxx != yyy.begin() + m_nonZeroCohort
it goes from 0 to m_nonZeroCohort - 1

2. See species.h++ for allometric functions and demographic functions

*/
void Population::cohortDynamics(double const t, double const delta_t, double const height_star, Environment const & env)
{
	// Convert height_star (h*) to species-specific dbh_star (demography is parameterised for dbh only)
	double const dbh_star = std::exp(1/m_species->b*(std::log10(height_star) - m_species->a)*std::log(10));

	// Compute local seed bank
	this->seedProduction(height_star);

	// Age cohort
	this->euler(t, delta_t, dbh_star, env);

	// Compute total trunk area, and total density
	this->totalDensity_totalTrunkArea();
}

void Population::euler(double const t, double const delta_t, double const dbh_star, Environment const & env)
{
	cohort_it it;
	cohort_it lim_it; // limit iterator

	// Check there is space to create a new cohort and merge/delete if required
	if (m_nonZeroCohort == m_maxCohorts)
		this->mergeCohorts(m_delta_s/8, 0.001);

	if (m_maxCohorts < m_nonZeroCohort)
		throw(Except_Population(m_maxCohorts, m_nonZeroCohort));
	
	// m_compReprod_ofs << t + delta_t << " "; // Added delta_t because, t is one step behind (Euler explicit)
	lim_it = m_cohortsVec.begin() + m_nonZeroCohort; // It might involve segmentation fault if maxCohort < nonZero

	// Integration within the size space Omega
	for (it = m_cohortsVec.begin(); it != lim_it; ++it)
		it->euler(t, delta_t, dbh_star, env, &Cohort::ODE_II);
}

void Population::recruitment(double const t, double const delta_t, double const dbh_star, Environment const & env, bool& isPopulated)
{
	cohort_it recruitment_it = m_cohortsVec.begin() + m_nonZeroCohort; // Position of recruitment (which becomes a new cohort if size > threshold)
	
	// Update iteration (recruitment is the last thing done)
	++m_currentIter;

	if (recruitment_it >= m_cohortsVec.end())
		throw(Except_Population(m_maxCohorts, m_nonZeroCohort, t));

	// Aging boundary cohort
	recruitment_it->euler(t, delta_t, dbh_star, env, m_localSeedBank, &Cohort::ODE_V);

	// If it reaches the threshold, the boundary cohort is released within Omega
	if (recruitment_it->m_mu > m_delta_s)
	{
		++m_nonZeroCohort;
		recruitment_it->m_birthIteration = m_currentIter;
		// Initialise properly the dbh (see equations 16 and 17 paper)
		recruitment_it->m_mu /= recruitment_it->m_lambda; // π = (η - a)λ Equation 16. Here, a = 0 and π and η have the same name: recruitment_it->m_mu
		
		if (!isPopulated) // If isPopulated is originally false, change it to true due to the newly created cohort
			isPopulated = true;
	}

	// Reset local seed bank to 0 (they have been transferred to the boundary condition by euler and ODE_V)
	m_localSeedBank = 0;
}

/* Seed production:
Only canopy trees can reproduce, hence we stop at the last tree that can
reproduce (i.e., the last tree taller than s*).
In Purves2008, the reproduction is proportional to the total sun-exposed crown
area (which is 0 for understorey trees). That is to say:
	Seed production = fecundity x crownArea(s, s*)
The fecundity parameter is from table S4 Purves2008. It is a difficult parameter
to estimate (check his Appendix2: Parameter estimation)
*/
void Population::seedProduction(double const height_star)
{
	double popReprod = 0;
	c_cohort_it cohort_it = m_cohortsVec.cbegin();

	while (cohort_it != m_cohortsVec.end() && cohort_it->m_height > height_star) // sum_l F * lambda (eq 27 article), sum_k is managed by Forest
	{
		popReprod += cohort_it->crownArea(height_star) * cohort_it->m_lambda;
		++cohort_it;
	}

	m_localProducedSeeds = popReprod * m_species->fecundity;
}

// ****************************************************************************************************************
/* To compute the total trunk area of the community and the 'number of individuals'
The total trunk area is expresses in metres. I, therefore, coerce the diameter
to metres
The number of individuals is supposed to be (spatial case):
	\int_{x \in Plot} \int_{0}^{+Inf} N(s, x, t) ds dx
In my case, the density is per m per square meter (the first meter is for the
size s, the square meter is for the ground area). I therefore conclude:
	\int_{x \in Plot} \int_{0}^{+Inf} N(s, x, t) ds dx
		= \int_{x \in Plot} 1 dx \times \int_{0}^{+Inf} N(s, x, t) ds
		= plotArea \times \int_{0}^{+Inf} N(s, x, t) ds
*/
void Population::totalDensity_totalTrunkArea()
{
	double sumTrunkArea = 0;
	double currentDensity = 0;

	c_cohort_it it = m_cohortsVec.cbegin();
	c_cohort_it it_lim = m_cohortsVec.cbegin() + m_nonZeroCohort;
	for (; it != it_lim; ++it)
	{
		sumTrunkArea += (it->m_mu)*(it->m_mu) * (it->m_lambda); // dbh² \times density per m per m². The += is the 'integral over s'
		currentDensity += it->m_lambda; //
	}

	sumTrunkArea *= M_PI/(4*1000); // π dbh^2/(4 x 1000), the 1000 is to convert from mm to metres
	m_sumTrunkArea = sumTrunkArea;
	m_totalDensity = currentDensity;
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream& operator<<(std::ostream& os, Population const &pop)
{
	c_cohort_it lim_it = pop.m_cohortsVec.cbegin() + pop.m_nonZeroCohort;
	for (c_cohort_it it = pop.m_cohortsVec.cbegin(); it != lim_it; ++it)
		os << pop.m_currentIter << " " << *it << std::endl;
	return os;
}

/************************************************/
/******        Sorting & organising        ******/
/************************************************/
void Population::sort(bool const decreasingOrder)
{
	cohort_it first = m_cohortsVec.begin();
	cohort_it last = m_cohortsVec.end();
	if (decreasingOrder)
		std::sort(first, last, std::greater<Cohort>());
	else
		std::sort(first, last);
}

/* Merge/delete cohorts
There is only a limited amount of cohorts, set by the user (m_maxCohorts in the
population's constructor). Therefore, to avoid segmentation fault, it is neces-
sary to release some space. There are two ways:
	- Merge similar cohorts
	- Remove negligeable cohorts
*/
bool Population::mergeCohorts(double const thresholdSimilarity, double const thresholdDensity)
{
	cohort_it it = m_cohortsVec.begin() + 1; // moving iterator
	cohort_it first = m_cohortsVec.begin();
	cohort_it ref_it; // reference iterator
	cohort_it lim_it = first + m_nonZeroCohort; // limit iterator

	bool mergedCohorts = false;
	double mu, lambda;

	// First, merging cohorts taller than m_s_inf
	while ((it != m_cohortsVec.end()) && (it->m_mu > m_s_inf))
	{
		// Weighted (by density m_lambda) sum for mu
		first->m_mu = (first->m_lambda * first->m_mu + it->m_lambda * it->m_mu)/(first->m_lambda + it->m_lambda);
		first->m_lambda += it->m_lambda;
		this->resetCohorts(it);
		mergedCohorts = true;
		++it;
	}

	ref_it = it;

	// Second (if required), merge similar cohorts
	while (ref_it != lim_it)
	{
		mu = 0;
		lambda = 0;
		it = ref_it + 1;
		// abs val for float = fabs, useless in the next while loop because population is sorted (decreasing order)
		while ((it != lim_it) && (ref_it->m_mu - it->m_mu < thresholdSimilarity)) // BUG <--- TO CHECK, DID I FORGET TO REMOVE IT
		{
			mu += it->m_lambda * it->m_mu;
			lambda += it->m_lambda;
			this->resetCohorts(it);
			mergedCohorts = true;
			++it;
		} // END BUG <--- TO CHECK, DID I FORGET TO REMOVE IT
		ref_it->m_mu = (ref_it->m_lambda * ref_it->m_mu + mu)/(ref_it->m_lambda + lambda);
		ref_it->m_lambda += lambda;
		ref_it = it;
	}

	// Reseting cohorts considered negligible
	for (it = m_cohortsVec.begin(); it != m_cohortsVec.end(); ++it)
	{
		if ((0 < it->m_lambda) && (it->m_lambda < thresholdDensity))
			this->resetCohorts(it);
	}
	// Reordering, to put the zero cohorts at the end
	this->sort(true);
	std::cout << "mergeCohorts ended" << std::endl;
	return mergedCohorts;
}

void Population::resetCohorts(cohort_it const it)
{
	it->m_mu = 0;
	it->m_lambda = 0;
	if (m_nonZeroCohort > 0)
		--m_nonZeroCohort;
}

void Population::printNonZero() const
{
	c_cohort_it it = m_cohortsVec.cbegin();
	c_cohort_it lim_it = m_cohortsVec.cbegin() + m_nonZeroCohort;
	
	std::cout << "Species: " << m_species->m_speciesName << std::endl;
	for (; it != lim_it; ++it)
		std::cout << *it << std::endl;
	
	std::cout << "Number of zeros and boundary cohorts: " << m_cohortsVec.size() - m_nonZeroCohort << std::endl;
}

void Population::saveResults()
{
	std::ofstream summary(m_summary_ofs, std::ofstream::app);

	if(!summary.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR (from Population::saveResults): cannot open output file <" << m_summary_ofs << ">";
		throw (std::runtime_error (ss.str()));
	}

	std::ofstream popDyn(m_popDyn_ofs, std::ofstream::app);
	if(!popDyn.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR (from Population::saveResults): cannot open output file <" << m_popDyn_ofs << ">";
		throw (std::runtime_error (ss.str()));
	}

	summary << m_currentIter << " " << m_localProducedSeeds << " " << m_localSeedBank << " " <<
		m_sumTrunkArea << " " << m_totalDensity << std::endl;

	popDyn << *this;

	// Close ofstreams
	summary.close();
	popDyn.close();
}

#endif
