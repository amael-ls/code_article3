
// Official headers
#include <algorithm> // std::sort
#include <exception> // std::throw
#include <stdexcept> // exceptions (bad_alloc, bad_cast, out_of_range, ...)
#include <iomanip> // std::setw, std::left
#include <fstream>
#include <sstream>
#include <string>
// #include <chrono> // std::chrono::system_clock::now
#include <cmath>

// My headers
#include "Error_classes.h++"
#include "Population.h++"

// Define typedef shortcuts
typedef std::vector<Cohort>::iterator cohort_it;
typedef std::vector<Cohort>::const_iterator c_cohort_it;

/****************************************/
/******        Constructors        ******/
/****************************************/
/*
	The constructors are build in this order:
	First constructor: set λ and μ to values provided separately by the user
		λ, μ
		calculate s_star
		set m_crown_area, which is the area of the crown cross-section at s* for each cohort

	Second constructor: set λ and μ to values provided by the user with a vector of cohorts

*/
Population::Population(unsigned int const maxCohorts, Species* const sp,
	std::vector<double> const & lambda, std::vector<double> const & mu, Environment* const env):
	m_maxCohorts(maxCohorts), m_nonZeroCohort(lambda.size()), m_s_inf(sp->maxDiameter), m_delta_s(m_s_inf/maxCohorts),
	m_cohortsVec(maxCohorts), m_species(sp), m_env(env)
{
	if (m_maxCohorts < m_nonZeroCohort)
	{
		std::stringstream ss;
		ss << "Error (from constructor): dimensions mismatch. maxCohort = " << m_maxCohorts;
		ss << " nonZeroCohort = " << m_nonZeroCohort << std::endl;
		throw(std::out_of_range (ss.str()));
	}

	try
	{
		for (unsigned int i = 0; i < m_nonZeroCohort; ++i)
		{
			m_cohortsVec[i].m_lambda = lambda[i];
			m_cohortsVec[i].m_mu = mu[i];
			m_cohortsVec[i].m_species = m_species;
		}

		for (unsigned int i = m_nonZeroCohort; i < m_maxCohorts; ++i)
		{
			m_cohortsVec[i].m_lambda = 0;
			m_cohortsVec[i].m_mu = 0;
			m_cohortsVec[i].m_species = m_species;
		}
	}
	catch(const std::out_of_range& ex)
	{
		std::stringstream ss;
		ss << "Error (from Constructor): maximum population size exceded" << ex.what();
		throw (std::runtime_error (ss.str()));
	}

	double tallest_tree = std::max_element(m_cohortsVec.cbegin(), m_cohortsVec.cend())->m_mu;
	if (m_s_inf < tallest_tree)
	{
		std::stringstream ss;
		ss << "Error (from constructor): s_inf = " << m_s_inf << " must be larger than the tallest tree (" << tallest_tree << ")";
		throw(std::out_of_range (ss.str()));
	}

	this->sort(true); // true to sort by decreasing size
	this->competition();
}

Population::Population(unsigned int const maxCohorts, Species* const sp,
	std::vector<Cohort> const & cohorts, Environment* const env):
	m_maxCohorts(maxCohorts), m_nonZeroCohort(cohorts.size()), m_s_inf(sp->maxDiameter), m_delta_s(m_s_inf/maxCohorts),
	m_cohortsVec(cohorts), m_species(sp), m_env(env)
{
	double tallest_tree = std::max_element(m_cohortsVec.cbegin(), m_cohortsVec.cend())->m_mu;
	if (m_s_inf < tallest_tree)
	{
		std::stringstream ss;
		ss << "Error (from constructor): s_inf = " << m_s_inf << " must be larger than the tallest tree (" << tallest_tree << ")";
		throw(std::out_of_range (ss.str()));
	}

	if (m_maxCohorts < m_nonZeroCohort)
	{
		std::stringstream ss;
		ss << "Error (from constructor): maxCohorts = " << m_maxCohorts << " must be larger than";
		ss << " m_nonZeroCohort = " << m_nonZeroCohort << std::endl;
		throw(std::out_of_range (ss.str()));
	}

	 // Set species to each cohort
	for (cohort_it it = m_cohortsVec.begin(); it != m_cohortsVec.end(); ++it)
		it->m_species = m_species;

	this->sort(true); // true to sort by decreasing size
	this->competition();
}

Population::Population(unsigned int const maxCohorts, Species* const sp,
	std::string const& fileName, Environment* const env):
	m_maxCohorts(maxCohorts), m_s_inf(sp->maxDiameter), m_delta_s(m_s_inf/maxCohorts),
	m_species(sp), m_env(env)
{
	// std::cout << "From constr" << std::endl;
	// std::cout << *env << std::endl;

	std::ifstream inputFile(fileName);
	if(!inputFile.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR: cannot open file <" << fileName << ">";
		throw (std::runtime_error (ss.str()));
	}

	std::string line;
	getline(inputFile, line);
	if (line.find("density", 0) == std::string::npos || line.find("dbh", 0) == std::string::npos)
	{
		std::stringstream ss;
		ss << "*** ERROR: file <" << fileName << "> must have density and dbh headers";
		throw (std::runtime_error (ss.str()));
	}
	
	double density(0), dbh(0);
	size_t afterNumVal; // This value is set by std::stod to position of the next character in str after the numerical value
	m_nonZeroCohort = 0;

	while(inputFile.good())
	{
		getline(inputFile, line);
		if (line.empty())
			continue;
		density = std::stod(line, &afterNumVal);
		dbh = std::stod(line.substr(afterNumVal));

		m_cohortsVec.emplace_back(Cohort(density, dbh, sp));
		++m_nonZeroCohort;
		if (m_maxCohorts < m_nonZeroCohort)
			throw(Except_Population(m_maxCohorts, fileName));
	}

	std::cout << "non zero cohorts = " << m_nonZeroCohort << std::endl;
	std::cout << "maximum cohorts allowed = " << m_maxCohorts << std::endl;

	// Fill with zero cohorts up to m_maxCohorts. No problem if m_cohortsVec full
	std::cout << "Size cohorts = " << m_cohortsVec.size() << std::endl;
	for (int count = m_nonZeroCohort; count < m_maxCohorts; ++count)
		m_cohortsVec.emplace_back(Cohort(sp));
	
	double tallest_tree = std::max_element(m_cohortsVec.cbegin(), m_cohortsVec.cend())->m_mu;
	if (m_s_inf < tallest_tree)
		throw(Except_Population(m_s_inf, tallest_tree, fileName));

	this->sort(true); // true to sort by decreasing size
	this->competition();
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
*/
void Population::euler(unsigned int n_t, double t0, double t_max, std::string const& compReprodFile, std::string const& popTimeFile)
{
	double t;
	double delta_t = (t_max - t0)/(n_t - 1);
	double popReprod = 0;
	cohort_it it;
	cohort_it lim_it; // limit iterator

	std::ofstream outputCompReprod (compReprodFile);
	std::ofstream outputPopTime (popTimeFile);
	if(!outputCompReprod.is_open() || !outputPopTime.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR (from Population::euler): cannot open file";
		throw (std::runtime_error (ss.str()));
	}

	outputCompReprod << "time reproduction competition basalArea totalDensity" << std::endl;
	outputPopTime << "density dbh" << std::endl;

	std::cout << "Pop::euler starting" << std::endl;
    for (unsigned int i = 0; i < n_t; ++i) // time loop
	{
		// Check there is space to create a new cohort and merge/delete if required
		if (m_nonZeroCohort == m_maxCohorts)
			this->mergeCohorts(m_delta_s/8, 0.001);

		if (m_maxCohorts < m_nonZeroCohort)
		{
			std::stringstream ss;
			ss << "Error (from Euler): Index = " << m_nonZeroCohort << " out of boundaries. maxCohort = " << m_maxCohorts << std::endl;
			throw(std::out_of_range (ss.str()));
		}

		t = t0 + i*delta_t; // i starts at 0, hence it is explicit Euler
		outputCompReprod << t + delta_t << " "; // I write t_{n + 1}, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
		lim_it = m_cohortsVec.begin() + m_nonZeroCohort; // It might involve segmentation fault if maxCohort < nonZero

		// Integration within the size space Omega
		for (it = m_cohortsVec.begin(); it != lim_it; ++it)
			it->euler(t, delta_t, m_s_star, *m_env, &Cohort::ODE_II);

		// Dynamics at the boundary condition in size
		popReprod = this->reproduction();
		outputCompReprod << popReprod << " ";

		// New cohort of species m_species, lambda = mu = 0
		it->euler(t, delta_t, m_s_star, *m_env, popReprod, &Cohort::ODE_V);
		if (it->m_mu > m_delta_s) // If it reach the threshold, it becomes a cohort within Omega
			m_nonZeroCohort += 1;

		// std::cout << "m_nonZeroCohort = " << m_nonZeroCohort << std::endl;
		// Compute competition, basal area, and total density
		this->competition();
		this->totalDensity_basalArea();
		outputCompReprod << m_s_star << " " << m_basalArea << " " << m_totalDensity << std::endl;

		outputPopTime << *this;
    }
	std::cout << "Pop::euler finishing" << std::endl;
}

/* Reproduction:
Only canopy trees can reproduce, hence we stop at the last tree that can
reproduce (i.e., the last tree taller than s*).
In Purves2008, the reproduction is proportional to the total sun-exposed crown
area (which is 0 for understorey trees). That is to say:
	Reproduction = fecundity x crownArea(s, s*)
The fecundity parameter is from table S4 Purves2008. It is a difficult parameter
to estimate (check his Appendix2: Parameter estimation)
*/
double Population::reproduction()
{
	double popReprod = 0;
	c_cohort_it it = m_cohortsVec.cbegin();
	double G0 = m_species->v(0, m_s_star, m_env->annual_mean_temperature, m_env->annual_precipitation);
	double fecundity = m_species->fecundity;

	while (it->m_mu > m_s_star && it != m_cohortsVec.end())
	{
		popReprod += it->crownArea(m_s_star) * it->m_lambda;
		++it;
	}

	popReprod *= fecundity/G0;
	return popReprod;
}

/* Competition calculation:
To calculate the competition, I start from the biggest tree, and calculate the
its crown area. As long as I have not found s* such that the sum of the crown
areas of trees above s* equals the ground area, I decrement s* of delta_s. If
such an s* cannot be found, I set s* = 0. It means there is a gap in the canopy.

Remarks:
1. This calculation is costly, and has to be done at each time step!
2. In the case of the flat-top model, there is no need to integrate over s, we
can directly jump from cohorts to cohorts as the crown area is a step function
of size s in this particular case. The flat-top case is also discontinuous, and
there is no guarantee of having a solution for summedArea == 1.
*/
void Population::competition()
{
	double summedArea = 0;
	c_cohort_it it_cohort = m_cohortsVec.cbegin();
	m_s_star = it_cohort->m_mu;

	while ((m_s_star > 0) && (summedArea < 1)) // needs a plot_area later? Is it a density and thus 1?
	{
		// Reset the values
		summedArea = 0;
		it_cohort = m_cohortsVec.cbegin();
		while ((it_cohort != m_cohortsVec.cend()) && (it_cohort->m_mu >= m_s_star))
		{
			// crown area evaluated at s* given cohort of size m_mu, multiplied by density m_lambda of the cohort
			summedArea = it_cohort->crownArea(m_s_star)*it_cohort->m_lambda;
			++it_cohort;
		}
		m_s_star = m_s_star - m_delta_s;
	}
	if (m_s_star < 0)
		m_s_star = 0;
}

/* To compute the basal area of the community and the 'number of individuals'
The basal area is supposed to be in m²/ha (square meter per hectar). However, my
plot surface unit is square meter, and my diameter unit is millimiter. I there-
fore coerce the output to fit currently used units
The number of individuals is supposed to be (spatial case):
	\int_{x \in Plot} \int_{0}^{+Inf} N(s, x, t) ds dx
In my case, the density is per m per square meter (the first meter is for the
size s, the square meter is for the ground area). I therefore conclude:
	\int_{x \in Plot} \int_{0}^{+Inf} N(s, x, t) ds dx
		= \int_{x \in Plot} 1 dx \times \int_{0}^{+Inf} N(s, x, t) ds
		= plotArea \times \int_{0}^{+Inf} N(s, x, t) ds
*/
void Population::totalDensity_basalArea()
{
	double current_BA = 0;
	double currentDensity = 0;

	c_cohort_it it = m_cohortsVec.cbegin();
	c_cohort_it it_lim = m_cohortsVec.cbegin() + m_nonZeroCohort;
	for (; it != it_lim; ++it)
	{
		current_BA += it->m_mu*it->m_mu*it->m_lambda; // dbh² \times density per m per m². The += is the 'integral over s'
		currentDensity += it->m_lambda; //
	}

	double conversion = 10000.0/1000000.0; // Conversion factor to have it in m²/ha: 10 000 to get ha, 1 000 000 to get m² from mm²
	current_BA *= M_PI/4*conversion; // π dbh^2/4
	m_basalArea = current_BA;
	m_totalDensity = currentDensity;
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream& operator<<(std::ostream& os, Population const &pop)
{
	for (c_cohort_it it = pop.m_cohortsVec.begin(); it != pop.m_cohortsVec.end(); ++it)
		os << *it << std::endl;
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
	// std::cout << "mergeCohorts called" << std::endl;
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
	// std::cout << "Here 2 " << *it << std::endl;
	// std::cout << *this << std::endl;
	// std::cout << "-------------" << m_nonZeroCohort << std::endl;

	ref_it = it;
	// Second (if required), merge similar cohorts
	while (ref_it != lim_it)
	{
		mu = 0;
		lambda = 0;
		it = ref_it + 1;
		// abs val for float = fabs, useless in the next while loop because population is sorted
		while ((it != lim_it) && (ref_it->m_mu - it->m_mu < thresholdSimilarity)) // BUG
		{
			mu += it->m_lambda * it->m_mu;
			lambda += it->m_lambda;
			this->resetCohorts(it);
			mergedCohorts = true;
			++it;
		} // END BUG
		ref_it->m_mu = (ref_it->m_lambda * ref_it->m_mu + mu)/(ref_it->m_lambda + lambda);
		ref_it->m_lambda += lambda;
		ref_it = it;
	}
	// std::cout << "Here 3" << std::endl;
	// std::cout << *this << std::endl;
	// std::cout << "-------------" << m_nonZeroCohort << std::endl;

	// Reseting cohorts considered negligeable
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
	for (; it != lim_it; ++it)
		std::cout << *it << std::endl;
	std::cout << "Canopy height s*: " << this->m_s_star << std::endl;
	std::cout << "Number of zeros and boundary cohorts: " << m_cohortsVec.size() - m_nonZeroCohort << std::endl;
}
