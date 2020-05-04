
// Official headers
#include <algorithm> // std::sort
#include <exception> // std::throw
#include <iomanip> // std::setw, std::left
#include <string>
#include <cmath>

// My headers
#include "Population.h++"

// Define typedef shortcuts
typedef std::vector<Cohort>::iterator cohort_it;
typedef std::vector<Cohort>::const_iterator c_cohort_it;

/****************************************/
/******        Constructors        ******/
/****************************************/
Population::Population(unsigned int const maxCohorts, double const s_inf):
	m_maxCohorts(maxCohorts), m_nonZeroCohort(0), m_s_inf(s_inf),
	m_delta_s(s_inf/maxCohorts), m_cohortsVec(maxCohorts), m_s_star(0)
{
	// m_lastReproducer = m_cohortsVec.begin();
	m_nonZeroCohort = maxCohorts - 1;
	std::cout << m_maxCohorts << std::endl;
	std::cout << m_nonZeroCohort << std::endl;
	try
	{
		for (unsigned int i = 0; i < m_nonZeroCohort; ++i)
		{
			m_cohortsVec[i].m_lambda = std::exp(-m_delta_s*i); // /!\ i is unsigned int, so do not try exp(-i*delta_s)
			m_cohortsVec[i].m_mu = m_delta_s*i;
			m_cohortsVec[i].m_crownArea = m_cohortsVec[i].crownArea();
		}
	}
	catch(const std::out_of_range& ex)
	{
		std::stringstream ss;
		ss << "Warning: maximum population size exceded" << ex.what();
		throw (std::runtime_error (ss.str()));
	}
	this->sort(true);
	this->competition();
	this->lastReproducer();
}

Population::Population(unsigned int const maxCohorts, double const s_inf, std::vector<double> const & lambda,
	std::vector<double> const & mu): m_maxCohorts(maxCohorts), m_nonZeroCohort(lambda.size()),
	m_s_inf(s_inf), m_delta_s(s_inf/maxCohorts), m_cohortsVec(maxCohorts)
{
	try
	{
		for (unsigned int i = 0; i < m_nonZeroCohort; ++i)
		{
			m_cohortsVec[i].m_lambda = lambda[i];
			m_cohortsVec[i].m_mu = mu[i];
			m_cohortsVec[i].m_crownArea = m_cohortsVec[i].crownArea();
		}
	}
	catch(const std::out_of_range& ex)
	{
		std::stringstream ss;
		ss << "Warning: maximum population size exceded" << ex.what();
		throw (std::runtime_error (ss.str()));
	}
	this->sort(true);
	this->competition();
	this->lastReproducer();
}

Population::Population(unsigned int const maxCohorts, double const s_inf, std::vector<Cohort> const & cohorts):
	m_maxCohorts(maxCohorts), m_nonZeroCohort(cohorts.size()), m_s_inf(s_inf),
	m_delta_s(s_inf/maxCohorts), m_cohortsVec(cohorts)
{
	if (m_maxCohorts < m_nonZeroCohort)
	{
		m_maxCohorts = m_nonZeroCohort + 1;
		m_cohortsVec.resize(m_maxCohorts);
	}
	this->sort(true);
	this->competition();
	this->lastReproducer();
}

/********************************************/
/******        Euler & dynamics        ******/
/********************************************/
void Population::euler(unsigned int n_t, double t0, double t_max)
{
	double t;
	double delta_t = (t_max - t0)/(n_t - 1);

    for (unsigned int i = 0; i < n_t; ++i) // time loop
	{
		t = t0 + i*delta_t;
		for (cohort_it it = m_cohortsVec.begin(); it != m_cohortsVec.end(); ++it)
			it->euler(t, delta_t, m_s_star);
    }
}

void Population::euler2(unsigned int n_t, double t0, double t_max)
{
	double t;
	double delta_t = (t_max - t0)/(n_t - 1);
	double popReprod = 0;
	cohort_it it;
	cohort_it lim_it; // limit iterator

    for (unsigned int i = 0; i < n_t; ++i) // time loop
	{
		// Check there is one space to create a new cohort and merge/delete if required
		if (m_nonZeroCohort == m_maxCohorts)
			this->mergeCohorts(m_delta_s/8, 0.001);

		lim_it = m_cohortsVec.begin() + m_nonZeroCohort; // Check if this does not go beyond m_cohortsVec.end()

		// Integration within the size space Omega
		for (it = m_cohortsVec.begin(); it != lim_it; ++it)
			it->euler(t, delta_t, m_s_star, &Cohort::ODE_II);

		// Dynamics at the boundary condition in size
		popReprod = this->reproduction();
		it->euler(t, delta_t, m_s_star, popReprod, &Cohort::ODE_V);
		if (it->m_mu > m_delta_s) // If it reach the threshold, it becomes a cohort within Omega
			m_nonZeroCohort += 1;

		this->competition();
		t = t0 + i*delta_t;
    }
}

double Population::reproduction()
{
	double popReprod = 0;
	c_cohort_it it = m_cohortsVec.cbegin();
	for(; it != m_lastReproducer; ++it)
		popReprod += std::exp(-it->m_mu); // For analytical test 0.0071 * it->m_crownArea;
	// std::exp(-it->m_mu); // For analytical test
	return popReprod;
}

void Population::competition()
{
	double summedArea = 0;
	c_cohort_it it = m_cohortsVec.cbegin();
	while ((it != m_cohortsVec.cend()) && (summedArea < 1)) // needs a plot_area later? Is it a density and thus 1?
	{
		summedArea += it->crownArea();
		++it;
	}

	if (it == m_cohortsVec.cend())
		m_s_star = 0;
	else
		m_s_star = it->m_mu;
	// m_s_star = 0; // For analytical test
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream& operator<<(std::ostream& os, Population const &pop)
{
	std::cout << "Population size: " << pop.m_cohortsVec.size() << std::endl;
	os << "lambda,mu" << std::endl;
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

void Population::lastReproducer()
{
	// std::cout << "lastReproducer called" <<std::endl;
	m_lastReproducer = m_cohortsVec.begin();
	while (m_lastReproducer->m_mu > m_s_star && m_lastReproducer != m_cohortsVec.end())
		++m_lastReproducer;
}

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
	std::cout << "Number of zeros and boundary cohorts: " << m_cohortsVec.size() - m_nonZeroCohort << std::endl;
}
