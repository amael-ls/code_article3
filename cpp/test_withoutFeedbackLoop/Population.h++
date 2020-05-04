
#ifndef POPULATION_H
#define POPULATION_H

// Official headers
#include <iostream>
#include <vector>

// My headers
#include "Cohort.h++"

class Population
{
	public :
		Population(unsigned int const maxCohorts, double const s_inf);
		Population(unsigned int const maxCohorts, double const s_inf,
			std::vector<double> const & lambda, std::vector<double> const & mu);
		Population(unsigned int const maxCohorts, double const s_inf,
			std::vector<Cohort> const & cohorts);

		void euler(unsigned int n_t, double t0, double t_max);
		void euler2(unsigned int n_t, double t0, double t_max);
		double reproduction();
		void competition();

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Population const &pop);

		// Sorting and organising
		void sort(bool const decreasingOrder);
		void lastReproducer();
		bool mergeCohorts(double const thresholdSimilarity, double const thresholdDensity);
		void resetCohorts(std::vector<Cohort>::iterator const it);
		void printNonZero() const;

	private :
		unsigned int m_maxCohorts; // maximal number of cohorts
		unsigned int m_nonZeroCohort; // Number of non empty cohorts
		std::vector<Cohort>::iterator m_lastReproducer; // Last cohort able to reproduce (given sorted by decreasing order)
		double const m_s_inf, m_delta_s;
		std::vector<Cohort> m_cohortsVec; // The population of cohorts
		double m_s_star; // Competition
		// Tree const *m_tree; // tree is constant, not the pointer so no problem with vector
};

#endif
