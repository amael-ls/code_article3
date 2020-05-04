
#ifndef COHORT_H
#define COHORT_H

#include <vector>

#include "Tree.h++"

class Cohort
{
	friend class Population;
	public :
		Cohort();
		Cohort(double const lambda, double const mu);

		double crownArea() const;
		double reproduction() const;

		std::vector<double> ODE_II(double const t, double const s_star);
		std::vector<double> ODE_V(double const t, double const s_star, double const popReprod);
		void euler(double const t, double const delta_t, double const s_star);
		void euler(double const t, double const delta_t, double const s_star,
			std::vector<double> (Cohort::*ode)(double, double));
		void euler(double const t, double const delta_t, double const s_star, double const popReprod,
			std::vector<double> (Cohort::*ode)(double, double, double));

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Cohort const &cohort);
		friend bool operator<(Cohort const& cohort1, Cohort const& cohort2);
		friend bool operator>(Cohort const& cohort1, Cohort const& cohort2);

	private :
		double m_lambda, m_mu;
		double m_crownArea;
		// Tree const *m_tree; // tree is constant, not the pointer so no problem with vector
};

// external functions
/*
bool operator<(Cohort const& cohort1, Cohort const& cohort2); // overload to sort cohorts by height (increase order by default)
std::ostream& operator<<(std::ostream& os, Cohort const &cohort);
*/
#endif
