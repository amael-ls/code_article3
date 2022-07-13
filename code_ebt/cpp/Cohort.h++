
/* Description of the class
The class Cohort represents an ensemble of trees of the same species, same size
at a given place. I organise it:
	- A density (number of individuals per unit size) and diameter
	- A species pointer

I list the functions here, but describe them in the associated c++ file:
	- Three constructors, default (for vector)
	- Characteristics (crown area, number of offspring produced)
	- ODEs and Euler (to solve the ODEs that represent my PDEs model)
	- Overloading for ease

Remarks:
1. A crown area associated to the cohort is only useful in the flat-top model.
Indeed, in this case, the crown area is the same whatever the distance to the
top of the tree (or 0 if above the top of the tree). Therefore, the crown area
does not need to be reevaluate at each size-step when searching for s*.
*/

#ifndef COHORT_H
#define COHORT_H

// Official headers
// Official headers
#include <functional> // std::multiplies
#include <algorithm> // std::transform
#include <iomanip> // std::setw, std::left
#include <vector>
#include <cmath> // for log, M_PI

// My headers
#include "Environment.h++"
#include "Species.h++"

class Cohort
{
	friend class Population;
	friend class Patch;

	public :
		// Constructors
		Cohort();
		Cohort(Species const *sp, unsigned int birthIteration);
		Cohort(Cohort const& cohort, unsigned int birthIteration);
		Cohort(double const lambda, double const mu, Species const *sp, unsigned int birthIteration);

		// Others
		double crownArea(double const height_star) const;
		double reproduction() const;

		// Dynamics
		std::vector<double> ODE_II(double const s_star, Environment const& env);
		std::vector<double> ODE_V(double const s_star, Environment const& env, double const popReprod);

		void euler(double const t, double const delta_t, double const s_star, Environment const& env,
			std::vector<double> (Cohort::*ode)(double, Environment const&));
		void euler(double const t, double const delta_t, double const s_star, Environment const& env,
			double const popReprod, std::vector<double> (Cohort::*ode)(double, Environment const&, double));

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Cohort const &cohort);
		friend bool operator<(Cohort const& cohort1, Cohort const& cohort2);
		friend bool operator>(Cohort const& cohort1, Cohort const& cohort2);
		friend bool greaterCohortPtr(Cohort const* cohort1, Cohort const* cohort2);

	private :
		double m_lambda, m_mu, m_height; // Density, dbh, height
		Species const *m_species; // Species is constant, not the pointer so no problem with vector
		unsigned int m_birthIteration; // Iteration at which the cohort appeared
};

// External function
bool greaterCohortPtr(Cohort const* cohort1, Cohort const* cohort2);

#endif
