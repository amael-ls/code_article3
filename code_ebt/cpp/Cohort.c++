
/*
	Remarks:
		1. Ref: Purves 2008, Predicting and understanding forest dynamics using a simple tractable model. PNAS
		2. In Purves 2008, the dbh is in mm
		3. Ref: LeSquin2022, Development of a spatial Escalator Boxcar Train algorithm for sessile species: application to the
			boreal-temperate forest ecotone. ???
*/

#ifndef COHORT_C
#define COHORT_C

// My headers
#include "Cohort.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
/*
	The constructors are build in this order:
	First constructor (default): Everything is set to 0 or NULL pointer

	Second constructor: Everything is set to 0, except the species

	Third constructor: Density λ, diameter μ, and species set
*/
Cohort::Cohort() :
	m_lambda(0), m_mu(0), m_height(0), m_species(NULL)
{}


Cohort::Cohort(Species const *sp, unsigned int birthIteration) :
	m_lambda(0), m_mu(0), m_height(0), m_species(sp), m_birthIteration(birthIteration)
{}

Cohort::Cohort(Cohort const& cohort, unsigned int birthIteration) :
	m_lambda(cohort.m_lambda), m_mu(cohort.m_mu), m_species(cohort.m_species), m_birthIteration(birthIteration)
{
	// Convert dbh to height. If dbh is in mm, then height is in m. Trick: x^n = Exp[n Log[x]]
	m_height = std::exp((m_species->a - m_species->b + m_species->b*std::log10(m_mu))*std::log(10));
}

Cohort::Cohort(double const lambda, double const mu, Species const *sp, unsigned int birthIteration) :
	m_lambda(lambda), m_mu(mu), m_species(sp), m_birthIteration(birthIteration)
{
	// Convert dbh to height. If dbh is in mm, then height is in m. Trick: x^n = Exp[n Log[x]]
	m_height = std::exp((m_species->a - m_species->b + m_species->b*std::log10(m_mu))*std::log(10));
}

/*******************************************/
/******        Characteristics        ******/
/*******************************************/
/* The crown area is from Purves2008's formulae.
Let an individual be of height h. To compute the crown area of this individual
at a given height h*, it is necessary to calculate the distance h - h*. If this
is negative, then the area is 0 (because the tree is smaller than h*). Otherwise,
it is Purves2008's formulae defined in his appendices 2 and 3. Note that because
we work with diameter rather than height, it is necessary to first convert it to
a height equivalent (using allometries also in Purves2008).
All the required equations and values are described in the function, and where to
find them in Purves2008's appendices.
*/
double Cohort::crownArea(double const height_star) const
{
	return std::exp(-m_mu); // It is exp(-s) in the article, with s being the dbh!
}


/************************************/
/******        Dynamics        ******/
/************************************/
// Second order scheme, from de Roos (1988), ODE (II). This is described in the paper
std::vector<double> Cohort::ODE_II(double const s_star, Environment const& env)
{
	/* DESCRIPTION:
		lambda = y[0], number of individuals in the cohort
		mu = y[1], averaged size (the i-state)
		This function represents the dynamics of a cohort along its characteristics
	*/
	double temperature_growth = env.annual_mean_temperature;
	double precipitation_growth = env.annual_precipitation;
	double temperature_mortality = env.min_temperature_of_coldest_month;
	double precipitation_mortality = env.precipitation_of_driest_quarter;

	std::vector<double> y (2);
	y[0] = -m_species->d(m_mu, s_star, temperature_mortality, precipitation_mortality) * m_lambda;
	y[1] = m_species->v(m_mu, s_star, temperature_growth, precipitation_growth);

	return y;
}

// Second order scheme, from de Roos (1988), ODE (V). This is described in the paper
std::vector<double> Cohort::ODE_V(double const s_star, Environment const& env, double const popReprod)
{
	/* DESCRIPTION:
		lambda = y[0], number of individuals in the cohort
		pi = y[1], averaged size (the i-state)
		This function represents the dynamics of a cohort at the boundary condition (i.e., newborns)
	*/
	std::vector<double> y (2);

	double pi = m_lambda*m_mu; // See equation 15 in reference Le Squin et. al. 2022 (See README.md)

	double temperature_growth = env.annual_mean_temperature;
	double precipitation_growth = env.annual_precipitation;
	double temperature_mortality = env.min_temperature_of_coldest_month;
	double precipitation_mortality = env.precipitation_of_driest_quarter;

	y[0] = -m_species->d(0, s_star, temperature_mortality, precipitation_mortality) * m_lambda -
		m_species->dd_ds(0, s_star, temperature_mortality, precipitation_mortality) * pi + popReprod;

	y[1] = m_species->v(0, s_star, temperature_growth, precipitation_growth) *
		m_lambda + m_species->dv_ds(0, s_star, temperature_growth, precipitation_growth) * pi -
		m_species->d(0, s_star, temperature_mortality, precipitation_mortality) * pi;

	return y;
}

// Euler method for ODE II
void Cohort::euler(double const t, double const delta_t, double const s_star, Environment const& env,
	std::vector<double> (Cohort::*ode)(double, Environment const&))
{
	std::vector<double> fy = (this->*ode)(s_star, env);
	m_lambda = m_lambda + delta_t * fy[0];
	m_mu = m_mu + delta_t * fy[1];

	// Update height
	m_height = std::exp((m_species->a - m_species->b + m_species->b*std::log10(m_mu))*std::log(10));
}

// Euler method for ODE V
void Cohort::euler(double const t, double const delta_t, double const s_star, Environment const& env,
	double const popReprod, std::vector<double> (Cohort::*ode)(double, Environment const&, double))
{
	std::vector<double> fy = (this->*ode)(s_star, env, popReprod);
	m_lambda = m_lambda + delta_t * fy[0];
	m_mu = m_mu + delta_t * fy[1];

	// Update height
	m_height = std::exp((m_species->a - m_species->b + m_species->b*std::log10(m_mu))*std::log(10));
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream &operator<<(std::ostream &os, Cohort const& cohort)
{
	// os << std::setw(10);
	// os << std::setprecision(5); 
	os << cohort.m_birthIteration << " " << cohort.m_lambda << " " << cohort.m_mu << " " << cohort.m_height;
	// os << "B: " << cohort.m_birthIteration << std::endl << "λ " << cohort.m_lambda << std::endl << "μ " << cohort.m_mu << std::endl << "H " << cohort.m_height << "------";
	return os;
}

bool operator<(Cohort const& cohort1, Cohort const& cohort2) // Compare height, not dbh, which is required for the multi-species case
{
	return (cohort1.m_height < cohort2.m_height);
}

bool operator>(Cohort const& cohort1, Cohort const& cohort2) // Compare height, not dbh, which is required for the multi-species case
{
	return (cohort1.m_height > cohort2.m_height);
}

bool greaterCohortPtr(Cohort const* cohort1, Cohort const* cohort2) // Compare height, not dbh, which is required for the multi-species case
{
	return (cohort1->m_height > cohort2->m_height);
}

#endif
