
/*
	Remarks:
		1. Ref: Purves 2008, Predicting and understanding forest dynamics using a simple tractable model. PNAS
		2. In Purves 2008, the dbh is in mm
*/

#ifndef COHORT_C
#define COHORT_C

// Official headers
#include <iomanip> // std::setw, std::left
#include <cmath> // for log, M_PI

// My headers
#include "Species.h++"
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
	m_lambda(0), m_mu(0), m_species(NULL)
{}


Cohort::Cohort(Species const *sp) :
	m_lambda(0), m_mu(0), m_species(sp)
{}

Cohort::Cohort(Cohort const& cohort) :
	m_lambda(cohort.m_lambda), m_mu(cohort.m_mu), m_species(cohort.m_species)
{}

Cohort::Cohort(double const lambda, double const mu, Species const *sp) :
	m_lambda(lambda), m_mu(mu), m_species(sp)
{}

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
double Cohort::crownArea(double const s_star) const
{
	// Access parameters from *Species to calculate the crown radius
	// Table S2 Purves 2008
	double R0_C0 = m_species->R0_C0;
	double R0_C1 = m_species->R0_C1;

	double R40_C0 = m_species->R40_C0;
	double R40_C1 = m_species->R40_C1;

	double M_C0 = m_species->M_C0;
	double M_C1 = m_species->M_C1;

	double B_C0 = m_species->B_C0;
	double B_C1 = m_species->B_C1;

	double a = m_species->a;
	double b = m_species->b;
	double T_param = m_species->T_param;

	// Appendix S3, Eq S3.3 (erroneously denoted S2.3 in the article)
	double R0 = (1 - T_param)*R0_C0 + T_param*R0_C1;
	double R40 = (1 - T_param)*R40_C0 + T_param*R40_C1;

	double M = (1 - T_param)*M_C0 + T_param*M_C1;
	double B = (1 - T_param)*B_C0 + T_param*B_C1;

	// Calculate potential max radius Eq S1.6, watch out dbh in cm in Purves 2008!
	double Rp_max = R0 + (R40 - R0)*m_mu/400;

	// Convert dbh to height. If dbh is in mm, then height is in m. Trick: x^n = Exp[n Log[x]]
	double height = std::exp((a - b + b*std::log10(m_mu))*std::log(10));
	double height_star = std::exp((a - b + b*std::log10(s_star))*std::log(10));

	if (height == 0)
		return 0;

	// Calculate distance to the top
	double distToTop = height - height_star;

	// Vector crown radius
	double crownRadius = Rp_max * std::exp(B*std::log(std::min(distToTop, height*M) / (height*M))); // R_{i, y}^p

	return M_PI*crownRadius*crownRadius;
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
		mu = y[1], averaged size (the i-state)
		This function represents the dynamics of a cohort at the boundary condition (i.e., newborns)
	*/
	std::vector<double> y (2);

	double pi = m_lambda*m_mu;

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
}

// Euler method for ODE V
void Cohort::euler(double const t, double const delta_t, double const s_star, Environment const& env,
	double const popReprod, std::vector<double> (Cohort::*ode)(double, Environment const&, double))
{
	std::vector<double> fy = (this->*ode)(s_star, env, popReprod);
	m_lambda = m_lambda + delta_t * fy[0];
	m_mu = m_mu + delta_t * fy[1];
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream &operator<<(std::ostream &os, Cohort const& cohort)
{
	// os << std::setw(10);
	// os << std::setprecision(5);
	os << cohort.m_lambda << " " << cohort.m_mu;
	return os;
}

bool operator<(Cohort const& cohort1, Cohort const& cohort2)
{
    return (cohort1.m_mu < cohort2.m_mu);
}

bool operator>(Cohort const& cohort1, Cohort const& cohort2)
{
    return (cohort1.m_mu > cohort2.m_mu);
}

#endif
