
/*
	Remarks:
		1. Ref: Purves 2008, Predicting and understanding forest dynamics using a simple tractable model. PNAS
		2. In Purves 2008, the dbh is in mm
*/

#ifndef COHORT_C
#define COHORT_C

// Official headers
#include <functional> // std::multiplies
#include <algorithm> // std::transform
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
	// // Access parameters from *Species to calculate the crown radius
	// // Table S2 Purves 2008
	// double R0_C0 = m_species->R0_C0;
	// double R0_C1 = m_species->R0_C1;

	// double R40_C0 = m_species->R40_C0;
	// double R40_C1 = m_species->R40_C1;

	// double M_C0 = m_species->M_C0;
	// double M_C1 = m_species->M_C1;

	// double B_C0 = m_species->B_C0;
	// double B_C1 = m_species->B_C1;

	// double a = m_species->a;
	// double b = m_species->b;
	// double T_param = m_species->T_param;

	// // Appendix S3, Eq S3.3 (erroneously denoted S2.3 in the article)
	// double R0 = (1 - T_param)*R0_C0 + T_param*R0_C1;
	// double R40 = (1 - T_param)*R40_C0 + T_param*R40_C1;

	// double M = (1 - T_param)*M_C0 + T_param*M_C1;
	// double B = (1 - T_param)*B_C0 + T_param*B_C1;

	// // Calculate potential max radius Eq S1.6, watch out dbh in cm in Purves 2008!
	// double Rp_max = R0 + (R40 - R0)*m_mu/400;

	// // Convert dbh to height. If dbh is in mm, then height is in m. Trick: x^n = Exp[n Log[x]]
	// double height = std::exp((a - b + b*std::log10(m_mu))*std::log(10));
	// double height_star = std::exp((a - b + b*std::log10(s_star))*std::log(10));

	// if (height == 0)
	// 	return 0;

	// // Calculate distance to the top
	// double distToTop = height - height_star;

	// // Vector crown radius
	// double crownRadius = Rp_max * std::exp(B*std::log(std::min(distToTop, height*M) / (height*M))); // R_{i, y}^p

	// return M_PI*crownRadius*crownRadius;

	return (std::exp(-m_mu*m_mu/2.0 - m_mu));
}


/************************************/
/******        Dynamics        ******/
/************************************/
// Second order scheme, from de Roos (1988), ODE (II). This is described in the paper
std::vector<double> Cohort::ODE_II(double const s_star, Environment const& env)
{
	/* DESCRIPTION:
		m_lambda = y[0], number of individuals in the cohort
		m_mu = y[1], averaged size (the i-state)
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

// Second order scheme, from de Roos (1988), ODE (II) for RK4 methods
std::vector<double> Cohort::ODE_II(double const lambda, double const mu, double const s_star, Environment const& env)
{
	/* DESCRIPTION:
		m_lambda = y[0], number of individuals in the cohort
		m_mu = y[1], averaged size (the i-state)
		This function represents the dynamics of a cohort along its characteristics
	*/
	double temperature_growth = env.annual_mean_temperature;
	double precipitation_growth = env.annual_precipitation;
	double temperature_mortality = env.min_temperature_of_coldest_month;
	double precipitation_mortality = env.precipitation_of_driest_quarter;

	std::vector<double> y (2);
	y[0] = -m_species->d(m_mu + mu, s_star, temperature_mortality, precipitation_mortality) * (m_lambda + lambda);
	y[1] = m_species->v(m_mu + mu, s_star, temperature_growth, precipitation_growth);

	return y;
}

// Second order scheme, from de Roos (1988), ODE (V). This is described in the paper
std::vector<double> Cohort::ODE_V(double const s_star, Environment const& env, double const popReprod)
{
	/* DESCRIPTION:
		m_lambda = y[0], number of individuals in the cohort
		m_mu = y[1], averaged size (the i-state)
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

// Second order scheme, from de Roos (1988), ODE (V) for RK4
std::vector<double> Cohort::ODE_V(double const lambda, double const mu, double const s_star, Environment const& env, double const popReprod)
{
	/* DESCRIPTION:
		m_lambda = y[0], number of individuals in the cohort
		m_mu = y[1], averaged size (the i-state)
		This function represents the dynamics of a cohort at the boundary condition (i.e., newborns)
	*/
	std::vector<double> y (2);

	double pi = (m_lambda + lambda)*(m_mu + mu);

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

// Runge-Kutta 4 method for ODE II. Note that this function is only for autonomous systems
void Cohort::rk4(double const t, double const delta_t, double const s_star, Environment const& env,
	std::vector<double> (Cohort::*ode)(double, double, double, Environment const&))
{
	std::vector<double> fy = (this->*ode)(0, 0, s_star, env);
	std::vector<double> k1(2), k2(2), k3(2), k4(2);
	
	// k1 = Δt f(t_i, y_i)
	std::transform(fy.begin(), fy.end(), k1.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// k2 = Δt f(t_i + 1/2 Δt, y_i + 1/2 k1)
	k2 = (this->*ode)(1.0/2.0 * k1[0], 1.0/2.0 * k1[1], s_star, env);
	std::transform(k2.begin(), k2.end(), k2.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// k3 = Δt f(t_i + 1/2 Δt, y_i + 1/2 k2)
	k3 = (this->*ode)(1.0/2.0 * k2[0], 1.0/2.0 * k2[1], s_star, env);
	std::transform(k3.begin(), k3.end(), k3.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// k4 = Δt f(t_i + Δt, y_i + k3)
	k4 = (this->*ode)(k3[0], k3[1], s_star, env);
	std::transform(k4.begin(), k4.end(), k4.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// y_n+1 = y_n + 1/6 k1 + 1/3 k2 + 1/3 k3 + 1/6 k4
	for (unsigned int i = 0; i < k1.size(); ++i)
		fy[i] = 1.0/6.0 * k1[i] + 1.0/3.0 * k2[i] + 1.0/3.0 * k3[i] + 1.0/6.0 * k4[i];

	m_lambda = m_lambda + fy[0]; // k1, ..., k4 are already multiplied by delta_t
	m_mu = m_mu + fy[1]; // k1, ..., k4 are already multiplied by delta_t
}

// Runge-Kutta 4 method for ODE V. Note that this function is only for autonomous systems
void Cohort::rk4(double const t, double const delta_t, double const s_star, Environment const& env,
	double const popReprod, std::vector<double> (Cohort::*ode)(double, double, double, Environment const&, double))
{
	std::vector<double> fy = (this->*ode)(0, 0, s_star, env, popReprod);

		std::vector<double> k1(2), k2(2), k3(2), k4(2);
	
	// k1 = Δt f(t_i, y_i)
	std::transform(fy.begin(), fy.end(), k1.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// k2 = Δt f(t_i + 1/2 Δt, y_i + 1/2 k1)
	k2 = (this->*ode)(1.0/2.0 * k1[0], 1.0/2.0 * k1[1], s_star, env, popReprod);
	std::transform(k2.begin(), k2.end(), k2.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// k3 = Δt f(t_i + 1/2 Δt, y_i + 1/2 k2)
	k3 = (this->*ode)(1.0/2.0 * k2[0], 1.0/2.0 * k2[1], s_star, env, popReprod);
	std::transform(k3.begin(), k3.end(), k3.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// k4 = Δt f(t_i + Δt, y_i + k3)
	k4 = (this->*ode)(k3[0], k3[1], s_star, env, popReprod);
	std::transform(k4.begin(), k4.end(), k4.begin(),
		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

	// y_n+1 = y_n + 1/6 k1 + 1/3 k2 + 1/3 k3 + 1/6 k4
	for (unsigned int i = 0; i < k1.size(); ++i)
		fy[i] = 1.0/6.0 * k1[i] + 1.0/3.0 * k2[i] + 1.0/3.0 * k3[i] + 1.0/6.0 * k4[i];

	m_lambda = m_lambda + fy[0]; // k1, ..., k4 are already multiplied by delta_t
	m_mu = m_mu + fy[1]; // k1, ..., k4 are already multiplied by delta_t
}

// // Runge-Kutta-Dormand-Prince RK-45 method for ODE II. Note that this function is only for autonomous systems
// double Cohort::rk45(double const t, double const delta_t, double const s_star, Environment const& env,
// 	std::vector<double> (Cohort::*ode)(double, double, double, Environment const&))
// {
// 	std::vector<double> fy = (this->*ode)(0, 0, s_star, env);
// 	std::vector<double> k1, k2, k3, k4, k5, k6, k7;

// 	// k1 = Δt f(t_i, y_i)
// 	std::transform(fy.begin(), fy.end(), k1.begin(),
// 		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

// 	// k2 = Δt f(t_i + 1/5 Δt, y_i + 1/5 k1)
// 	k2 = (this->*ode)(1.0/4.0 * k1[0], 1.0/5.0 * k1[1], s_star, env);
// 	std::transform(k2.begin(), k2.end(), k2.begin(),
// 		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

// 	// k3 = Δt f(t_i + 3/10 Δt, y_i + 3/40 k1 + 9/40 k2)
// 	k3 = (this->*ode)(3.0/40.0 * k1[0] + 9.0/40.0 * k2[0],
// 		3.0/40.0 * k1[1] + 9.0/40.0 * k2[1], s_star, env);
// 	std::transform(k3.begin(), k3.end(), k3.begin(),
// 		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

// 	// k4 = Δt f(t_i + 4/5 Δt, y_i + 44/45 k1 - 56/15 k2 + 32/9 k3)
// 	k4 = (this->*ode)(44.0/45.0 * k1[0] - 56.0/15.0 * k2[0] + 32.0/9.0 * k3[0],
// 		44.0/45.0 * k1[1] - 56.0/15.0 * k2[1] + 32.0/9.0 * k3[1], s_star, env);
// 	std::transform(k4.begin(), k4.end(), k4.begin(),
// 		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

// 	// k5 = Δt f(t_i + 8/9 Δt, y_i + 19372/6561 k1 - 25360/2187 k2 + 64448/6561 k3 - 212/729 k4)
// 	k5 = (this->*ode)(19372.0/6561.0 * k1[0] - 25360.0/2187.0 * k2[0] + 64448.0/6561.0 * k3[0] - 212.0/729.0 * k4[0],
// 		19372.0/6561.0 * k1[1] - 25360.0/2187.0 * k2[1] + 64448.0/6561.0 * k3[1] - 212.0/729.0 * k4[1], s_star, env);
// 	std::transform(k5.begin(), k5.end(), k5.begin(),
// 		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

// 	// k6 = Δt f(t_i + Δt, y_i + 9017/3168 k1 - 355/33 k2 + 46732/5247 k3 + 49/176 k4 - 5103/18656 k5)
// 	k6 = (this->*ode)(9017.0/3168.0 * k1[0] - 355.0/33.0 * k2[0] + 46732.0/5247.0 * k3[0] + 49.0/176.0 * k4[0] - 5103.0/18656.0 * k5[0],
// 		9017.0/3168.0 * k1[1] - 355.0/33.0 * k2[1] + 46732.0/5247.0 * k3[1] + 49.0/176.0 * k4[1] - 5103.0/18656.0 * k5[1], s_star, env);
// 	std::transform(k6.begin(), k6.end(), k6.begin(),
// 		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

// 	// k7 = Δt f(t_i + Δt, y_i + 35/384 k1 + 0 k2 + 500/1113 k3 + 125/192 k4 - 2187/6784 k5 + 11/84 k6)
// 	k7 = (this->*ode)(35.0/384.0 * k1[0] + 500.0/1113.0 * k3[0] + 125.0/192.0 * k4[0] - 2187.0/6784.0 * k5[0] + 11.0/84.0 * k6[0],
// 		35.0/384.0 * k1[1] + 500.0/1113.0 * k3[1] + 125.0/192.0 * k4[1] - 2187.0/6784.0 * k5[1] + 11.0/84.0 * k6[1], s_star, env);
// 	std::transform(k7.begin(), k7.end(), k7.begin(),
// 		std::bind(std::multiplies<double>(), std::placeholders::_1, delta_t));

// 	std::vector<double> y_next(2), z_next (2);

// 	// y_next = 35/384 k1 + 0 k2 + 500/1113 k3 + 125/192 k4 - 2187/6784 k5 + 11/84 k6 + 0 k7
// 	for (unsigned int i = 0; i < k1.size(); ++i)
// 		y_next[i] = fy[i] + 35.0/384.0 * k1[i] + 500.0/1113.0 * k3[i] + 125.0/192.0 * k4[i] - 2187.0/6784.0 * k5[i] + 11.0/84.0 * k6[i];

// 	// z_next = 5179/57600 k1 + 0 k2 + 7571/16695 k3 + 393/640 k4 - 92097/339200 k5 + 187/2100 k6 + 1/40 k7
// 	for (unsigned int i = 0; i < k1.size(); ++i)
// 		z_next[i] = fy[i] + 5179.0/57600.0 * k1[i] + 7571.0/16695.0 * k3[i] + 393.0/640.0 * k4[i] - 92097.0/339200.0 * k5[i] + 187.0/2100.0 * k6[i] + 1.0/40.0 * k7[i];

// 	double err_vector; // RESTART HERE!!! this method. For now there might be a problem: I do not know how many steps I will use!
// 	// Before doing any adaptive method, I will program RK4

// 	m_lambda = m_lambda + delta_t * fy[0];
// 	m_mu = m_mu + delta_t * fy[1];
// }

// // Runge-Kutta-Dormand-Prince RK-45 method for ODE V. Note that this function is only for autonomous systems
// void Cohort::rk45(double const t, double const delta_t, double const s_star, Environment const& env,
// 	double const popReprod, std::vector<double> (Cohort::*ode)(double, Environment const&, double))
// {
// 	std::vector<double> fy = (this->*ode)(s_star, env, popReprod);
// 	m_lambda = m_lambda + delta_t * fy[0];
// 	m_mu = m_mu + delta_t * fy[1];
// }

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
