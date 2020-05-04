
#ifndef COHORT_C
#define COHORT_C

#include <iomanip> // std::setw, std::left
#include <cmath> // for log

#include "Demography.h++"
#include "Cohort.h++"


/****************************************/
/******        Constructors        ******/
/****************************************/
Cohort::Cohort() : //, Tree const *tree) :
	m_lambda(0), m_mu(0), m_crownArea(0) //, m_tree(tree)
{
//	double alpha(m_tree->getAlpha()), beta(m_tree->getBeta());
//	m_height = alpha*exp(beta*log(mu)); // Allometry from PNAS Predicting and understanding forest dynamics using a simple tractable model
}


Cohort::Cohort(double const lambda, double const mu) : //, Tree const *tree) :
	m_lambda(lambda), m_mu(mu) //, m_tree(tree)
{
	m_crownArea = this->crownArea();
	std::cout << "constructor: " << m_crownArea << std::endl;
//	double alpha(m_tree->getAlpha()), beta(m_tree->getBeta());
//	m_height = alpha*exp(beta*log(mu)); // Allometry from PNAS Predicting and understanding forest dynamics using a simple tractable model
}

/*******************************************/
/******        Characteristics        ******/
/*******************************************/
double Cohort::crownArea() const
{
	// Access parameters from *Tree, then calculate the area
	return std::exp(-m_mu);
}


/************************************/
/******        Dynamics        ******/
/************************************/
// Reproduction
double Cohort::reproduction() const // Is it useful?? Never called function 2019-09-18
{
	return 2;
}
// Second order scheme, from de Roos (1988), ODE (II)
std::vector<double> Cohort::ODE_II(double const t, double const s_star)
{
	/* DESCRIPTION:
		lambda = y[0], number of individuals in the cohort
		mu = y[1], averaged size (the i-state)
		This function represents the dynamics of a cohort along its characteristics
	*/
	std::vector<double> y (2);

	y[0] = -demography::d(t, m_mu, s_star) * m_lambda;
	y[1] = demography::v(t, m_mu, s_star);

	return y;
}


// Second order scheme, from de Roos (1988), ODE (V)
std::vector<double> Cohort::ODE_V(double const t, double const s_star, double const popReprod)
{
	/* DESCRIPTION:
		lambda = y[0], number of individuals in the cohort
		mu = y[1], averaged size (the i-state)
		This function represents the dynamics of a cohort along its characteristics
	*/
	std::vector<double> y (2);

	double pi = m_lambda*m_mu;

/*	std::cout << "d: " << demography::d(t, 0, s_star) << std::endl;
	std::cout << "v: " << demography::v(t, 0, s_star) << std::endl;
	std::cout << "l: " << m_lambda << std::endl;
	std::cout << "m: " << m_mu << std::endl;
*/
	y[0] = -demography::d(t, 0, s_star) * m_lambda - demography::dd_ds(t, 0, s_star) * pi +
		popReprod;
	y[1] = demography::v(t, m_mu, s_star) * m_lambda + demography::dv_ds(t, 0, s_star) * pi -
		demography::d(t, 0, s_star) * pi;

	return y;
}


void Cohort::euler(double const t, double const delta_t, double const s_star)
{
	std::vector<double> fy = this->ODE_II(t, s_star);
	m_lambda = m_lambda + delta_t * fy[0];
	m_mu = m_mu + delta_t * fy[1];
}

void Cohort::euler(double const t, double const delta_t, double const s_star,
	std::vector<double> (Cohort::*ode)(double, double))
{
	std::vector<double> fy = (this->*ode)(t, s_star);
	m_lambda = m_lambda + delta_t * fy[0];
	m_mu = m_mu + delta_t * fy[1];
}

void Cohort::euler(double const t, double const delta_t, double const s_star, double const popReprod,
	std::vector<double> (Cohort::*ode)(double, double, double))
{
	std::vector<double> fy = (this->*ode)(t, s_star, popReprod);
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
	os << cohort.m_lambda << "," << cohort.m_mu;
	return os;
}


bool operator<(Cohort const& cohort1, Cohort const& cohort2)
{
	// std::cout << cohort1.m_mu << "\t" << cohort2.m_mu << std::endl;
	// std::cout << std::boolalpha;
	// bool aa = cohort1.m_mu < cohort2.m_mu;
	// bool bb = cohort2.m_mu < cohort1.m_mu;
	//
	// std::cout << aa << "\t" << bb << std::endl;

    return (cohort1.m_mu < cohort2.m_mu);
}

bool operator>(Cohort const& cohort1, Cohort const& cohort2)
{
    return (cohort1.m_mu > cohort2.m_mu);
}

#endif
