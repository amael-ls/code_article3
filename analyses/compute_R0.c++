/*
	This script is to compute R0 from Le Squin et. al, 2020, Climate-induced variation in the demography of 14 tree species is not sufficient to explain their distribution in eastern North America

	To compile the program, use g++ -I /usr/local/boost_1_79_0 compute_R0.c++ -o test
	To execute: ./test

	I must admit that I did it quick and dirty, with global variables... All of them are for Acer saccharum
*/

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm> // std::min
#include <iostream>
#include <cmath> // M_PI

// Global variables, the values correspond to the climate condition used in Orford.
// --- Mortality
double scaling_dbh_mu_M = 187.137;
double scaling_dbh_sd_M = 106.148;
double beta_0_M = -4.57597;
double beta_0_M_understorey = -4.96793;
double beta_dbh_M = -0.35367;
double beta_dbh_sq_M = 0.108749;

// --- Growth
double scaling_dbh_mu_G = 189.029;
double scaling_dbh_sd_G = 105.28;
double scaling_G_mu = 0.527844;
double scaling_G_sd = 0.868286;
double beta_0_G = 0.022408;
double beta_0_G_understorey = -0.223798;
double beta_1 = 0.160547;
double beta_2 = -0.0456582;

// --- Allometry crown
double R0_C0 = 0.503;
double R0_C1 = 3.126;

double R40_C0 = 0.5;
double R40_C1 = 10;

double M_C0 = 0.95;
double M_C1 = 0.95;

double B_C0 = 0.196;
double B_C1 = 0.511;

double T_param = 0.56;

double height_a = 0.4742; // For Acer saccharum
double height_b = 0.5743; // For Acer saccharum

// --- Height threshold (in m) and dbh threshold (in mm) at equilibrium for Abies balsamea in Orford
double height_star_abies = 31.72; // value based on a run with tmax = 2000, and nIter = 9001. Transect number 6, patch_id in [1600, 2601]
double dbh_star_acer = 614.4632; // 10^(1/b * (b - a + log(h)/log(10))), with h = height_star_abies

// double height_star_abies = 20.03259; // Actually, this is for Acer, but I was too lazy to change any name in the program...
// double dbh_star_acer = 276.0255; // 10^(1/b * (b - a + log(h)/log(10))), with h = height_star_acer

// Namespace
using namespace boost::math::quadrature;

// Functions
// --- Individual growth rate, dbh is the diameter
double growth_understorey(double dbh)
{
	dbh = (dbh - scaling_dbh_mu_G)/scaling_dbh_sd_G;

	// Growth function
	return std::exp(scaling_G_mu + scaling_G_sd * (beta_0_G_understorey + beta_1*dbh + beta_2*dbh*dbh));
}

double growth_overstorey(double dbh)
{
	dbh = (dbh - scaling_dbh_mu_G)/scaling_dbh_sd_G;

	// Growth function
	return std::exp(scaling_G_mu + scaling_G_sd * (beta_0_G + beta_1*dbh + beta_2*dbh*dbh));
}

// --- Individual death rate, dbh is the diameter
double mortality_understorey(double dbh)
{
	dbh = (dbh - scaling_dbh_mu_M)/scaling_dbh_sd_M;

	return std::exp(beta_0_M_understorey + beta_dbh_M*dbh + beta_dbh_sq_M*dbh*dbh);
}

double mortality_overstorey(double dbh)
{
	dbh = (dbh - scaling_dbh_mu_M)/scaling_dbh_sd_M;

	return std::exp(beta_0_M + beta_dbh_M*dbh + beta_dbh_sq_M*dbh*dbh);
}

double crownArea(double dbh, double height_star) // dbh in mm, height_star in m
{
	if (dbh == 0)
		return 0;
	
	// Appendix S3, Eq S3.3 (erroneously denoted S2.3 in the article)
	double R0 = (1 - T_param)*R0_C0 + T_param*R0_C1;
	double R40 = (1 - T_param)*R40_C0 + T_param*R40_C1;

	double M = (1 - T_param)*M_C0 + T_param*M_C1;
	double B = (1 - T_param)*B_C0 + T_param*B_C1;

	// Calculate potential max radius Eq S1.6, watch out dbh in cm in Purves 2008!
	double Rp_max = R0 + (R40 - R0)*dbh/400;

	// Calculate distance to the top
	double height = std::exp((height_a - height_b + height_b*std::log10(dbh))*std::log(10));
	double distToTop = height - height_star;

	// Vector crown radius
	double crownRadius = Rp_max * std::exp(B*std::log(std::min(distToTop, height*M) / (height*M))); // R_{i, y}^p

	return M_PI*crownRadius*crownRadius;
}

double recruitment_ratio_G(double dbh)
{
	return 0.0071*crownArea(dbh, height_star_abies)/growth_overstorey(dbh);
}

// --- Proportions survives from inf_bound to x
double ratio_understorey(double dbh)
{
	return mortality_understorey(dbh)/growth_understorey(dbh);
}

double ratio_overstorey(double dbh)
{
	return mortality_overstorey(dbh)/growth_overstorey(dbh);
}

double prop(double x, double inf_bound, bool isUnderstorey)
{
	double error;
	if (isUnderstorey)
		return std::exp(-gauss_kronrod<double, 15>::integrate(ratio_understorey, inf_bound, x, 5, 1e-9, &error));
	
	return std::exp(-gauss_kronrod<double, 15>::integrate(ratio_overstorey, inf_bound, x, 5, 1e-9, &error));
}

double integrand(double dbh)
{
	return recruitment_ratio_G(dbh) * prop(dbh, dbh_star_acer, false);
}

int main()
{
	double error = 0;
	double prop_under = prop(dbh_star_acer, 0, true);
	double integ = gauss_kronrod<double, 15>::integrate(integrand, dbh_star_acer, 1500, 5, 1e-9, &error);
	double total = prop_under * integ;
	std::cout << prop_under << std::endl;
	std::cout << integ << std::endl;
	std::cout << "R0 = " << total << std::endl;
	return 0;
}
