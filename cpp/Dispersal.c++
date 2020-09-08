
#ifndef DISPERSAL_C
#define DISPERSAL_C

// My headers
#include "Dispersal.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Dispersal::Dispersal(Species const* const sp, std::string const climateFilename):
	m_species(sp)
{
	/**** Read landscape parameters from file climateFilename ****/
	par::Params climateParams(climateFilename.c_str(), "=");

	m_nRow_land = climateParams.get_val<unsigned int>("nRow");
	m_nCol_land = climateParams.get_val<unsigned int>("nCol");
	m_dim_land = m_nRow_land*m_nCol_land;

	m_deltaLon = climateParams.get_val<unsigned int>("deltaLon");
	m_deltaLat = climateParams.get_val<unsigned int>("deltaLat");

	/**** Compute the integral for all possible distances in landscape ****/
	
}
/*************************************************/
/******        Dispersal integration        ******/
/*************************************************/
void Dispersal::Kintegrand_lon(double x, double xminusa, double bminusx, double &y, void *ptr) const
{
	double *param = (double *) ptr; // casting
	double const z = param[0];

	y = x*x*cos(sqrt(x - 3)) + 2 - cos(z)*sin(z) + exp(-z*z/2.0); // m_species->K(x, z); // 
}

void Dispersal::Kintegral_lon(double z, double xminusa, double bminusx, double &y, void *ptr) const // Watch out between integrand and integral!
{
	// y(z) = integral(integrand(x), from c to d)
	double *param = (double *) ptr; // casting
	// Integral bounds
	double c = param[0];
	double d = param[1];

	alglib::autogkstate s;
	alglib::autogkreport rep;
	alglib::autogksmooth(c, d, s);
	alglib::autogkintegrate(s, Dispersal::wrapper_To_Call_Kintegral_lon, &z);
	alglib::autogkresults(s, y, rep);
}

void Dispersal::wrapper_To_Call_Kintegral(double x, double xminusa, double bminusx, double &y, void *ptr)
{
	// explicitly cast global variable <pt2Object> to a pointer to Dispersal
	// warning: <pt2Object> MUST point to an appropriate object!
	Dispersal* mySelf = (Dispersal*) pt2Object;

	// call member
	mySelf->Kintegral_lon(x, xminusa, bminusx, y, ptr);
}

void Dispersal::wrapper_To_Call_Kintegral_lon(double x, double xminusa, double bminusx, double &y, void *ptr)
{
	// explicitly cast global variable <pt2Object> to a pointer to Dispersal
	// warning: <pt2Object> MUST point to an appropriate object!
	Dispersal* mySelf = (Dispersal*) pt2Object;

	// call member
	mySelf->Kintegrand_lon(x, xminusa, bminusx, y, ptr);
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream& operator<<(std::ostream& os, Dispersal const &dispersal)
{
	// Species
	dispersal.m_species->printName(os);

	// Landscape
	os << "Dimensions (row x col): " << dispersal.m_nRow_land << " x " << dispersal.m_nCol_land << std::endl;
	os << "Resolution (lon x lat): " << dispersal.m_deltaLon << " x " << dispersal.m_deltaLat << std::endl;

	// Sum integrals
	std::map<double, double>::const_iterator map_it = (dispersal.m_distance_integral).cbegin();
	double totIntegral = 0;
	for (; map_it != (dispersal.m_distance_integral).cend(); ++map_it)
		totIntegral += map_it->second;
	
	os << "Integral on Γ: " << totIntegral << std::endl;

	os << "Dimension of the map: " << (dispersal.m_distance_integral).size();
	return os;
}

#endif
