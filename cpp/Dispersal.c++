
#ifndef DISPERSAL_C
#define DISPERSAL_C

// My headers
#include "Dispersal.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Dispersal::Dispersal(Species const* const sp, std::string const climateFilename):
	m_species(sp), m_dispersalProbaThreshold(sp->dispersalProbaThreshold), m_min_dispersalProba(sp->min_dispersalProba),
	m_dispersalDistThreshold(sp->dispersalDistThreshold), m_max_dispersalDist(sp->max_dispersalDist)
{
	/**** Read landscape parameters from file climateFilename ****/
	par::Params climateParams(climateFilename.c_str(), "=");

	m_nRow_land = climateParams.get_val<unsigned int>("nRow");
	m_nCol_land = climateParams.get_val<unsigned int>("nCol");
	m_dim_land = m_nRow_land*m_nCol_land;

	m_deltaLon = climateParams.get_val<unsigned int>("deltaLon");
	m_deltaLat = climateParams.get_val<unsigned int>("deltaLat");
}
/*************************************************/
/******        Dispersal integration        ******/
/*************************************************/
void Dispersal::kernel(double x, double xminusa, double bminusx, double &y, void *ptr) const
{
	double *param = (double *) ptr; // casting
	double const z = param[0];

	y = m_species->K(x, z); //y = x*x*cos(sqrt(x - 3)) + 2 - cos(z)*sin(z) + exp(-z*z/2.0);
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
	mySelf->kernel(x, xminusa, bminusx, y, ptr);
}

void landscapeIntegrals(Dispersal& disp)
{
	double distanceToZero;
	double currentLongitude;
	double currentLatitude;

	double deltaLon(disp.m_deltaLon);
	double deltaLat(disp.m_deltaLat);

	pt2Object = (void*) &disp;

	double value2d = 0;

	/**** Compute the integral for all possible distances in landscape ****/
	// cf Remark 2 in the header file for comments (especially why distance to 0 only).
	for (unsigned int row = 0; row < disp.m_nRow_land; ++row) // latitude direction
	{
		for (unsigned int col = 0; col < disp.m_nCol_land; ++col) // longitude direction
		{
			currentLongitude = col*deltaLon;
			currentLatitude = row*deltaLat;
			distanceToZero = std::sqrt(currentLongitude*currentLongitude + currentLatitude*currentLatitude);
			
			double arrayParams[2] = {currentLongitude, currentLongitude + deltaLon};
			double (*params)[2] = &arrayParams;

			alglib::autogkstate ss;
			alglib::autogkreport reprep;
			alglib::autogksmooth(currentLatitude, currentLatitude + deltaLat, ss);
			alglib::autogkintegrate(ss, Dispersal::wrapper_To_Call_Kintegral, params);
			alglib::autogkresults(ss, value2d, reprep);

			// if (distanceToZero < 20000)
				// std::cout << distanceToZero << std::endl;

			if ((disp.m_dispersalProbaThreshold) &&(value2d > disp.m_min_dispersalProba)) // If using proba threshold
			{
				if (disp.m_map_distance_integral.find(distanceToZero) == disp.m_map_distance_integral.end())
					disp.m_map_distance_integral[distanceToZero] = value2d;
			}

			if ((disp.m_dispersalDistThreshold) && (distanceToZero < disp.m_dispersalDistThreshold)) // If using distance threshold
			{
				// std::cout << distanceToZero << "    " << disp.m_max_dispersalDist << std::endl;
				if (disp.m_map_distance_integral.find(distanceToZero) == disp.m_map_distance_integral.end())
					disp.m_map_distance_integral[distanceToZero] = value2d;
			}
		}
	}
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream& operator<<(std::ostream& os, Dispersal const &dispersal)
{
	// Species
	dispersal.m_species->printName(os);

	// Landscape
	os << std::endl << "Dimensions (row x col): " << dispersal.m_nRow_land << " x " << dispersal.m_nCol_land << std::endl;
	os << "Resolution (lon x lat): " << dispersal.m_deltaLon << " x " << dispersal.m_deltaLat << std::endl;

	// Sum integrals
	std::map<double, double>::const_iterator map_it = (dispersal.m_map_distance_integral).cbegin();
	double totIntegral = 0;
	for (; map_it != (dispersal.m_map_distance_integral).cend(); ++map_it)
		totIntegral += map_it->second;
	
	os << "Integral on Γ: " << totIntegral << std::endl;

	os << "Dimension of the map: " << (dispersal.m_map_distance_integral).size();
	return os;
}

#endif
