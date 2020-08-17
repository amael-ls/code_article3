
#ifndef SPECIES_C
#define SPECIES_C

// Official headers

// ALGLIB headers
#include "integration.h"

// My headers
#include "Dispersal.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Dispersal::Dispersal(Species const* const sp, Landscape const* const land, Environment const* const popEnv):
	m_species(sp), m_landscape(land), m_popEnv(popEnv)
{
	// If the species has a maximum dispersal distance, compute all distances from source
	// and keep only the patches within maximal distance
	double proportion = 0;
	if (sp->max_dispersalDist)
	{
		double distanceSourceReceiver = 0;
		unsigned int i = 0;
		for (; i < land->m_envVec.size(); ++i)
		{
			distanceSourceReceiver = popEnv->distance(*(land->m_envVec[i]));
			if (distanceSourceReceiver < sp->dispersalDistThreshold)
			{
				m_indices.push_back(i);
				proportion = 0;
				m_proportions.push_back(2);
			}
		}
	}
	if (sp->min_dispersalProba)
	{
		std::vector<double> distanceSourceReceiver;
	}
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

#endif
