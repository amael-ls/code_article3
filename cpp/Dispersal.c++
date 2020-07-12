
#ifndef SPECIES_C
#define SPECIES_C

// Official headers

// My headers
#include "Dispersal.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Dispersal::Dispersal(Species const* const sp, Landscape const* const land, double const longitude, double const latitude):
	m_species(sp), m_landscape(land), m_longitude(longitude), m_latitude(latitude)
{
	std::vector<Environment*>::const_iterator it_landscape = land->m_envVec.cbegin();

	// If the species has a maximum dispersal distance, compute all distances from source
	// and keep only the patches within maximal distance
	if (sp->max_dispersalDist)
	{
		std::vector<double> distanceSourceReceiver;
	}
	for (; it_landscape != land->m_envVec.cend(); ++it_landscape)
		;
}

#endif
