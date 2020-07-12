
/* Description of the class
The class Dispersal represents the set of patches accessible from a population at a
given place. I organise it:
	- A vector of unsigned integers storing the indices of patches reachable by the population
	- A vector of doubles computing the proportion seeds dispersed 

I list the functions here, but describe them in the associated c++ file:
	- A constructor
*/

#ifndef DISPERSAL_H
#define DISPERSAL_H

// Official headers
#include <vector>

// My headers
#include "landscape.h++"
#include "Species.h++"

class Dispersal
{
	friend class Population;

	public :
		// Constructors
		Dispersal(Species const* const sp, Landscape const* const land, double const longitude, double const latitude);

	private :
	Species const* const m_species;
	Landscape const* const m_landscape;
	double const m_longitude;
	double const m_latitude;
	std::vector<unsigned int> m_indices;
	std::vector<double> m_proportions;
};

#endif
