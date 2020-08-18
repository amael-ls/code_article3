
#ifndef FOREST_H
#define FOREST_H

// Official headers
#include <vector>

// My headers
#include "Landscape.h++"
#include "Population.h++"

class Forest
{
	public :
		void spatialDynamics();

		// Return pointers to the neighbour cells
		std::vector<Population*> neighbours(unsigned int const target); //, double const deltaX, double const deltaY, double const maxDispersalDist

	private :
		Landscape* m_land;
		std::vector<Population> m_popVec; // de taille inf a land. Chq element -> 1 elem de land
};

// Notes
// push_back dans fec, + order by Rlang

#endif
