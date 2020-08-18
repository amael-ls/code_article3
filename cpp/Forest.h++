
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
		void neighbours_indices(unsigned int const target, std::vector<int>& boundingBox) const;

	private :
		Landscape* m_land;
		std::vector<Population> m_popVec; // de taille inf a land. Chq element -> 1 elem de land
		Species *m_sp;
};

// Notes
// push_back dans fec, + order by Rlang

#endif
