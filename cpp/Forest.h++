
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
		Forest(Landscape* land, Species* sp, std::string const init_path, std::string const init_filenamePattern, unsigned int const maxCohort);
		void spatialDynamics();

		// Return pointers to the neighbour cells
		void neighbours_indices(unsigned int const target, std::vector<int>& boundingBox) const;

		void print() const;

	private :
		Landscape* m_land;
		std::vector<Population> m_popVec; // size(m_popVec) <= size(m_land)
		Species* m_sp;

		// Landscape dimensions
		unsigned int m_nRow_land, m_nCol_land, m_dim_land;
};

// Notes
// push_back dans fec, + order by Rlang

#endif
