
#ifndef FOREST_C
#define FOREST_C

// Official headers
#include <filesystem> // To list files from folder, experimental/filesystem is now deprecated

// My headers
#include "Forest.h++"

Forest::Forest() {}

void Forest::spatialDynamics()
{
	std::vector<Population>::iterator target_it;
	for (target_it = m_popVec.begin(); target_it != m_popVec.end(); ++target_it)
	{
		;
	}
}

std::vector<Population*> Forest::neighbours(unsigned int const target)
{
	int topLeft, topRight, bottomLeft;
	, double const deltaX, double const deltaY, double const maxDispersalDist
}

#endif
