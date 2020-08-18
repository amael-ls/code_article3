
#ifndef FOREST_C
#define FOREST_C

// Official headers
#include <filesystem> // To list files from folder, experimental/filesystem is now deprecated
#include <cmath>

// My headers
#include "Forest.h++"

Forest::Forest() {}

void Forest::spatialDynamics()
{
	std::vector<Population>::iterator target_it;
	std::vector<int> boundingBox;

	unsigned int counter = 0;
	for (target_it = m_popVec.begin(); target_it != m_popVec.end(); ++target_it)
	{
		;
		++counter;
	}
}

void Forest::neighbours_indices(unsigned int const target, std::vector<int>& boundingBox) const
{
	int topLeft_r, topLeft_c, topRight_r, topRight_c, bottomLeft_r, bottomLeft_c, bottomRight_r, bottomRight_c;

	int nRow(m_land->m_nRow), nCol(m_land->m_nCol), dim(m_land->m_dim);
	double const deltaX (m_land->m_deltaLon);
	double const deltaY(m_land->m_deltaLat);

	double maxDispersalDist;

	if (m_sp->max_dispersalDist)
		maxDispersalDist = m_sp->dispersalDistThreshold;
	else if (m_sp->min_dispersalProba)
		maxDispersalDist = 100; // To compute, use a dichotomy... while (integral from d to + inf > minDispProba) {++d}

	int influenceRadius_x = std::ceil(maxDispersalDist/deltaX);
	int influenceRadius_y = std::ceil(maxDispersalDist/deltaY);

	int col_ind = target % nCol;
	int row_ind = (int) (target - col_ind)/nCol;

	topLeft_r = std::max(row_ind - influenceRadius_y, 0);
	topLeft_c = std::max(col_ind - influenceRadius_x, 0);
	topRight_c = std::min(col_ind + influenceRadius_x, nCol);
	bottomLeft_r = std::min(row_ind + influenceRadius_y, nRow);

	boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};
}

#endif
