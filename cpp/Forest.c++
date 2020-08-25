
#ifndef FOREST_C
#define FOREST_C

// Official headers
#include <filesystem> // To list files from folder, experimental/filesystem is now deprecated
#include <string>
#include <cmath>

// My headers
#include "Forest.h++"

Forest::Forest(Landscape* land, Species *sp, std::string const init_path, std::string const init_filenamePattern, unsigned int const maxCohort) :
	m_land(land), m_sp(sp)
{
	m_nRow_land = m_land->m_nRow;
	m_nCol_land = m_land->m_nCol;
	m_dim_land = m_land->m_dim;

	std::string pathFile(init_path + init_filenamePattern);
	std::string init_filename;

	bool initialisedPatch;

	for (unsigned int i = 0; i < m_land->m_dim; ++i)
	{
		initialisedPatch = m_land->m_initLoc[i];
		if (initialisedPatch)
		{
			init_filename = pathFile + std::to_string((m_land->m_envVec[i])->m_patchId) + ".txt";
			m_popVec.emplace_back(Population(maxCohort, m_sp, init_filename, m_land->m_envVec[i], 0U)); // 0U for unsigned int
		}
	}
}

void Forest::spatialDynamics()
{
	std::vector<Population>::iterator target_it;
	std::vector<int> boundingBox;

	for (target_it = m_popVec.begin(); target_it != m_popVec.end(); ++target_it)
	{
		neighbours_indices((target_it->m_env)->m_patchId, boundingBox);
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

// Need to order Forest according to Spatial order

void Forest::print() const
{
	std::vector<Environment*>::const_iterator it = m_land->m_envVec.cbegin();
	unsigned int counter(0), popCounter(0);
	for (; it != m_land->m_envVec.cend(); ++it)
	{
		std::cout << (*it)->m_patchId << "    " << m_land->m_initLoc[counter] << std::endl;
		if ((*it)->m_initPopulated)
		{
			std::cout << m_popVec[popCounter] << std::endl;
			++popCounter;
		}
		++counter;
	}
}

#endif
