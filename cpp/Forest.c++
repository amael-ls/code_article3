
#ifndef FOREST_C
#define FOREST_C

// Official headers
#include <filesystem> // To list files from folder, experimental/filesystem is now deprecated
#include <fstream>
#include <string>
#include <cmath>

// My headers
#include "Error_classes.h++"
#include "Forest.h++"
#include "Params.h++"

Forest::Forest(Landscape* land, Species *sp, std::string const forestParamsFilename) :
	m_forestParamsFilename(forestParamsFilename), m_land(land), m_sp(sp)
{
	try
	{
		/**** Read forest parameters from file forestParamsFilename ****/
		// Load params
		par::Params forestParams(forestParamsFilename.c_str(), " = ");

		// Initial condition
		m_initFilenamePattern = forestParams.get_val<std::string>("initFilenamePattern");
		m_initPath = forestParams.get_val<std::string>("initPath");

		// Saving options
		m_compReprodFilePattern = forestParams.get_val<std::string>("compReprodFilePattern");
		m_pathCompReprodFile = forestParams.get_val<std::string>("pathCompReprodFile");
		m_popTimeFile = forestParams.get_val<std::string>("popTimeFile");
		m_freqSave = forestParams.get_val<unsigned int>("freqSave");

		// Simulation parameters
		m_t0 = forestParams.get_val<double>("t0");
		m_tmax = forestParams.get_val<double>("tmax");
		m_nIter = forestParams.get_val<unsigned int>("nIter");
		m_maxCohorts = forestParams.get_val<unsigned int>("maxCohorts");

		// get rasterOrder_Rlang parameter
		std::string rasterOrder_Rlang = forestParams.get_val<std::string>("rasterOrder_Rlang");
		std::transform(rasterOrder_Rlang.begin(), rasterOrder_Rlang.end(), rasterOrder_Rlang.begin(),
			[](unsigned char c){ return std::tolower(c); });

		if (rasterOrder_Rlang == "true")
			m_rasterOrder_Rlang = true;

		// Get saveOnlyLast parameter
		std::string saveOnlyLast = forestParams.get_val<std::string>("saveOnlyLast");
		std::transform(saveOnlyLast.begin(), saveOnlyLast.end(), saveOnlyLast.begin(),
			[](unsigned char c){ return std::tolower(c); });

		if (saveOnlyLast == "true")
		{
			m_saveOnlyLast = true;
			std::cout << "Only the last iteration will be saved, despite freqSave = " << m_freqSave << std::endl;
		}

		if (m_freqSave > m_tmax && !m_saveOnlyLast)
			throw Except_Forest(m_freqSave, m_nIter);

		// Dimensions from landscape
		m_nRow_land = m_land->m_nRow;
		m_nCol_land = m_land->m_nCol;
		m_dim_land = m_land->m_dim;

		// Checking and creating folder m_pathCompReprodFile if necessary
		if (! std::filesystem::exists(m_pathCompReprodFile))
		{
			std::filesystem::create_directories(m_pathCompReprodFile);
			std::cout << "Directory <" << m_pathCompReprodFile << "> successfully created" << std::endl;
		}

		// Initialise cohorts existing at time 0
		std::string pathFile(m_initPath + m_initFilenamePattern);
		std::string init_filename;
		std::string compReprodFilename;

		bool initialisedPatch;

		for (unsigned int i = 0; i < m_land->m_dim; ++i)
		{
			initialisedPatch = m_land->m_initLoc[i];
			if (initialisedPatch)
			{
				init_filename = pathFile + std::to_string((m_land->m_envVec[i])->m_patchId) + ".txt";
				compReprodFilename = m_pathCompReprodFile + m_compReprodFilePattern + std::to_string((m_land->m_envVec[i])->m_patchId) + ".txt";
				m_popVec.emplace_back(Population(m_maxCohorts, m_sp, init_filename, m_land->m_envVec[i], 0U, compReprodFilename)); // 0U = unsigned int
			}
		}
		std::cout << "Forest constructed with success, using file <" << m_forestParamsFilename << ">" << std::endl;
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
}

void Forest::spatialDynamics()
{
	// Time variables
	double t;
	double delta_t = (m_tmax - m_t0)/(m_nIter - 1);

	// Iterators
	std::vector<Population>::iterator pop_it;
	std::vector<Environment*>::iterator env_it;

	std::vector<int> boundingBox;

	// Compute the total seed production for each population (the 2nd sum in eq 27 -> l index), and apply growth and mortality to each cohort of the pop (eq 14)
	for (unsigned int i = 1; i < m_nIter; ++i) // time loop
	{
		t = m_t0 + (i - 1)*delta_t; // i starts at 1, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
		for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
		{
			pop_it->reproduction(); // Compute local seed bank for each pop
			// Call euler
		}

		for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
		{
			neighbours_indices((pop_it->m_env)->m_patchId, boundingBox);
		}
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
	else if (m_sp->min_dispersalProba) // It is actually stupid to compute it everytime... sp should get a function to compute maxDispersalDist and run it at construction
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
void Forest::sort(bool const rasterOrder_Rlang)
{

}

// Print function
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

// Create output files (one per patch) containing local seed bank, s*, basal area, and total density
/*
void Forest::createOutputFiles() const
{
	std::vector<Environment*>::const_iterator env_it;
	std::string name;
	for (env_it = (m_land->m_envVec).cbegin(); env_it != (m_land->m_envVec).cend(); ++env_it)
	{
		std::ofstream outputCompReprod;
		name = m_pathCompReprodFile + m_compReprodFilePattern + std::to_string((*env_it)->m_patchId);
		
		outputCompReprod.open(name, std::ofstream::trunc); // Discard the content of the file if already exists

		if(!outputCompReprod.is_open())
		{
			std::stringstream ss;
			ss << "*** ERROR (from Population::euler): cannot open output files";
			throw (std::runtime_error (ss.str()));
		}

		outputCompReprod << "time reproduction competition basalArea totalDensity" << std::endl;
	}
}
*/

#endif
