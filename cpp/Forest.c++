
#ifndef FOREST_C
#define FOREST_C

// Official headers
#include <filesystem> // To list files from folder, experimental/filesystem is now deprecated
#include <algorithm> // std::sort
#include <fstream>
#include <string>
#include <cmath>

// My headers
#include "Error_classes.h++"
#include "Forest.h++"
#include "Params.h++"

// Define typedef shortcuts
typedef std::vector<Species*>::const_iterator c_species_it;
typedef std::vector<Patch>::const_iterator c_patch_it;
typedef std::vector<Patch>::iterator patch_it;

Forest::Forest(std::string const forestParamsFilename, std::vector<Species*> const speciesList, std::string const climateFilename) :
	m_forestParamsFilename(forestParamsFilename), m_speciesList(speciesList)
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
		m_popDynFilePattern = forestParams.get_val<std::string>("popDynFilePattern");
		m_pathPopDynFile = forestParams.get_val<std::string>("pathPopDynFile");
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
			throw Except_Forest(m_freqSave, m_nIter, m_dim_land, false);

		/**** Read landscape parameters from file climateFilename ****/
		par::Params climateParams(climateFilename.c_str(), "=");
		m_nRow_land = climateParams.get_val<unsigned int>("nRow");
		m_nCol_land = climateParams.get_val<unsigned int>("nCol");
		m_dim_land = m_nRow_land*m_nCol_land;
		std::string pathLandscape = climateParams.get_val<std::string>("path");
		std::string delimiter = climateParams.get_val<std::string>("delimiter");
		std::string climateFilenamePattern = climateParams.get_val<std::string>("filenamePattern");
		std::string isPopulated;

		/**** Creating folders for saving dynamics ****/
		// Checking and creating folder m_pathCompReprodFile if necessary
		c_species_it species_it = m_speciesList.cbegin();
		std::string path_compReprod, path_popDyn;
		bool folderCreated;

		std::cout << "List of species:" << std::endl;
		std::vector<std::string> speciesNames;

		for (; species_it != m_speciesList.cend(); ++species_it)
		{
			std::cout << "    - " << (*species_it)->m_speciesName << std::endl;
			speciesNames.emplace_back((*species_it)->m_speciesName);

			path_compReprod = m_pathCompReprodFile + (*species_it)->m_speciesName + "/";
			folderCreated = std::filesystem::create_directories(path_compReprod);
			if (folderCreated)
				std::cout << "Directory <" << path_compReprod << "> successfully created" << std::endl;

			// Checking and creating folder m_pathPopDynFile if necessary
			path_popDyn = m_pathPopDynFile + (*species_it)->m_speciesName + "/";
			folderCreated = std::filesystem::create_directories(path_popDyn);
			if (folderCreated)
				std::cout << "Directory <" << path_popDyn << "> successfully created" << std::endl;
		}

		/**** Creating forest of patches ****/
		// Fill the environment
		std::string climateFile;
		std::string initFile;
		unsigned int counterPatch = 0;
		bool foundAnInitFile = false;

		for(auto& p: std::filesystem::directory_iterator(pathLandscape))
		{
			climateFile = p.path().filename();
			if (climateFile.find(climateFilenamePattern) != std::string::npos)
			{
				climateFile = pathLandscape + "/" + climateFile;

				// Create environment
				Environment env(climateFile, delimiter);
				
				if (env.m_initPopulated)
				{
					for (species_it = m_speciesList.cbegin(); species_it != m_speciesList.cend(); ++species_it)
					{
						initFile = m_initPath + (*species_it)->m_speciesName + "/" + m_initFilenamePattern +
							std::to_string(env.m_patchId) + ".txt";
						
						if (std::filesystem::exists(initFile))
						{
							m_patchVec.emplace_back(Patch(env, *species_it, initFile,  m_maxCohorts));
							foundAnInitFile = true;
						}
						else
							m_patchVec.emplace_back(Patch(env, *species_it, m_maxCohorts));
					}
					if (!foundAnInitFile)
						throw Except_Forest(env.m_patchId, speciesNames);
					foundAnInitFile = false;
				}
				else
				{
					for (species_it = m_speciesList.cbegin(); species_it != m_speciesList.cend(); ++species_it)
							m_patchVec.emplace_back(Patch(env, *species_it, m_maxCohorts));
				}
				
				++counterPatch;
				if (counterPatch > m_dim_land)
					throw Except_Landscape(m_dim_land, climateFile);
			}
		}

		if (counterPatch != m_dim_land)
			throw Except_Landscape(m_dim_land, counterPatch);

		// Sort forest
		this->sort(m_rasterOrder_Rlang);

		// Compute competition for each patch

		// Compute basal area and density for each species and within each patch

		// Compute basal area and density for each patch

		std::cout << "Forest constructed with success, using file <" << m_forestParamsFilename << ">" << std::endl;
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
}

// void Forest::spatialDynamics()
// {
// 	// Time variables
// 	double t;
// 	double delta_t = (m_tmax - m_t0)/(m_nIter - 1);

// 	// Iterators
// 	std::vector<Population>::iterator pop_it;
// 	std::vector<Environment*>::iterator targetEnv_it;
// 	Environment* sourceEnv;

// 	// Variables for neighbours and recruitment
// 	std::vector<int> boundingBox;
// 	double targetSeedBank, seedContribution;

// 	// Others
// 	std::string compReprodFilename;
// 	std::string popDynFilename;

// 	// Compute the total seed production for each population (the 2nd sum in eq 27 -> l index), and apply growth and mortality to each cohort of the pop (eq 14)
// 	if (!m_saveOnlyLast)
// 	{
// 		for (unsigned int i = 1; i < m_nIter; ++i) // time loop
// 		{
// 			t = m_t0 + (i - 1)*delta_t; // i starts at 1, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
// 			for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
// 			{
// 				pop_it->reproduction(); // Compute local seed bank for each pop
// 				pop_it->euler(t, delta_t);
// 				(pop_it->m_currentIter)++;
				
// 				if (i % m_freqSave == 0)
// 				{
// 					pop_it->m_compReprod_ofs << t + delta_t << " " << pop_it->m_localProducedSeeds << " "
// 						<< pop_it->m_s_star << " " << pop_it->m_basalArea << " " << pop_it->m_totalDensity << std::endl;

// 					pop_it->m_popDyn_ofs << *pop_it;
// 				}
// 			}

// 			for (targetEnv_it = m_land->m_envVec.begin(); targetEnv_it != m_land->m_envVec.end(); ++targetEnv_it) // Starting reproduction prog
// 			{
// 				neighbours_indices((*targetEnv_it)->m_patchId, boundingBox); // boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};
// 				for (unsigned int row = boundingBox[0]; row < boundingBox[3]; ++row)
// 				{
// 					for (unsigned int col = boundingBox[1]; col < boundingBox[2]; ++col)
// 					{
// 						// index = row*m_nCol_land + col;
// 						// std::cout << index << std::endl;
// 					}
// 				}
// 			}
// 		}

// 		// Save final time
// 		for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
// 		{
// 			pop_it->m_compReprod_ofs << t + delta_t << " " << pop_it->m_localProducedSeeds << " "
// 				<< pop_it->m_s_star << " " << pop_it->m_basalArea << " " << pop_it->m_totalDensity << std::endl;

// 			pop_it->m_popDyn_ofs << *pop_it;
// 		}
// 	}
// 	else
// 	{
// 		for (unsigned int i = 1; i < 2; ++i) // m_nIter; ++i) // time loop
// 		{
// 			t = m_t0 + (i - 1)*delta_t; // i starts at 1, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
// 			for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
// 			{
// 				pop_it->reproduction(); // Compute local seed bank for each pop
// 				pop_it->euler(t, delta_t);	
// 				(pop_it->m_currentIter)++;
// 			}

// 			for (targetEnv_it = m_land->m_envVec.begin(); targetEnv_it != m_land->m_envVec.end(); ++targetEnv_it) // Starting reproduction prog
// 			{
// 				std::cout << "Bibibop: " << (*targetEnv_it)->m_patchId << std::endl;
// 				if ((*targetEnv_it)->m_patchId == 34)
// 					std::cout << "Bibibop2: " << (*targetEnv_it)->m_isPresent[m_sp]->m_env->m_patchId << std::endl;

// 				// Reset local seed bank
// 				targetSeedBank = 0;

// 				// Get neighbours
// 				neighbours_indices((*targetEnv_it)->m_patchId, boundingBox); // boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};

// 				// Cover all the sources within neighbours to collect dispersed seeds
// 				for (unsigned int row = boundingBox[0]; row <= boundingBox[3]; ++row) // Important: less than or equal to (<=)
// 				{
// 					for (unsigned int col = boundingBox[1]; col <= boundingBox[2]; ++col) // Important: less than or equal to (<=)
// 					{
// 						sourceEnv = m_land->m_envVec[row*m_nCol_land + col]; // index = row*m_nCol_land + col;
// 						if ((sourceEnv->m_isPresent).find(m_sp) != (sourceEnv->m_isPresent).end()) // if species is present in source
// 						{
// 							std::cout << "Entered" << std::endl;
// 							// Compute contribution from source to target seed bank
// 							seedContribution = (sourceEnv->m_isPresent)[m_sp]->m_localProducedSeeds; // x integral(K)
							
// 							// Update local seed bank used for target
// 							targetSeedBank += seedContribution;

// 							// Update source seed bank
// 							(sourceEnv->m_isPresent)[m_sp]->m_localProducedSeeds -= seedContribution;
// 							std::cout << "Exited" << std::endl;
// 						}
// 					}
// 				}
// 				std::cout << "Target: " << (*targetEnv_it)->m_patchId << ", " << targetSeedBank << std::endl;
// 				if ((*targetEnv_it)->m_isPresent.find(m_sp) != ((*targetEnv_it)->m_isPresent).end()) // if species is present in target
// 				{
// 					if ((*targetEnv_it)->m_patchId == 34)
// 						targetSeedBank = 1;

// 					if (targetSeedBank > 0)
// 					{
// 						std::cout << "Allo: " << (((*targetEnv_it)->m_isPresent)[m_sp]->m_env)->m_patchId << std::endl;
// 						// Update the target local seed bank
// 						((*targetEnv_it)->m_isPresent)[m_sp]->m_localSeedBank = targetSeedBank;

// 						// Apply recruitment to the population, which apply euler to the boundary cohort and reset the population seedBank to 0
// 						((*targetEnv_it)->m_isPresent)[m_sp]->recruitment(t, delta_t);
// 					}
// 				}
// 				else
// 				{
// 					if (targetSeedBank > 0)
// 					{
// 						// Create new Population with one cohort
// 						// --- Outputs' filename
// 						std::cout << "No pop in target" << std::endl;
// 						compReprodFilename = m_pathCompReprodFile + m_compReprodFilePattern + std::to_string((*targetEnv_it)->m_patchId) + ".txt";
// 						popDynFilename = m_pathPopDynFile + m_popDynFilePattern + std::to_string((*targetEnv_it)->m_patchId) + ".txt";

// 						// --- Add population to the forest
// 						m_popVec.emplace_back(Population(m_maxCohorts, m_sp, targetSeedBank, *targetEnv_it, i, compReprodFilename, popDynFilename));
// 						std::cout << "emplace_back done" << std::endl;

// 						// Update the forest

// 						// Sort the population vector in the same order than landscape

// 						// The presence of species sp in the target environment is already set via population constructor
// 					}
// 				}
				
// 			}
// 		}

// 		// Save final time
// 		for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
// 		{
// 			pop_it->m_compReprod_ofs << t + delta_t << " " << pop_it->m_localProducedSeeds << " "
// 				<< pop_it->m_s_star << " " << pop_it->m_basalArea << " " << pop_it->m_totalDensity << std::endl;
// 			pop_it->m_popDyn_ofs << *pop_it;
// 		}
// 	}
	
// 	// Close output files when simulation is done
// 	this->close_ofs();
// }

// void Forest::neighbours_indices(unsigned int const target, std::vector<int>& boundingBox) const
// {
// 	int topLeft_r, topLeft_c, topRight_r, topRight_c, bottomLeft_r, bottomLeft_c, bottomRight_r, bottomRight_c;

// 	int nRow(m_land->m_nRow), nCol(m_land->m_nCol), dim(m_land->m_dim);
// 	double const deltaX (m_land->m_deltaLon);
// 	double const deltaY(m_land->m_deltaLat);

// 	double maxDispersalDist;

// 	if (m_sp->max_dispersalDist)
// 		maxDispersalDist = m_sp->dispersalDistThreshold;
// 	else if (m_sp->min_dispersalProba) // It is actually stupid to compute it everytime... sp should get a function to compute maxDispersalDist and run it at construction
// 		maxDispersalDist = 100; // To compute, use a dichotomy... while (integral from d to + inf > minDispProba) {++d}
// 	else
// 		maxDispersalDist = 100; // default value

// 	int influenceRadius_x = std::ceil(maxDispersalDist/deltaX);
// 	int influenceRadius_y = std::ceil(maxDispersalDist/deltaY);

// 	int col_ind = target % nCol;
// 	int row_ind = (int) (target - col_ind)/nCol;

// 	// Default neighbour is itself
// 	topLeft_c = col_ind;
// 	topRight_c = col_ind;
// 	topLeft_r = row_ind;
// 	bottomLeft_r = row_ind;

// 	if (influenceRadius_x >= deltaX) // If more than itself is covered in longitude direction
// 	{
// 		topLeft_c = std::max(col_ind - influenceRadius_x, 0);
// 		topRight_c = std::min(col_ind + influenceRadius_x, nCol);
// 	}

// 	if (influenceRadius_y >= deltaY) // If more than itself is covered in latitude direction
// 	{
// 		topLeft_r = std::max(row_ind - influenceRadius_y, 0);
// 		bottomLeft_r = std::min(row_ind + influenceRadius_y, nRow);
// 	}
// 	boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};
// }

// // Need to order Forest according to Spatial order
// void Forest::sort(bool const rasterOrder_Rlang)
// {

// }

// // Print function
// void Forest::print() const
// {
// 	std::vector<Environment*>::const_iterator it = m_land->m_envVec.cbegin();
// 	unsigned int counter(0), popCounter(0);
// 	for (; it != m_land->m_envVec.cend(); ++it)
// 	{
// 		std::cout << (*it)->m_patchId << "    " << m_land->m_initLoc[counter] << std::endl;
// 		if ((*it)->m_initPopulated)
// 		{
// 			std::cout << m_popVec[popCounter] << std::endl;
// 			++popCounter;
// 		}
// 		++counter;
// 	}
// }

// // Close ofstreams
// void Forest::close_ofs()
// {
// 	std::vector<Population>::iterator pop_it;
// 	for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
// 		pop_it->m_compReprod_ofs.close();
// }

/***********************************/
/******        Sorting        ******/
/***********************************/
void Forest::sort(bool const rasterOrder_Rlang)
{
	// Sorting
	std::vector<Patch>::iterator first = m_patchVec.begin();
	std::vector<Patch>::iterator last = m_patchVec.end();

	if (rasterOrder_Rlang)
		std::sort(first, last);
	else
		std::sort(first, last, std::greater<Patch>());
}

// /************************************/
// /******        Overload        ******/
// /************************************/
std::ostream& operator<<(std::ostream& os, Forest const &forest)
{
	c_patch_it it = forest.m_patchVec.cbegin();
	for (; it != forest.m_patchVec.cend(); ++it)
		os << *it << std::endl;
	return os;
}

#endif
