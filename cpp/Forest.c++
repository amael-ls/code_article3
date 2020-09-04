
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

		if (m_nRow_land < 1 || m_nCol_land < 1)
			throw Except_Forest(m_nRow_land, m_nCol_land);

		/**** Creating folders for saving dynamics ****/
		// Checking and creating folder m_pathCompReprodFile if necessary
		c_species_it species_it = m_speciesList.cbegin();
		std::string path_compReprod, path_popDyn;
		bool folderCreated;

		std::cout << "List of species:" << std::endl;
		for (; species_it != m_speciesList.cend(); ++species_it)
		{
			std::cout << "    - " << (*species_it)->m_speciesName << std::endl;

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

				// Create Patch which initialise the populations
				m_patchVec.emplace_back(Patch(env, m_speciesList, m_initPath, m_initFilenamePattern, m_maxCohorts));

				++counterPatch;
				if (counterPatch > m_dim_land)
					throw Except_Landscape(m_dim_land, climateFile);
			}
		}

		if (counterPatch != m_dim_land)
			throw Except_Landscape(m_dim_land, counterPatch);

		// Sort forest
		this->sort(true);

		std::vector<Cohort *> test;
		m_patchVec[34].getAllNonZeroCohorts(test);

		// Compute Δlongitude and Δlatitude. No need to compute max and min of lon and lat: landscape is sorted
		// We assume the lattice is regular. Otherwise Δlongitude and Δlatitude should both be in Environment.
		c_patch_it patch = m_patchVec.cbegin();
		if (m_nCol_land == 1 && m_nRow_land == 1)
		{
			m_deltaLon = sqrt(patch->m_env.plotArea);
			m_deltaLat = m_deltaLon;
		}
		else if (m_nCol_land == 1 && m_nRow_land > 1)
		{
			m_deltaLat = (patch->m_env).distance(std::next(patch)->m_env); // The next patch is also the next latitude
			m_deltaLon = patch->m_env.plotArea/m_deltaLat;
		}
		else if (m_nCol_land > 1) // Whatever m_nRow
		{
			m_deltaLon = (patch->m_env).distance(std::next(patch)->m_env);
			m_deltaLat = patch->m_env.plotArea/m_deltaLon;
		}

		// Compute basal area and density for each patch

		std::cout << "Forest constructed with success, using file <" << m_forestParamsFilename << ">" << std::endl;
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
}

void Forest::patchDynamics(double const t, double const delta_t)
{
	patch_it targetPatch;
	for (targetPatch = m_patchVec.begin(); targetPatch != m_patchVec.end(); ++targetPatch)
	{
		if (targetPatch->m_isPopulated) // Do the following computation only when necessary
			targetPatch->populationDynamics(t, delta_t); // Compute local seed bank, local competition, and age local population for each species
	}
}

void Forest::spatialDynamics()
{
	// Time variables
	double t;
	double delta_t = (m_tmax - m_t0)/(m_nIter - 1);

	// Iterators
	patch_it targetPatch_it, neighbour_it;

	// Variables for neighbours and recruitment
	std::vector<int> boundingBox;
	double targetSeedBank, seedContribution;

	// Others
	std::string compReprodFilename;
	std::string popDynFilename;

	// // Compute the total seed production for each population (the 2nd sum in eq 27 -> l index), and apply growth and mortality to each cohort of the pop (eq 14)
	// if (!m_saveOnlyLast)
	// {
	// 	for (unsigned int i = 1; i < m_nIter; ++i) // time loop
	// 	{
	// 		t = m_t0 + (i - 1)*delta_t; // i starts at 1, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
	// 		for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
	// 		{
	// 			pop_it->reproduction(); // Compute local seed bank for each pop
	// 			pop_it->euler(t, delta_t);
	// 			(pop_it->m_currentIter)++;
				
	// 			if (i % m_freqSave == 0)
	// 			{
	// 				pop_it->m_compReprod_ofs << t + delta_t << " " << pop_it->m_localProducedSeeds << " "
	// 					<< pop_it->m_s_star << " " << pop_it->m_basalArea << " " << pop_it->m_totalDensity << std::endl;

	// 				pop_it->m_popDyn_ofs << *pop_it;
	// 			}
	// 		}

	// 		for (targetEnv_it = m_land->m_envVec.begin(); targetEnv_it != m_land->m_envVec.end(); ++targetEnv_it) // Starting reproduction prog
	// 		{
	// 			neighbours_indices((*targetEnv_it)->m_patchId, boundingBox); // boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};
	// 			for (unsigned int row = boundingBox[0]; row < boundingBox[3]; ++row)
	// 			{
	// 				for (unsigned int col = boundingBox[1]; col < boundingBox[2]; ++col)
	// 				{
	// 					// index = row*m_nCol_land + col;
	// 					// std::cout << index << std::endl;
	// 				}
	// 			}
	// 		}
	// 	}

	// 	// Save final time
	// 	for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
	// 	{
	// 		pop_it->m_compReprod_ofs << t + delta_t << " " << pop_it->m_localProducedSeeds << " "
	// 			<< pop_it->m_s_star << " " << pop_it->m_basalArea << " " << pop_it->m_totalDensity << std::endl;

	// 		pop_it->m_popDyn_ofs << *pop_it;
	// 	}
	// }
	// else
	{
		for (unsigned int iter = 1; iter < m_nIter; ++iter) // time loop, starts at 1 because the initial condition is considered the 0th iteration
		{
			t = m_t0 + (iter - 1)*delta_t; // iter starts at 1, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
			this->patchDynamics(t, delta_t);

			for (targetPatch_it = m_patchVec.begin(); targetPatch_it != m_patchVec.end(); ++targetPatch_it) // Starting reproduction prog
			{
				// Reset local seed bank
				targetSeedBank = 0;

				// Get neighbours
				neighbours_indices((targetPatch_it->m_env).m_patchId, boundingBox); // boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};

				// Cover all the sources within neighbours to collect dispersed seeds

				for (unsigned int row = boundingBox[0]; row <= boundingBox[3]; ++row) // Important: less than or equal to (<=)
				{
					for (unsigned int col = boundingBox[1]; col <= boundingBox[2]; ++col) // Important: less than or equal to (<=)
					{
						sourceEnv = m_land->m_envVec[row*m_nCol_land + col]; // index = row*m_nCol_land + col;
						if ((sourceEnv->m_isPresent).find(m_sp) != (sourceEnv->m_isPresent).end()) // if species is present in source
						{
							std::cout << "Entered" << std::endl;
							// Compute contribution from source to target seed bank
							seedContribution = (sourceEnv->m_isPresent)[m_sp]->m_localProducedSeeds; // x integral(K)
							
							// Update local seed bank used for target
							targetSeedBank += seedContribution;

							// Update source seed bank
							(sourceEnv->m_isPresent)[m_sp]->m_localProducedSeeds -= seedContribution;
							std::cout << "Exited" << std::endl;
						}
					}
				}
				std::cout << "Target: " << (*targetEnv_it)->m_patchId << ", " << targetSeedBank << std::endl;
				if ((*targetEnv_it)->m_isPresent.find(m_sp) != ((*targetEnv_it)->m_isPresent).end()) // if species is present in target
				{
					if ((*targetEnv_it)->m_patchId == 34)
						targetSeedBank = 1;

					if (targetSeedBank > 0)
					{
						std::cout << "Allo: " << (((*targetEnv_it)->m_isPresent)[m_sp]->m_env)->m_patchId << std::endl;
						// Update the target local seed bank
						((*targetEnv_it)->m_isPresent)[m_sp]->m_localSeedBank = targetSeedBank;

						// Apply recruitment to the population, which apply euler to the boundary cohort and reset the population seedBank to 0
						((*targetEnv_it)->m_isPresent)[m_sp]->recruitment(t, delta_t);
					}
				}
				else
				{
					if (targetSeedBank > 0)
					{
						// Create new Population with one cohort
						// --- Outputs' filename
						std::cout << "No pop in target" << std::endl;
						compReprodFilename = m_pathCompReprodFile + m_compReprodFilePattern + std::to_string((*targetEnv_it)->m_patchId) + ".txt";
						popDynFilename = m_pathPopDynFile + m_popDynFilePattern + std::to_string((*targetEnv_it)->m_patchId) + ".txt";

						// --- Add population to the forest
						m_popVec.emplace_back(Population(m_maxCohorts, m_sp, targetSeedBank, *targetEnv_it, i, compReprodFilename, popDynFilename));
						std::cout << "emplace_back done" << std::endl;

						// Update the forest

						// Sort the population vector in the same order than landscape

						// The presence of species sp in the target environment is already set via population constructor
					}
				}
				
			}
		}

		// Save final time
		for (pop_it = m_popVec.begin(); pop_it != m_popVec.end(); ++pop_it)
		{
			pop_it->m_compReprod_ofs << t + delta_t << " " << pop_it->m_localProducedSeeds << " "
				<< pop_it->m_s_star << " " << pop_it->m_basalArea << " " << pop_it->m_totalDensity << std::endl;
			pop_it->m_popDyn_ofs << *pop_it;
		}
	}
	
	// Close output files when simulation is done
	this->close_ofs();
}

void Forest::neighbours_indices(unsigned int const target, std::vector<int>& boundingBox, Species const* species) const
{
	int topLeft_r, topLeft_c, topRight_r, topRight_c, bottomLeft_r, bottomLeft_c, bottomRight_r, bottomRight_c;
	double maxDispersalDist;

	if (species->max_dispersalDist)
		maxDispersalDist = species->dispersalDistThreshold;
	else if (species->min_dispersalProba) // It is actually stupid to compute it everytime... sp should get a function to compute maxDispersalDist and run it at construction
		maxDispersalDist = 100; // To compute, use a dichotomy... while (integral from d to + inf > minDispProba) {++d}
	else
		maxDispersalDist = 100; // default value

	unsigned int influenceRadius_x = std::ceil(maxDispersalDist/m_deltaLon);
	unsigned int influenceRadius_y = std::ceil(maxDispersalDist/m_deltaLat);

	unsigned int col_ind = target % m_nCol_land;
	unsigned int row_ind = (int) (target - col_ind)/m_nCol_land;

	// Default neighbour is itself
	topLeft_c = col_ind;
	topRight_c = col_ind;
	topLeft_r = row_ind;
	bottomLeft_r = row_ind;

	if (influenceRadius_x >= m_deltaLon) // If more than itself is covered in longitude direction
	{
		topLeft_c = std::max(col_ind - influenceRadius_x, 0U); // 0U for unsigned int
		topRight_c = std::min(col_ind + influenceRadius_x, m_nCol_land); // Casting required
	}

	if (influenceRadius_y >= m_deltaLat) // If more than itself is covered in latitude direction
	{
		topLeft_r = std::max(row_ind - influenceRadius_y, 0U);
		bottomLeft_r = std::min(row_ind + influenceRadius_y, m_nRow_land); // Casting required
	}
	boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};
}

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
