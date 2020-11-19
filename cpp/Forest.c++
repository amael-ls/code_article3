
#ifndef FOREST_C
#define FOREST_C

// My headers
#include "Forest.h++"

// Define typedef shortcuts
typedef std::vector<Species*>::const_iterator c_species_it;
typedef std::vector<Patch>::const_iterator c_patch_it;
typedef std::vector<Patch>::iterator patch_it;

Forest::Forest(par::Params const& forestParameters, std::vector<Species*> const speciesList, std::string const climateFilename) :
	m_speciesList(speciesList)
{
	/**** Read forest parameters ****/
	// Initial condition
	m_initFilenamePattern = forestParameters.get_val<std::string>("initFilenamePattern");
	m_initPath = forestParameters.get_val<std::string>("initPath");
	checkPath(m_initPath, "initPath");

	// Saving options
	m_summaryFilePattern = forestParameters.get_val<std::string>("summaryFilePattern");
	m_summaryFilePath = forestParameters.get_val<std::string>("summaryFilePath");
	checkPath(m_summaryFilePath, "pathSummaryFile");

	m_popDynFilePattern = forestParameters.get_val<std::string>("popDynFilePattern");
	m_popDynFilePath = forestParameters.get_val<std::string>("popDynFilePath");
	checkPath(m_popDynFilePath, "popDynFilePath");

	m_freqSave = forestParameters.get_val<unsigned int>("freqSave");

	// Simulation parameters
	m_t0 = forestParameters.get_val<double>("t0");
	m_tmax = forestParameters.get_val<double>("tmax");
	m_nIter = forestParameters.get_val<unsigned int>("nIter");
	m_maxCohorts = forestParameters.get_val<unsigned int>("maxCohorts");

	// get rasterOrder_Rlang parameter
	std::string rasterOrder_Rlang = forestParameters.get_val<std::string>("rasterOrder_Rlang");
	std::transform(rasterOrder_Rlang.begin(), rasterOrder_Rlang.end(), rasterOrder_Rlang.begin(),
		[](unsigned char c){ return std::tolower(c); });

	if (rasterOrder_Rlang == "true")
		m_rasterOrder_Rlang = true;
	else
		m_rasterOrder_Rlang = false;

	// Get saveOnlyLast parameter
	std::string saveOnlyLast = forestParameters.get_val<std::string>("saveOnlyLast");
	std::transform(saveOnlyLast.begin(), saveOnlyLast.end(), saveOnlyLast.begin(),
		[](unsigned char c){ return std::tolower(c); });

	if (saveOnlyLast == "true")
	{
		m_saveOnlyLast = true;
		std::cout << "Only the last iteration will be saved, despite frequency = " << m_freqSave << std::endl;
	}
	else
	{
		m_saveOnlyLast = false;
		std::cout << "Populations' densities will be saved at a frequency = " << m_freqSave << std::endl;
	}

	if (m_freqSave > m_nIter && !m_saveOnlyLast)
		throw Except_Forest(m_freqSave, m_nIter, m_dim_land, false);

	// Check if last iteration is included when saving with a frequency
	m_lastIncludedInFreq = (((m_nIter - 1) % m_freqSave) == 0); // -1 because it goes from 0 to nIter - 1

	/**** Read landscape parameters from file climateFilename ****/
	par::Params climateParams(climateFilename.c_str(), "=");
	m_nRow_land = climateParams.get_val<unsigned int>("nRow");
	m_nCol_land = climateParams.get_val<unsigned int>("nCol");
	m_dim_land = m_nRow_land*m_nCol_land;
	m_deltaLon = climateParams.get_val<double>("deltaLon");
	m_deltaLat = climateParams.get_val<double>("deltaLat");
	std::string pathLandscape = climateParams.get_val<std::string>("path");
	checkPath(pathLandscape, "path (for landscape)");
	std::string distanceType = climateParams.get_val<std::string>("distance");

	std::string delimiter = climateParams.get_val<std::string>("delimiter");
	std::string climateFilenamePattern = climateParams.get_val<std::string>("filenamePattern");
	std::string isPopulated;

	if (m_nRow_land < 1 || m_nCol_land < 1)
		throw Except_Forest(m_nRow_land, m_nCol_land);

	/**** Creating species-specific folders (for output) and dispersal objects ****/
	// Checking and creating folder if necessary
	c_species_it species_it = m_speciesList.cbegin();
	std::string path_summary, path_popDyn;
	bool folderCreated;

	std::cout << "List of species:" << std::endl;
	for (; species_it != m_speciesList.cend(); ++species_it)
	{
		std::cout << "    - " << (*species_it)->m_speciesName << std::endl;

		// Checking summary's path and creating folder summary if necessary
		path_summary = m_summaryFilePath + (*species_it)->m_speciesName + "/";
		folderCreated = std::filesystem::create_directories(path_summary);
		if (folderCreated)
			std::cout << "Directory <" << path_summary << "> successfully created" << std::endl;

		// Checking popDyn's path and creating folder popDyn if necessary
		path_popDyn = m_popDynFilePath + (*species_it)->m_speciesName + "/";
		folderCreated = std::filesystem::create_directories(path_popDyn);
		if (folderCreated)
			std::cout << "Directory <" << path_popDyn << "> successfully created" << std::endl;

		// Creating dispersal object
		try
		{
			m_map_dispersal.emplace(*species_it, Dispersal(*species_it, climateFilename));

			if (m_map_dispersal.find(*species_it) != m_map_dispersal.end())
			{
				landscapeIntegrals(m_map_dispersal.at(*species_it));
				std::cout << m_map_dispersal.at(*species_it) << std::endl;
				std::cout << "------------------" << std::endl;
			}
		}
		catch(const std::exception& e)
		{
			std::cerr << e.what() << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	std::cout << "Species and dispersal object created" <<  std::endl;

	/**** Creating forest of patches ****/
	// Fill the environment
	std::string climateFile;
	std::string initFile;
	unsigned int counterPatch = 0;

	for (auto& p: std::filesystem::directory_iterator(pathLandscape))
	{
		climateFile = p.path().filename();
		if (climateFile.find(climateFilenamePattern) != std::string::npos)
		{
			climateFile = pathLandscape + climateFile;

			// Create environment
			Environment env(climateFile, delimiter, distanceType);
			
			// Precision set at 1e-3, should be fine for the numbers I am comparing: I neither expect microscopic nor gargantuan areas
			if (std::fabs(env.plotArea - m_deltaLon*m_deltaLat) > 1e-3)
				throw Except_Landscape(env.plotArea, m_deltaLon, m_deltaLat, climateFile);

			// Create Patch which initialise the populations
			try
			{
				m_patchVec.emplace_back(Patch(env, m_speciesList, m_initPath, m_initFilenamePattern, m_summaryFilePath, m_summaryFilePattern,
					m_popDynFilePath, m_popDynFilePattern, m_maxCohorts));
			}
			catch(const std::exception& e)
			{
				std::cerr << "An error occurred with the patch linked to <" << climateFile << ">" << "\n" << "Maybe check your initialisation path <" << m_initPath << ">" << e.what() << std::endl;
				exit(EXIT_FAILURE);
			}

			++counterPatch;
			if (counterPatch > m_dim_land)
				throw Except_Landscape(m_dim_land, climateFile);
		}
	}

	if (counterPatch != m_dim_land)
		throw Except_Landscape(m_dim_land, counterPatch);

	// Sort forest
	this->sort(true);

	std::cout << "Forest constructed with success, using following parameters " << std::endl << forestParameters << std::endl;
}

void Forest::patchDynamics(double const t, double const delta_t)
{
	std::for_each(std::execution::par_unseq, m_patchVec.begin(), m_patchVec.end(), [=] (Patch& patch){patch.populationDynamics(t, delta_t);});
	/*  Explanation lambda functions: https://docs.microsoft.com/en-us/cpp/cpp/lambda-expressions-in-cpp?view=vs-2019
		Square brackets [] are for the capture clause, i.e., how to access variables in the enclosing scope.
		In my case, the variables are t and delta_t. Because I put [=], I access them by value.
		If I had put [&], they were be accessed by reference. I could have put an hybrid such as [&t, delta_t]
	*/
}

void Forest::recruitment(double const t, double const delta_t)
{
	// Variables for neighbours
	std::vector<unsigned int> boundingBox;
	Patch* sourcePatch;

	// Iterator
	c_species_it sp_it;

	/*
	std::execution::seq (Williams2019, p. 330)
	The sequenced policy imposes few requirements on the iterators, values, and callable objects used with the algorithm:
	they may freely use synchronization mechanisms, and may rely on all operations being invoked on the same thread,
	though they cannot rely on the order of these operations.
	*/
	std::for_each(std::execution::seq, m_patchVec.begin(), m_patchVec.end(),
	[&] (Patch& patch){
		for (sp_it = m_speciesList.cbegin(); sp_it != m_speciesList.cend(); ++sp_it)
		{
			// Get neighbours, and store it in bounding box (using reference '&')
			neighbours_indices((patch.m_env).m_patchId, boundingBox, *sp_it); // boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};
			// Cover all the sources within neighbours to collect dispersed seeds
			for (unsigned int row = boundingBox[0]; row <= boundingBox[3]; ++row) // Important: less than or equal to (<=)
			{
				for (unsigned int col = boundingBox[1]; col <= boundingBox[2]; ++col) // Important: less than or equal to (<=)
				{
					sourcePatch = &m_patchVec[row*m_nCol_land + col];
					// Compute dispersal from source to target, and update the seed banks:
					patch.dispersal(sourcePatch, *sp_it, (m_map_dispersal.at(*sp_it)).m_map_distance_integral, m_deltaLat, m_deltaLon);
				}
			}
			// patch.recruitment(*sp_it, t, delta_t); // Compute recruitment for target patch and reset its seed bank
		}
	});
	std::for_each(std::execution::par_unseq, m_patchVec.begin(), m_patchVec.end(),
	[&] (Patch& patch)
	{
		for (sp_it = m_speciesList.cbegin(); sp_it != m_speciesList.cend(); ++sp_it)
			patch.recruitment(*sp_it, t, delta_t); // Compute recruitment for target patch and reset its seed bank
	});

}

void Forest::dynamics()
{
	// Time variables
	double t;
	double delta_t = (m_tmax - m_t0)/(m_nIter - 1);

	// Others
	std::string compReprodFilename;
	std::string popDynFilename;

	// Time loop
	if (!m_saveOnlyLast)
	{
		for (unsigned int iter = 1; iter < m_nIter; ++iter) // time loop, starts at 1 because the initial condition is considered the 0th iteration
		{
			t = m_t0 + (iter - 1)*delta_t; // iter starts at 1, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
			this->patchDynamics(t, delta_t);
			this->recruitment(t, delta_t);

			this->summary();

			if (iter % m_freqSave == 0)
				this->saveForest();
		}
		if (!m_lastIncludedInFreq)
			this->saveForest();
	}
	else
	{
		for (unsigned int iter = 1; iter < m_nIter; ++iter) // time loop, starts at 1 because the initial condition is considered the 0th iteration
		{
			t = m_t0 + (iter - 1)*delta_t; // iter starts at 1, but remember explicit Euler y_{n + 1} = y_n + delta_t f(t_n, y_n)
			this->patchDynamics(t, delta_t);
			this->recruitment(t, delta_t);
			this->summary();
		}
		this->saveForest();
	}

	std::cout << "Simulation done. Files saved in folders: <" << m_summaryFilePath << "*> and <" << m_popDynFilePath << "*>" << std::endl;
}

void Forest::neighbours_indices(unsigned int const target, std::vector<unsigned int>& boundingBox, Species const* species) const
{
	unsigned int topLeft_r, topLeft_c, topRight_c, bottomLeft_r;
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
		if (col_ind < influenceRadius_x) // To avoid overflows of unsigned int
			topLeft_c = 0;
		else
			topLeft_c = col_ind - influenceRadius_x;
		
		topRight_c = std::min(col_ind + influenceRadius_x, m_nCol_land);
	}

	if (influenceRadius_y >= m_deltaLat) // If more than itself is covered in latitude direction
	{
		if (row_ind < influenceRadius_y) // To avoid overflows of unsigned int
			topLeft_r = 0;
		else
			topLeft_r = row_ind - influenceRadius_y;
		
		bottomLeft_r = std::min(row_ind + influenceRadius_y, m_nRow_land);
	}

	// Particular cases: if max rows/cols reached, remove last row/col, because Forest::recruitment uses <= instead of <
	if (bottomLeft_r == m_nRow_land) // To prevent segmentation fault due to <= operator in Forest::recruitment
		--bottomLeft_r;
	if (topRight_c == m_nCol_land) // To prevent segmentation fault due to <= operator in Forest::recruitment
		--topRight_c;

	boundingBox = {topLeft_r, topLeft_c, topRight_c, bottomLeft_r};
}

/***********************************/
/******        Writing        ******/
/***********************************/
void Forest::saveForest() const
{
	std::for_each(std::execution::par_unseq, m_patchVec.cbegin(), m_patchVec.cend(),
		[=] (Patch const& patch){patch.savePatch();});
}

void Forest::summary() const
{
	std::for_each(std::execution::par_unseq, m_patchVec.cbegin(), m_patchVec.cend(),
		[=] (Patch const& patch){patch.summary();});
}

/************************************************/
/******        Sorting and managing        ******/
/************************************************/
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

void checkPath(std::string const& path, std::string const& key)
{
	if (path.empty())
		throw Except_Forest(path, key);

	if (path.back() != '/')
		throw Except_Forest(path, key);
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
