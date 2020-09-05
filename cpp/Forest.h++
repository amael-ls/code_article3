
#ifndef FOREST_H
#define FOREST_H

// Official headers
#include <vector>

// My headers
// #include "Landscape.h++"
#include "Species.h++"
#include "Patch.h++"

class Forest
{
	public :
		Forest(std::string const forestParamsFilename, std::vector<Species*> const speciesList, std::string const climateFilename);
		void dynamics();		
	
	// Overloading
		friend std::ostream& operator<<(std::ostream& os, Forest const &forest);

	private :
	// Forest parameters file
		std::string const m_forestParamsFilename;
		std::string m_initFilenamePattern;
		std::string m_initPath;

	// Landscape
	// --- Dimensions
		unsigned int m_nRow_land, m_nCol_land, m_dim_land; // Landscape

	// --- Discretisation
		double m_deltaLon, m_deltaLat; // Correspond to Δx and Δy respectively
	
	// --- Order; if true, then same order than a raster in R language
		bool m_rasterOrder_Rlang;
	
	// Species and population
		std::vector<Patch> m_patchVec;
		unsigned int m_maxCohorts;

		std::vector<Species*> m_speciesList;

	// Parameters related to dynamics
		double  m_t0, m_tmax;
		unsigned int m_nIter;
	
	// Functions related to dynamics
		void patchDynamics(double const t, double const delta_t);
		void recruitment(double const t, double const delta_t);

	// Bounding box of the neighbours of target patch
		void neighbours_indices(unsigned int const target, std::vector<unsigned int>& boundingBox, Species const* species) const;

	// Ordering
		void sort(bool const rasterOrder_Rlang);

	// Saving options
	// --- Variables
		std::string m_summaryFilePattern;
		std::string m_pathSummaryFile;
		std::string m_popDynFilePattern;
		std::string m_pathPopDynFile;
		bool m_saveOnlyLast;
		unsigned int m_freqSave;

	// --- Functions
		void saveResults() const; // To write results

	// // Private functions
	// 	void createOutputFiles() const;
	// 	void close_ofs();
};

// External functions
void checkPath(std::string const& path, std::string const& key);

#endif
