
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
		void patchDynamics(double const t, double const delta_t);
		void spatialDynamics();

	// Return pointers to the neighbour cells
		void neighbours_indices(unsigned int const target, std::vector<int>& boundingBox, Species const* species) const;

	// // Print function
	// 	void print() const;

	// // Managing
		void sort(bool const rasterOrder_Rlang);
		// void addSpecies();
	
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

	// Dynamics parameters
		double  m_t0, m_tmax;
		unsigned int m_nIter;

	// // Reading parameters
	// 	std::string m_initFilenamePattern;
	// 	std::string m_initPath;

	// Saving options
		std::string m_compReprodFilePattern;
		std::string m_pathCompReprodFile;
		std::string m_popDynFilePattern;
		std::string m_pathPopDynFile;
		bool m_saveOnlyLast;
		unsigned int m_freqSave;

	// // Private functions
	// 	void createOutputFiles() const;
	// 	void close_ofs();
};

// Notes
// push_back dans fec, + order by Rlang

#endif
