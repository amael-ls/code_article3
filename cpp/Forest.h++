
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
		Forest(Landscape* land, Species* sp, std::string const forestParamsFilename);
		void spatialDynamics();

	// Return pointers to the neighbour cells
		void neighbours_indices(unsigned int const target, std::vector<int>& boundingBox) const;

	// Print function
		void print() const;

	// Sorting function
		void sort(bool const rasterOrder_Rlang);

	// Debug Fct
		void fct1(); // Print address map isPresent;

	private :
	// Forest parameters file
		std::string const m_forestParamsFilename;

	// Landscape
		Landscape* m_land;
		unsigned int m_nRow_land, m_nCol_land, m_dim_land; // Landscape dimensions
	
	// Order (see Landscape for more details)
		bool m_rasterOrder_Rlang; // if true, then same order than a raster in R language
	
	// Species and population
		std::vector<Population> m_popVec; // size(m_popVec) <= size(m_land)
		unsigned int m_maxCohorts;

		Species* m_sp;

	// Dynamics parameters
		double  m_t0, m_tmax;
		unsigned int m_nIter;

	// Reading parameters
		std::string m_initFilenamePattern;
		std::string m_initPath;

	// Saving options
		std::string m_compReprodFilePattern;
		std::string m_pathCompReprodFile;
		std::string m_popDynFilePattern;
		std::string m_pathPopDynFile;
		bool m_saveOnlyLast;
		unsigned int m_freqSave;

	// Private functions
		void createOutputFiles() const;
		void close_ofs();
};

// Notes
// push_back dans fec, + order by Rlang

#endif
