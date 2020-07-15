
#ifndef LANDSCAPE_H
#define LANDSCAPE_H

// Official headers
#include <string>
#include <vector>

// My headers
#include "Environment.h++"

/***************************************/
/******      CLASS LANDSCAPE      ******/
/***************************************/

class Landscape
{
	friend class Population;
	friend class Dispersal;
	friend class Cohort;

	public :
		// Constructors
		Landscape(std::string const& metadataFile);

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Landscape const &land);
		Environment* operator[] (int const i);

		// Sorting
		void sort(bool const rasterOrder_Rlang);

	private :
		// Files' name
		std::string const m_filenamePattern;
		std::string const m_metadataFile;
		std::string m_path;

		// Dimensions
		unsigned int m_nRow, m_nCol, m_dim;

		// Spatial discretisation
		double m_deltaLon, m_deltaLat; // Correspond to Δx and Δy respectively

		// Vector of environment
		std::vector<Environment*> m_envVec;
};

#endif
