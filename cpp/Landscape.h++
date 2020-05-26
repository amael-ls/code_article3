
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
	friend class Cohort;
	friend class Population;
	public :
		// constructors
		Landscape(std::string const& metadataFile);

		// Overloading
		Environment& operator[] (int const i);
		Environment* operator() (int const i);

	private :
		// Files' name
		std::string const m_filenamePattern;
		std::string const m_metadataFile;
		std::string m_path;

		// Dimensions
		unsigned int m_nRow, m_nCol, m_dim;

		// Vector of environment
		std::vector<Environment> m_envVec;
};

#endif
