
#ifndef LANDSCAPE_C
#define LANDSCAPE_C

// Official headers
// #include <iomanip> // std::setw, std::left, std::setprecision
#include <experimental/filesystem> // To list files from folder
#include <iostream>

// My headers
#include "Landscape.h++"
#include "Params.h++"

// Shortcut namespace
namespace fs = std::experimental::filesystem;

/****************************************/
/******        Constructors        ******/
/****************************************/
Landscape::Landscape(std::string const& filenamePattern, std::string const& metadataFile):
	m_filenamePattern(filenamePattern), m_metadataFile(metadataFile)
{
	// Read metadata from file
	par::Params metadata(m_metadataFile.c_str(), "=");
	m_nRow = metadata.get_val<unsigned int>("nRow");
	m_nCol = metadata.get_val<unsigned int>("nCol");
	m_dim = m_nRow*m_nCol;
	m_path = metadata.get_val<std::string>("path");
	std::string delimiter = metadata.get_val<std::string>("delimiter");

	// Resize the vector of environment
	m_envVec.resize(m_dim);

	// Fill the environment
	std::string filename;
	std::vector<Environment>::iterator counter_it = m_envVec.begin();

	for(auto& p: fs::directory_iterator(m_path))
	{
		filename = p.path().filename();
		if (filename.find(filenamePattern) != std::string::npos)
		{
			m_envVec.emplace(counter_it, Environment(filename, delimiter));
			++counter_it;
		}
	}
}


#endif
