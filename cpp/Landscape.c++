
#ifndef LANDSCAPE_C
#define LANDSCAPE_C

// Official headers
// #include <iomanip> // std::setw, std::left, std::setprecision
#include <experimental/filesystem> // To list files from folder
#include <stdexcept>
#include <iostream>

// My headers
#include "Error_classes.h++"
#include "Landscape.h++"
#include "Params.h++"

// Shortcut namespace
namespace fs = std::experimental::filesystem;

/****************************************/
/******        Constructors        ******/
/****************************************/
Landscape::Landscape(std::string const& metadataFile):
	m_metadataFile(metadataFile)
{
	// Read metadata from file
	par::Params metadata(m_metadataFile.c_str(), "=");
	m_nRow = metadata.get_val<unsigned int>("nRow");
	m_nCol = metadata.get_val<unsigned int>("nCol");
	m_dim = m_nRow*m_nCol;
	m_path = metadata.get_val<std::string>("path");
	std::string delimiter = metadata.get_val<std::string>("delimiter");
	std::string filenamePattern = metadata.get_val<std::string>("filenamePattern");

	std::cout << metadata << std::endl;
	
	// Fill the environment
	std::string filename;
	unsigned int counter = 0;
	std::vector<Environment>::iterator counter_it = m_envVec.begin();

	for(auto& p: fs::directory_iterator(m_path))
	{
		filename = p.path().filename();
		if (filename.find(filenamePattern) != std::string::npos)
		{
			filename = m_path + "/" + filename;
			std::cout << filename << std::endl;
			m_envVec.emplace(counter_it, Environment(filename, delimiter));
			++counter_it;
			++counter;
			if (counter > m_dim)
				throw except_Landscape(m_dim, filename);
		}
	}

	if (counter < m_dim)
		throw except_Landscape(m_dim);
}


#endif
