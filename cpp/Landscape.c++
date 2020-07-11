
#ifndef LANDSCAPE_C
#define LANDSCAPE_C

// Official headers
// #include <iomanip> // std::setw, std::left, std::setprecision
#include <filesystem> // To list files from folder, experimental/filesystem is now deprecated
#include <stdexcept>
#include <iostream>

// My headers
#include "Error_classes.h++"
#include "Landscape.h++"
#include "Params.h++"

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
	
	// Fill the environment
	std::string filename;
	unsigned int counter = 0;
	// std::vector<Environment*>::iterator counter_it = m_envVec.begin();

	for(auto& p: std::filesystem::directory_iterator(m_path))
	{
		filename = p.path().filename();
		if (filename.find(filenamePattern) != std::string::npos)
		{
			filename = m_path + "/" + filename;
			std::cout << filename << std::endl;
			Environment* env = new Environment(filename, delimiter);
			m_envVec.emplace_back(env);
			// ++counter_it;
			++counter;
			if (counter > m_dim)
				throw Except_Landscape(m_dim, filename);
		}
	}

	if (counter < m_dim)
		throw Except_Landscape(m_dim, counter);

	std::cout << "Landscape constructed with success" << std::endl;
}

/************************************/
/******        Overload        ******/
/************************************/
Environment* Landscape::operator[] (int const i)
{
	if (i > m_dim)
		throw Except_Landscape(m_dim, i);

	return m_envVec[i];
}

#endif
