
#ifndef LANDSCAPE_C
#define LANDSCAPE_C

// Official headers
// #include <iomanip> // std::setw, std::left, std::setprecision
#include <filesystem> // To list files from folder, experimental/filesystem is now deprecated
#include <stdexcept>
#include <algorithm> // std::sort
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
	std::string rasterOrder_Rlang = metadata.get_val<std::string>("rasterOrder_Rlang");

	if (rasterOrder_Rlang == "true")
		m_rasterOrder_Rlang = true;
	
	// Fill the environment
	std::string filename;
	unsigned int counter = 0;

	for(auto& p: std::filesystem::directory_iterator(m_path))
	{
		filename = p.path().filename();
		if (filename.find(filenamePattern) != std::string::npos)
		{
			filename = m_path + "/" + filename;
			Environment* env = new Environment(filename, delimiter, counter);
			m_envVec.emplace_back(env);
			++counter;
			if (counter > m_dim)
				throw Except_Landscape(m_dim, filename);
		}
	}

	if (counter < m_dim)
		throw Except_Landscape(m_dim, counter);

	this->sort(m_rasterOrder_Rlang); // true to sort respecting raster order from R language

	// Compute Δlongitude and Δlatitude. No need to compute max and min of lon and lat: landscape is sorted
	// We assume the lattice is regular. Otherwise Δlongitude and Δlatitude should both be in Environment.
	std::vector<Environment*>::const_iterator it_first = m_envVec.cbegin();
	m_deltaLon = (*it_first)->distance(**std::next(it_first));
	if (std::next(it_first, m_nCol) < m_envVec.cend())
		m_deltaLat = (*it_first)->distance(**std::next(it_first, m_nCol));
	else
		m_deltaLat = 0;

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

std::ostream& operator<<(std::ostream& os, Landscape const &land)
{
	std::vector<Environment*>::const_iterator it = land.m_envVec.begin();
	for (; it != land.m_envVec.end(); ++it)
		(*it)->printCoordinates(os);
	return os;
}

/***********************************/
/******        Sorting        ******/
/***********************************/
void Landscape::sort(bool const rasterOrder_Rlang)
{
	// Sorting
	std::vector<Environment*>::iterator first = m_envVec.begin();
	std::vector<Environment*>::iterator last = m_envVec.end();
	if (rasterOrder_Rlang)
		std::sort(first, last, lessThan);
	else
		std::sort(first, last, greaterThan);
	
	// Renumbering
	int patch_id = 0;
	std::vector<Environment*>::iterator env_it;
	for (env_it = m_envVec.begin(); env_it != m_envVec.end(); ++env_it)
	{
		(*env_it)->m_patchId = patch_id;
		++patch_id;
	}
}

#endif
