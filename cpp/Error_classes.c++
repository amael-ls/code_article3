
#ifndef ERROR_CLASSES_C
#define ERROR_CLASSES_C

#include "Error_classes.h++"

/**********************************/
/******        Forest        ******/
/**********************************/
Except_Forest::Except_Forest(unsigned int const freqSave, unsigned int const nIter):
	m_freqSave(freqSave)
{
	m_error_msg += "the frequency of saving (" + std::to_string(m_freqSave) + ") is beyond tmax (" + std::to_string(nIter) + "). No output would be saved\n";
}

const char* Except_Forest::what() const throw()
{
	return (m_error_msg.c_str());
}

/*************************************/
/******        Landscape        ******/
/*************************************/
Except_Landscape::Except_Landscape(int const dim, int const i):
	m_dimLandscape(dim)
{
	if (i < dim)
		m_error_msg += "not enough files provided. Length = " + std::to_string(m_dimLandscape) + " and number of read files = " + std::to_string(i) + "\n";
	
	if (i > dim)
		m_error_msg += "Length = " + std::to_string(m_dimLandscape) + " and index = " + std::to_string(i) + "\n";
}

Except_Landscape::Except_Landscape(int const dim, std::string const& filename):
	m_dimLandscape(dim)
{
	m_error_msg = "Error from Landscape: more files provided than dimension of the vector. Currently dim = " + std::to_string(m_dimLandscape) + ". ";
	m_error_msg = "Last file provided: " + filename + "\n";
}

const char* Except_Landscape::what() const throw()
{
	return (m_error_msg.c_str());
}

/**************************************/
/******        Population        ******/
/**************************************/
Except_Population::Except_Population(int const s_inf, int const tallestTree, std::string const& filename):
	m_maxCohorts(-1), m_s_inf(s_inf)
{
	m_error_msg += "tallest tree = " + std::to_string(tallestTree) + " from file = " + filename + " exceeds maximum size = " + std::to_string(m_s_inf) + "\n";
}

Except_Population::Except_Population(int const maxCohorts, std::string const& filename):
	m_maxCohorts(maxCohorts), m_s_inf(-1)
{
	m_error_msg += "number of cohorts from file = " + filename + " exceeds maximum number = " + std::to_string(m_maxCohorts) + "\n";
}

Except_Population::Except_Population(int const maxCohorts, int const nbCohorts):
	m_maxCohorts(maxCohorts), m_s_inf(-1)
{
	m_error_msg += "number of cohorts = " + std::to_string(nbCohorts) + " exceeds maximum number = " + std::to_string(m_maxCohorts) + "\n";
}

const char* Except_Population::what() const throw()
{
	return (m_error_msg.c_str());
}

#endif