
#ifndef ERROR_CLASSES_C
#define ERROR_CLASSES_C

#include "Error_classes.h++"

/**********************************/
/******        Forest        ******/
/**********************************/
Except_Forest::Except_Forest(unsigned int const freqSave, unsigned int const nIter, unsigned int const dimLandscape, bool const overPopulated):
	m_freqSave(freqSave), m_dimLandscape(dimLandscape), m_overPopulated(overPopulated)
{
	if (m_freqSave > nIter)
		m_error_msg += "the frequency of saving (" + std::to_string(m_freqSave) + ") is beyond tmax (" + std::to_string(nIter) + "). No output would be saved\n";
	
	if (m_overPopulated)
		m_error_msg += "Landscape is overpopulated. Dimension = " + std::to_string(m_dimLandscape) + " cells.\n";
}

Except_Forest::Except_Forest(unsigned int const nRow, unsigned int const nCol)
{
	if (nRow < 1)
		m_error_msg += "Wrong number of rows (" + std::to_string(nRow) + "). Should be at least 1 for longitude";
	if (nCol < 1)
		m_error_msg += "Wrong number of columns (" + std::to_string(nCol) + "). Should be at least 1 for latitude";
}

const char* Except_Forest::what() const throw()
{
	return (m_error_msg.c_str());
}

/*********************************/
/******        Patch        ******/
/*********************************/
Except_Patch::Except_Patch(unsigned int const patch_id, std::vector<std::string> const& speciesNames)
{
	m_error_msg += "No population file found despite the patch " + std::to_string(patch_id) + " is initially populated.\n";
	m_error_msg += "List of species provided:\n";
	for (unsigned int i = 0; i < speciesNames.size(); ++i)
		m_error_msg += "    - " + speciesNames[i] + "\n";
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

Except_Population::Except_Population(int const maxCohorts, int const nbCohorts, double const t):
	m_maxCohorts(maxCohorts), m_s_inf(-1)
{
	m_error_msg += "number of cohorts = " + std::to_string(nbCohorts) + " exceeds maximum number = " +
		std::to_string(m_maxCohorts) + " at time t = " + std::to_string(t)  + "\n";
}

const char* Except_Population::what() const throw()
{
	return (m_error_msg.c_str());
}

#endif