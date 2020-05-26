
#ifndef ERROR_CLASSES_C
#define ERROR_CLASSES_C

#include "Error_classes.h++"

except_Landscape::except_Landscape(int const dim):
	m_dimLandscape(dim)
{
	m_error_msg = "Error from Landscape";
	if (dim < 0)
		m_error_msg += ": dimension must be positive. Currently dim = " + std::to_string(m_dimLandscape) + "\n";
}

except_Landscape::except_Landscape(int const dim, std::string const& filename):
	m_dimLandscape(dim)
{
	m_error_msg = "Error from Landscape: more files provided than dimension of the vector. Currently dim = " + std::to_string(m_dimLandscape) + ". ";
	m_error_msg = "Last file provided: " + filename + "\n";
}

const char* except_Landscape::what() const throw()
{
	return (m_error_msg.c_str());
}

#endif