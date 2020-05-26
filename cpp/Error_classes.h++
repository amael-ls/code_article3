
#include <stdexcept>
#include <vector>
#include <string>

#include "Environment.h++"

#ifndef ERROR_CLASSES_H
#define ERROR_CLASSES_H

class except_Landscape : public std::exception
{
	public:
		except_Landscape(int const dim);
		except_Landscape(int const dim, std::string const& filename);
		except_Landscape(int const dim, std::vector<Environment>::iterator const current, std::vector<Environment>::iterator const endVec);
		const char* what() const throw();

	private:
		int m_dimLandscape;
		std::string m_error_msg;
};

#endif
