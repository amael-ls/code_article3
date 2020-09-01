
#include <stdexcept>
#include <vector>
#include <string>

#include "Environment.h++"

#ifndef ERROR_CLASSES_H
#define ERROR_CLASSES_H

/**********************************/
/******        Forest        ******/
/**********************************/
// class Except_Forest : public std::exception
// {
// 	public:
// 		Except_Forest(unsigned int const freqSave, unsigned int const nIter, unsigned int const dimLandscape, bool const overPopulated);
// 		const char* what() const throw();

// 	private:
// 		std::string m_error_msg = "Error from Forest: ";
// 		unsigned int const m_freqSave;
// 		unsigned int const m_dimLandscape;
// 		bool const m_overPopulated;
// };

/*************************************/
/******        Landscape        ******/
/*************************************/
// class Except_Landscape : public std::exception
// {
// 	public:
// 		Except_Landscape(int const dim);
// 		Except_Landscape(int const dim, int const i);
// 		Except_Landscape(int const dim, std::string const& filename);
// 		const char* what() const throw();

// 	private:
// 		int m_dimLandscape;
// 		std::string m_error_msg = "Error from Landscape: ";
// };

/**************************************/
/******        Population        ******/
/**************************************/
class Except_Population : public std::exception
{
	public:
		Except_Population(int const s_inf, int const tallestTree, std::string const& filename);
		Except_Population(int const maxCohorts, std::string const& filename);
		Except_Population(int const maxCohorts, int const nbCohorts, double const t);
		Except_Population(int const maxCohorts, int const nbCohorts);
		const char* what() const throw();

	private:
		int m_maxCohorts;
		int m_s_inf;
		std::string m_error_msg = "Error from Population: ";
};

#endif
