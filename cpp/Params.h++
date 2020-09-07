
/*
	This class is adapted from Matthew Talutto.
*/

#ifndef PARAMS_H
#define PARAMS_H

// Official headers
#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>

namespace par {

	std::vector<std::string> split(const std::string &s, const std::vector<char> &delim, const std::string &comment = "#");
	std::vector<std::string> split(const std::string &s, std::string const &delim, const std::string &comment = "#");

	template<typename T>
	T str_convert(const std::string &s)
	{
		T result;
		std::istringstream val(s); // create stream from the string
		if(!(val >> result)) {
			std::stringstream ss;
			ss << "Cannot convert value <" << s << "> from string into requested type";
			throw( std::runtime_error (ss.str() ));
		}
		return result;
	}

	template<> inline
	std::string str_convert(const std::string &s)
	{
		return s;
	}

	template<> inline
	std::vector<std::string> str_convert(const std::string &s)
	{
		std::vector<std::string> results;
		results.push_back(s);
		return results;
	}

	class Params
	{
		std::map<std::string, std::vector<std::string> > data;
		const char * source;
		std::vector<char> delimiters;
		std::string delimiter_str;
		bool const m_strBool; // To know whether delimiters or delimiter_str should be used

		const std::string comment;

		void read_file();
		void get_lines(std::ifstream &file);

	public:
		Params(const char * f, const std::string &delim = " \t", bool const strBool = true, const std::string &c = "#");
		void printParams(std::ostream &os) const;

		template<typename T>
		T get_val(const std::string &key) const
		{
			T result;
			try {
				const std::vector<std::string> vals = data.at(key);
				result = str_convert<T>(vals[0]);
			}

			catch(const std::out_of_range& ex) {
				std::stringstream ss;
				ss << "Warning: parameter parser tried to access unknown parameter <" << key << ">\t" << ex.what();
				throw (std::runtime_error (ss.str()));
			}

			return result;
		}
	};	
} // end namespace par

// Specialising the template, the second separator is necessarily a comma followed by a space
template<> inline
std::vector<std::string> par::Params::get_val(const std::string &key) const
{
	const std::vector<std::string> vals = data.at(key);
	std::vector<std::string> lineData = split(vals[0], ", ");
	return lineData;
}

std::ostream &operator<<(std::ostream &os, par::Params const& params);

#endif /* defined(PARAMS_H) */
