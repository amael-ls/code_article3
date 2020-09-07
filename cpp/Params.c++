
#ifndef PARAMS_C
#define PARAMS_C

// Offical headers
#include <fstream>
#include <iomanip> // std::setw, std::left

// My headers
#include "Params.h++"

par::Params::Params(const char *f, const std::string &delim, bool const strBool, const std::string &c):
	source(f), delimiter_str(delim), m_strBool(strBool), comment(c)
{
	for(unsigned int i = 0; i < delim.length(); i++)
		delimiters.push_back(delim[i]);

	read_file();
}

void par::Params::read_file()
{
	std::ifstream inputFile (source);
	if(!inputFile.is_open())
	{
		std::stringstream ss;
		ss << "*** ERROR: could not open file <" << source << ">";
		throw (std::runtime_error (ss.str()));
	}
	get_lines(inputFile);
	inputFile.close();
}

void par::Params::get_lines(std::ifstream &file)
{
	while(file.good())
	{
		static size_t lineno = 0;	// for reporting the line number of errors
		std::string line;
		getline(file, line);
		lineno++;

		std::vector<std::string> lineData;
		std::string varName;
		try {
			if (!m_strBool)
				lineData = split(line, delimiters);
			else
				lineData = split(line, delimiter_str);
		}
		catch(...) {
			std::cerr << "error: problem parsing input, line " << lineno << std::endl;
			throw;
		}

		if(lineData.size() == 0)
			continue; // skip lines that are empty or were commented

		varName = lineData[0];
		lineData.erase(lineData.begin());

		data[varName] = lineData;
	}
}

std::vector<std::string> par::split(const std::string &s, const std::vector<char> &delim, const std::string &comment)
{
	std::vector<std::string> dest;

	// ignore commented lines
	if( s.substr(0,comment.length()) == comment )
		return dest;

	std::stringstream ls(s);	// create stringstream out of input
	std::string dat;
	while(getline(ls, dat, delim.back()))
	{
		if(delim.size() > 1)
		{
			std::vector<std::string> fillDest = split(dat, std::vector<char> (delim.begin(), delim.end() - 1), comment);
			dest.insert(dest.end(), fillDest.begin(), fillDest.end());
		}
		else
		{
			dest.push_back(dat);
		}
	}
	return dest;
}



std::vector<std::string> par::split(const std::string &s, std::string const &delim, const std::string &comment)
{
	std::vector<std::string> dest;

	// ignore commented lines
	if( s.substr(0,comment.length()) == comment )
		return dest;

	std::stringstream ls(s);	// create stringstream out of input
	std::string key, dat;
	size_t posDelim(0);

	// I assume only one information per line, otherwise while(pos != std::string::npos) is necessary
	getline(ls, dat);
	posDelim = dat.find_first_of(delim);
	if (posDelim != std::string::npos)
	{
		key = dat.substr(0, posDelim);
		dat = dat.substr(posDelim + delim.length(), std::string::npos);
		dest.push_back(key);
		dest.push_back(dat);
	}
	return dest;
}



void par::Params::printParams(std::ostream &os) const
{
	std::map<std::string, std::vector<std::string> >::const_iterator it = data.cbegin();
	std::vector<std::string>::const_iterator it_vec;
	for (; it != data.cend(); ++it)
	{
		os << std::right << std::setw(10) << "<" << it->first << "> :" << std::endl;
		for (it_vec = (it->second).cbegin(); it_vec != (it->second).cend(); ++it_vec)
			os << std::right << std::setw(40) << *it_vec << std::endl;
		os << std::endl;
	}
}



std::ostream &operator<<(std::ostream &os, par::Params const& params)
{
	params.printParams(os);
	return os;
}



// ========== CRASH TEST ZONE
// std::vector<std::string> par::Params::get_val(const std::string &key) const
// {
// 	std::vector<std::string> result;
// 	try
// 	{
// 		const std::vector<std::string> vals = data.at(key);
// 		std::cout << "lalala: " << vals[0] << std::endl;
// 		std::vector<std::string> lineData = split(vals[0], ", ");
// 		for(int i = 0; i < lineData.size(); i++)
// 		{
// 			std::cout << lineData[i] << std::endl;
// 			// result.push_back(vals[i]);
// 		}
// 	}
// 	catch(const std::out_of_range& ex)
// 	{
// 		std::stringstream ss;
// 		ss << "Warning: parameter parser tried to access unknown parameter <" << key << ">\t" << ex.what();
// 		throw (std::runtime_error (ss.str()));
// 	}

// 	return result;
// }


#endif /* defined(PARAMS_C) */
