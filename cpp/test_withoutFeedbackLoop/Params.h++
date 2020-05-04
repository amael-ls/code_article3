
#ifndef PARAMS_H
#define PARAMS_H

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>

namespace par
{
	unsigned int spaceCounter(std::string const &str);

    std::vector<std::string> split(const std::string &s, const std::vector<char> &delim, const std::string &comment = "#");

    template<typename T>
    T str_convert(const std::string &s);

    template<> inline
	std::string str_convert<std::string>(const std::string &s);

	template <typename Type>
	unsigned int nbCompounds(Type number);

	template <> inline
	unsigned int nbCompounds(std::string str);

	void deleteSpace(std::string &str);

    class Params
    {
        std::map<std::string, std::vector<std::string> > data;
        const char *source;
        std::vector<char> delimiters;
        const std::string comment;

        void read_file();
        void get_lines(std::ifstream &file);

	public:
		Params(const char *f, const std::string &delim = " \t", const std::string &c = "#");

		template<typename T> inline
		std::vector<T> get_val(const std::string &key, bool kaamelott) const;

		std::string getKey(std::vector <std::string> const key) const;
		bool getKey(std::string const key) const;

		void printParams(std::ostream &os) const;

		bool isEmpty() const;
		void removeKey(std::string const &key);
    };


    // TEMPLATE FUNCTIONS
    template<typename T>
    std::vector<T> Params::get_val(const std::string &key, bool kaamelott) const
    {
		std::vector<T> result;
		try
		{
			std::vector<std::string> vals = data.at(key);
			if (!kaamelott)
			{
				for(unsigned int i = 0; i < vals.size(); i++)
					result.push_back(str_convert<T>(vals[i]));
			}
			else
			{
				for(unsigned int i = 0; i < vals.size(); i++)
				{
					deleteSpace(vals[i]);
					unsigned int j = 0;
					while (vals[i].size() != 0)
					{
						result.push_back(str_convert<T>(vals[i]));
						vals[i] = (vals[i]).substr(nbCompounds(result[j]), vals[i].size());
						deleteSpace(vals[i]);
						j++;
					}
				}
			}
		}
		catch(const std::out_of_range& ex)
		{
			std::stringstream ss;
			ss << "Warning: parameter parser tried to access unknown parameter <" << key << ">\t" << ex.what();
			throw (std::runtime_error (ss.str()));
		}

		return result;
    }



	template<typename T> inline
	T str_convert(const std::string &s)
	{
		T result;
		std::istringstream val(s); // create stream from the string
		if(!(val >> result))
		{
			std::stringstream ss;
			ss << "Cannot convert value <" << s << "> from string into requested type";
			throw( std::runtime_error (ss.str() ));
		}
		return result;
	}



	template<> inline
	std::string str_convert<std::string>(const std::string &s)
	{
//		std::cout << "Specialisation used here" << std::endl;
		return s;
	}



	template <typename Type>
	unsigned int nbCompounds(Type number)
	{
		int count = 0;
		if (number < 0)
			number = -number;
		do
		{
			number = number/10;
			count++;
		} while (number > 0);
		return count;
	}



	template <> inline
	unsigned int nbCompounds(std::string str)
	{
		return (str.length() - spaceCounter(str));
	}

	std::ostream &operator<<(std::ostream &os, par::Params const& params);
} // end namespace par

#endif /* defined(PARAMS_H) */
