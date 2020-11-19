
/* Description of the class
The class environment represents an ensemble of climatic conditions at a given
place. I organise it:
	- A file name (to get the data)
	- Climatic variables
	- Coordinates and projection string (in the proj4string format)

It contains:
	- annual_mean_temperature;
	- annual_precipitation;
	- min_temperature_of_coldest_month;
	- precipitation_of_driest_quarter;

I list the functions here, but describe them in the associated c++ file:
	- One constructors
	- Overloading for ease
*/

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

// Official headers
#include <algorithm> // for std::transform
#include <iomanip> // std::setw, std::left, std::setprecision
#include <limits> // for std::numeric_limits<double>::infinity()
#include <string>
#include <cmath> // for log, cos, sin, tg, M_PI
#include <map>

// My headers
#include "Params.h++"

// Forward declaration

class Environment
{
	friend class Distance;
	friend class Forest;
	friend class Cohort;
	friend class Patch;
	
	public :
		// Constructors
		Environment();
		Environment(std::string const filename, const std::string& delim, std::string const distance);

		// Geography
		double distance(Environment const Env2) const;
		friend double distance(Environment const& env1, Environment const& env2);
		std::ostream& printCoordinates(std::ostream& os) const;

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Environment const &env);
		friend bool operator<(Environment const& env1, Environment const& env2);
		friend bool operator>(Environment const& env1, Environment const& env2);

		// Others
		void printId(std::ostream& os) const;
		unsigned int printId() const;

	private :
		// File name
		std::string m_fileName;

		// Growth climate variables
		double annual_mean_temperature;
		double annual_precipitation;

		// Mortality climate variables
		double min_temperature_of_coldest_month;
		double precipitation_of_driest_quarter;

		// Initially populated
		bool m_initPopulated;

		// Plot area
		double plotArea;

		// Spatial coordinates
		unsigned int m_patchId;
		int m_row;
		int m_col;
		double longitude;
		double latitude;
		std::string proj4string;
		std::string m_distance;
};

// External functions
double distancePoints(double longitude1, double latitude1, double longitude2, double latitude2, std::string distanceType);
double distance(Environment const* const env1, Environment const* const env2);

#endif
