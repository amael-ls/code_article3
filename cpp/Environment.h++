
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
#include <string>
#include <map>

// My headers
#include "Params.h++"

// Forward declaration

class Environment
{
	friend class Population;
	friend class Landscape;
	friend class Forest;
	// friend class Dispersal;
	friend class Cohort;
	
	public :
		// Constructors
		Environment();
		Environment(std::string const filename, const std::string& delim);

		// Geography
		double distance(Environment const Env2) const;
		std::ostream& printCoordinates(std::ostream& os) const;

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Environment const &env);
		friend bool operator<(Environment const& env1, Environment const& env2);
		friend bool operator>(Environment const& env1, Environment const& env2);

		// Others
		void printId(std::ostream& os) const;

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
		double longitude;
		double latitude;
		std::string proj4string;
};

// External functions
double distancePoints(double longitude1, double latitude1, double longitude2, double latitude2);
bool lessThan(Environment* env1, Environment* env2);
bool greaterThan(Environment* env1, Environment* env2);

#endif
