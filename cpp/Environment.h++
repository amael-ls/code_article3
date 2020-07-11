
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

// My headers
#include "Params.h++"

class Environment
{
	friend class Cohort;
	friend class Population;
	public :
		// Constructors
		Environment();
		Environment(std::string const filename, const std::string& delim);

		// Geography
		double distance(Environment const Env2) const;

		// Overloading
		friend std::ostream& operator<<(std::ostream& os, Environment const &pop);

	private :
		// File name
		std::string m_fileName;

		// Growth climate variables
		double annual_mean_temperature;
		double annual_precipitation;

		// Mortality climate variables
		double min_temperature_of_coldest_month;
		double precipitation_of_driest_quarter;

		// Plot area
		double plotArea;

		// Spatial coordinates
		double longitude;
		double latitude;
		std::string proj4string;
};

// External functions
double distancePoints(double longitude1, double latitude1, double longitude2, double latitude2);

#endif
