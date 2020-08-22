
#ifndef ENVIRONMENT_C
#define ENVIRONMENT_C

// Official headers
#include <iomanip> // std::setw, std::left, std::setprecision
#include <limits> // for std::numeric_limits<double>::infinity()
#include <cmath> // for log, cos, sin, tg, M_PI

// My headers
#include "Environment.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Environment::Environment():
	m_fileName("")
{
	// Growth climate variables
	annual_mean_temperature = std::numeric_limits<double>::infinity();
	annual_precipitation = std::numeric_limits<double>::infinity();

	// Mortality cliPate variables
	min_temperature_of_coldest_month = std::numeric_limits<double>::infinity();
	precipitation_of_driest_quarter = std::numeric_limits<double>::infinity();

	// Plot area
	plotArea = std::numeric_limits<double>::infinity();

	// Spatial coordinates
	longitude = std::numeric_limits<double>::infinity();
	latitude = std::numeric_limits<double>::infinity();
	proj4string = "";
}

Environment::Environment(std::string const filename, const std::string& delim, unsigned int const patchId):
	m_patchId(patchId), m_fileName(filename)
{
	// Load parameters from files
	par::Params envParams(m_fileName.c_str(), delim, true);
	
	// Growth climate variables
	annual_mean_temperature = envParams.get_val<double>("annual_mean_temperature");
	annual_precipitation = envParams.get_val<double>("annual_precipitation");

	// Mortality cliPate variables
	min_temperature_of_coldest_month = envParams.get_val<double>("min_temperature_of_coldest_month");
	precipitation_of_driest_quarter = envParams.get_val<double>("precipitation_of_driest_quarter");

	// Plot area
	plotArea = envParams.get_val<double>("plotArea");

	// Spatial coordinates
	longitude = envParams.get_val<double>("longitude");
	latitude = envParams.get_val<double>("latitude");
	proj4string = envParams.get_val<std::string>("proj4string");
}

Environment::Environment(std::string const filename, const std::string& delim):
	Environment(filename, delim, 0) {}

/*************************************/
/******        Geography        ******/
/*************************************/
double distancePoints(double longitude1, double latitude1, double longitude2, double latitude2)
{
	/*
		It is assumed the longitudes and latitudes are in decimal degrees
	*/

	// Convert longitudes and latitudes to radians
	longitude1 *= M_PI/180;
	longitude2 *= M_PI/180;
	latitude1 *= M_PI/180;
	latitude2 *= M_PI/180;

	// Differences of lon
	double delta_lon = longitude2 - longitude1;
	double delta_lat = latitude2 - latitude1;

	// Compute dist in kilometers
	double radiusEarth = 6371;

	double haversine = sin(delta_lat/2)*sin(delta_lat/2) + cos(latitude1)*cos(latitude2)*sin(delta_lon/2)*sin(delta_lon/2);
	double angle = 2*atan2(sqrt(haversine), sqrt(1 - haversine));

	double dist = radiusEarth * angle;
	
	return dist;
}

double Environment::distance(Environment const Env2) const
{
	double dist = distancePoints(this->longitude, this->latitude, Env2.longitude, Env2.latitude);
	return dist;
}

std::ostream& Environment::printCoordinates(std::ostream& os) const
{
	os << this->m_patchId << "\t" << this->longitude << "\t" << this->latitude << std::endl;
	return os;
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream &operator<<(std::ostream &os, Environment const& env)
{
	os << std::setprecision(3);

	os << "Environment from file:" << std::endl;
	os << env.m_fileName << std::endl;
	os << std::endl;

	os << "Growth climate variables:" << std::endl;
	os << "T = " << env.annual_mean_temperature << " \t P = " << env.annual_precipitation << std::endl;
	os << std::endl;

	os << "Mortality climate variables:" << std::endl;
	os << "T = " << env.min_temperature_of_coldest_month << " \t P = " << env.precipitation_of_driest_quarter << std::endl;
	os << std::endl;

	os << "Plot area:" << std::endl;
	os << env.plotArea << std::endl;
	os << std::endl;

	os << "Longitude, latitude:" << std::endl;
	os << env.longitude << " \t " << env.latitude << std::endl;
	os << std::endl;

	os << "Proj4string:" << std::endl;
	os << env.proj4string << std::endl;
	os << std::endl;

	return os;
}

/*
	Remark on the way environment is sorted:
		- First by latitude, i.e., whatever the longitude, if an env is more to the south, it is smaller (i.e., the opposite direction of the latitude order)
		- Second, by longitude if same latitudes. In this case east is superior to west
	This is to keep the order R organise a raster. For example a 4 x 6 lattice is as follow:
		1  2  3  4
		5  6  7  8
		9  ...  12
		[...]   24
	We sort in the way to make the Id increasing when longitude increase and latitude decrease,
	where the Id is the index (from 1 to 24 in the R example. In C++ it would be shifted of -1)
*/
bool operator<(Environment const& env1, Environment const& env2)
{
	if (env1.latitude == env2.latitude)
		return (env1.longitude < env2.longitude);
	return (env1.latitude > env2.latitude); // Opposite direction of the latitude order, cf remark above
}

bool operator>(Environment const& env1, Environment const& env2)
{
	if (env1.latitude == env2.latitude)
		return (env1.longitude > env2.longitude);
	return (env1.latitude < env2.latitude); // Opposite direction of the latitude order, cf remark above
}

bool lessThan(Environment* env1, Environment* env2)
{
	return (*env1 < *env2);
}

bool greaterThan(Environment* env1, Environment* env2)
{
	return (*env1 > *env2);
}

#endif
