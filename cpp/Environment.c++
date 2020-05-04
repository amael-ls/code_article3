
#ifndef ENVIRONMENT_C
#define ENVIRONMENT_C

// Official headers
#include <iomanip> // std::setw, std::left, std::setprecision
#include <limits> // for std::numeric_limits<double>::infinity()
#include <cmath> // for log

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

Environment::Environment(std::string const filename, const std::string& delim):
	m_fileName(filename)
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

#endif
