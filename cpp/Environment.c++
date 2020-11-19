
#ifndef ENVIRONMENT_C
#define ENVIRONMENT_C

// My headers
#include "Environment.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Environment::Environment():
	m_fileName(""), m_initPopulated(false), m_distance("")
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
	m_patchId = std::numeric_limits<unsigned int>::infinity();
	m_row = std::numeric_limits<int>::infinity();
	m_col = std::numeric_limits<int>::infinity();
	longitude = std::numeric_limits<double>::infinity();
	latitude = std::numeric_limits<double>::infinity();
	proj4string = "";
}

Environment::Environment(std::string const filename, const std::string& delim, std::string const distance):
	m_fileName(filename), m_distance(distance)
{
	// Load parameters from files
	par::Params envParams(m_fileName.c_str(), delim, true);
	
	// Growth climate variables
	annual_mean_temperature = envParams.get_val<double>("annual_mean_temperature");
	annual_precipitation = envParams.get_val<double>("annual_precipitation");

	// Mortality climate variables
	min_temperature_of_coldest_month = envParams.get_val<double>("min_temperature_of_coldest_month");
	precipitation_of_driest_quarter = envParams.get_val<double>("precipitation_of_driest_quarter");

	// Initially populated
	std::string isPopulated = envParams.get_val<std::string>("isPopulated");
	std::transform(isPopulated.begin(), isPopulated.end(), isPopulated.begin(),
		[](unsigned char c){ return std::tolower(c); });
	
	if (isPopulated == "true")
		m_initPopulated = true;
	else
		m_initPopulated = false;

	// Plot area
	plotArea = envParams.get_val<double>("plotArea");

	// Spatial coordinates
	m_patchId = envParams.get_val<double>("patch_id");
	m_row = envParams.get_val<unsigned int>("row");
	m_col = envParams.get_val<unsigned int>("col");
	longitude = envParams.get_val<double>("longitude");
	latitude = envParams.get_val<double>("latitude");
	proj4string = envParams.get_val<std::string>("proj4string");
}

/*************************************/
/******        Geography        ******/
/*************************************/
double distancePoints(double longitude1, double latitude1, double longitude2, double latitude2, std::string const distanceType)
{
	if (distanceType == "euclidean")
	{
		return std::sqrt((longitude2 - longitude1)*(longitude2 - longitude1) + (latitude2 - latitude1)*(latitude2 - latitude1));
	}

	/*
		It is assumed the longitudes and latitudes are in decimal degrees
	*/
	if (distanceType == "orthodromic")
	{
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

		return radiusEarth * angle;
	}

	return (-1);
}

double Environment::distance(Environment const Env2) const
{
	double dist = distancePoints(longitude, latitude, Env2.longitude, Env2.latitude, m_distance);
	return dist;
}

double distance(Environment const& env1, Environment const& env2)
{
	return distancePoints(env1.longitude, env1.latitude, env2.longitude, env2.latitude, env1.m_distance);
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
	os << "Patch id: " << env.m_patchId << std::endl;
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

void Environment::printId(std::ostream& os) const
{
	os << m_patchId;
}

unsigned int Environment::printId() const
{
	return m_patchId;
}
#endif
