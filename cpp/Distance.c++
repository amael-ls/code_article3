
#ifndef DISTANCE_C
#define DISTANCE_C

#include "Distance.h++"

/***************************************/
/******        Constructor        ******/
/***************************************/
// latitude = row x deltaLat, longitude = col x deltaLon
Distance::Distance(int const row1, int const col1, int const row2, int const col2, double const deltaLat, double deltaLon):
	m_manhattan(2, -std::numeric_limits<int>::infinity())
{
	m_manhattan[0] = std::abs(col1 - col2); // x-direction (longitude)
	m_manhattan[1] = std::abs(row1 - row2); // y-direction (latitude)

	double x1 = col1*deltaLon;
	double x2 = col2*deltaLon;
	double y1 = row1*deltaLat;
	double y2 = row2*deltaLat;

	m_euclidean = std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

/************************************/
/******        Overload        ******/
/************************************/
bool operator<(Distance const d1, Distance const d2)
{
	if (d1.m_manhattan[1] == d2.m_manhattan[1]) // if (same latitude)
		return (d1.m_manhattan[0] < d2.m_manhattan[0]); // compare longitude
	
	return d1.m_manhattan[1] < d2.m_manhattan[1];
}

bool operator==(Distance const d1, Distance const d2)
{
	return (d1.m_manhattan == d2.m_manhattan); // This operator is overloaded for vectors
}

std::ostream& operator<<(std::ostream& os, Distance const &dist)
{
	os << "Manhattan: " << dist.m_manhattan[0] << " Δx + " << dist.m_manhattan[1] << " Δy" << std::endl;
	os << "Euclidean: " << dist.m_euclidean;
	return os;
}

#endif
