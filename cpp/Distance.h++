
#ifndef DISTANCE_H
#define DISTANCE_H

// Official headers
#include <iostream>
#include <vector>
#include <limits> // for std::numeric_limits
#include <cmath>

class Distance
{
	friend class Dispersal;
	public:
		// latitude = row x deltaLat, longitude = col x deltaLon
		Distance(int const row1, int const col1, int const row2, int const col2, double const deltaLat, double deltaLon);

		// Overload
		friend bool operator<(Distance const d1, Distance const d2);
		friend bool operator==(Distance const d1, Distance const d2);
		friend std::ostream& operator<<(std::ostream& os, Distance const &dist);

	private:
		// double m_orthodromic;
		double m_euclidean;
		std::vector<int> m_manhattan;
};

#endif
