
/* Description of the class
The class Dispersal computes the integral of K for a given distance and saves it in a map. It contains:
	- A species
	- Landscape dimensions
	- A map <distance, integral(K)>, recording all possible distances in a landscape and associating the
		the surface integral of the dispersal kernel
	- Functions to compute the 2D-integral

Remarks:
* --- 1
In order to integrate the Kernel K, I had to define a global variable and a wrapper function
to do a callback using alglib::autogkintegrate

If I could have modified the signature of the function, I would have done a safer wrapper,
but alglib::autogkintegrate accepts only on type of function which is:
	void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr)
	
Here is the signature of alglib::autogkintegrate
	void autogkintegrate(autogkstate &state,
		void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr),
		void *ptr = NULL, const xparams _xparams = alglib::xdefault);

The 2D-integral has been checked on Matlab with the function:
	f(x, y) = x^2 * cos(sqrt(x - 3)) + 2 - cos(y) * sin(y) + exp(-0.5 * y^2);

Answer Matlab: -17.2171...
Answer C++   : -17.2171

* --- 2
Because the kernel K is only function of the distance, we can assume the targeted patch is located
in (0, 0), and the other patches are (0 + m * Δx, 0 + n * Δy), where m, n are integers between 0 and
m_nRow_land, m_nCol_land respectively
*/

#ifndef DISPERSAL_H
#define DISPERSAL_H

// Official headers
#include <vector>
#include <string>
#include <map>

// ALGLIB headers
#include "integration.h"

// My headers
#include "Error_classes.h++"
#include "Species.h++"

extern void* pt2Object; // global variable which points to an arbitrary Dispersal object for kernel integration

class Dispersal
{
	friend class Forest;
	public :
	// Constructors
		Dispersal(Species const* const sp, std::string const climateFilename);

	// Wrapper for Kernel integral computation (which are private functions)
		static void wrapper_r_integral(double x, double xminusa, double bminusx, double &y, void *ptr);
		static void wrapper_theta_integral(double x, double xminusa, double bminusx, double &y, void *ptr);
		friend void landscapeIntegrals(Dispersal& disp);

	// Overloading
		friend std::ostream& operator<<(std::ostream& os, Dispersal const &dispersal);

	private :
	// Species-specific data (cf Species.h++ for their definitions)
		Species const* const m_species;

		double const m_dispersalProbaThreshold;
		bool const m_min_dispersalProba;
		double const m_dispersalDistThreshold;
		bool const m_max_dispersalDist;

	// Landscape
	// --- Dimensions
		unsigned int m_nRow_land, m_nCol_land, m_dim_land; // Landscape

	// --- Discretisation
		double m_deltaLon, m_deltaLat; // Correspond to Δx and Δy respectively

	// Integration data
		std::map<double, double> m_map_distance_integral; // map<distance, integral(K)>, for each distance, compute the integral of K
		double m_totalIntegral;

	// private functions to compute integral
	void kernel(double x, double xminusa, double bminusx, double &y, void *ptr) const; // This function calls the species-specific kernel
	void theta_integral(double x, double xminusa, double bminusx, double &y, void *ptr);
	void r_integral(double x, double xminusa, double bminusx, double &y, void *ptr);
};

void landscapeIntegrals(Dispersal& disp);

#endif
