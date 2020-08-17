
/* Description of the class
The class Dispersal represents the set of patches accessible from a population at a
given place. I organise it:
	- A vector of unsigned integers storing the indices of patches reachable by the population
	- A vector of doubles computing the proportion seeds dispersed 

I list the functions here, but describe them in the associated c++ file:
	- A constructor
*/

#ifndef DISPERSAL_H
#define DISPERSAL_H

// Official headers
#include <vector>

// My headers
#include "Landscape.h++"
#include "Species.h++"

/*
	In order to integrate the Kernel K, I had other choice to define a global variable and
	a wrapper function to do a callback using alglib::autogkintegrate
	If I could have modified the signature of the function, I would have done a safer wrapper,
	but alglib::autogkintegrate accepts only on type of function which is:
		void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr)
	
	Here is the signature of alglib::autogkintegrate
		void autogkintegrate(autogkstate &state,
    		void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr),
    		void *ptr = NULL, const xparams _xparams = alglib::xdefault);
*/
extern void* pt2Object; // global variable which points to an arbitrary Dispersal object for kernel integration

class Dispersal
{
	friend class Population;

	public :
		// Constructors
		Dispersal(Species const* const sp, Landscape const* const land, Environment const* const popEnv);

		// Wrapper for Kernel integral computation (which are private functions)
		static void wrapper_To_Call_Kintegral(double x, double xminusa, double bminusx, double &y, void *ptr);
		static void wrapper_To_Call_Kintegral_lon(double x, double xminusa, double bminusx, double &y, void *ptr);

	private :
	Species const* const m_species;
	Landscape const* const m_landscape;
	Environment const* const m_popEnv;
	std::vector<unsigned int> m_indices;
	std::vector<double> m_proportions;

	// private functions to compute integral
	void Kintegrand_lon(double x, double xminusa, double bminusx, double &y, void *ptr) const;
	void Kintegral_lon(double z, double xminusa, double bminusx, double &y, void *ptr) const;
};

#endif
