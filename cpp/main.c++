
#include <stdexcept> // exceptions (bad_alloc, bad_cast, out_of_range, ...)
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> // for high_resolution_clock
#include <cmath>

#include "Population.h++"
// #include "Landscape.h++"
// #include "Forest.h++"
#include "Params.h++"

// #include "integration.h" // You need to copy paste all the src files from ALGLIB and then compile. Example:
// g++ -I. -std=c++17 -o demo.out *.cpp main.c++ Species.c++ Params.c++ Landscape.c++ Environment.c++ Error_classes.c++
// Linking the library won't work
// g++ -Wall -std=c++17 -I ~/Documents/Cpp_libraries/alglib_v3-16/src main.c++ -o test

/*
	In order to integrate the Kernel K, I had other choice to define a global variable and
	a wrapper function to do a callback using alglib::autogkintegrate
	If I could have modified the signature of the function, I would have done a safer wrapper,
	but alglib::autogkintegrate accepts only a on type of function which is:
		void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr)
	
	Here is the signature of alglib::autogkintegrate
		void autogkintegrate(autogkstate &state,
    		void (*func)(double x, double xminusa, double bminusx, double &y, void *ptr),
    		void *ptr = NULL, const xparams _xparams = alglib::xdefault);
*/
void* pt2Object; // global variable which points to an arbitrary Dispersal object for kernel integration

// Compiler command: g++ -Wall -std=c++17 *.cpp -lstdc++fs -o test

// CODER RK4, JE PENSE QU'EULER EST TROP INSTABLE
int main(int argc, char *argv[])
{
	// auto start = std::chrono::high_resolution_clock::now();
	std::cout << "Program " << argv[0] << " is running" << std::endl;
	std::cout << "Number of arguments: " << argc - 1 << std::endl;

	if (argc != 2)
	{
		std::cout << "ERROR (from main): 1 argument is accepted only" << std::endl;
		std::cout << "prototype: ./" << argv[0] << " simulationParameters.txt" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	par::Params simulationParameters(argv[1], ": "); // Be extremely careful with the delimiter, especially white spaces
	// std::cout << simulationParameters << std::endl;
	
	// unsigned int maxCohorts = simulationParameters.get_val<unsigned int>("maxCohorts");
	// double n_t = simulationParameters.get_val<double>("n_t");
	// double t0 = simulationParameters.get_val<double>("t0");
	// double t_max = simulationParameters.get_val<double>("t_max");
	std::string climate_file = simulationParameters.get_val<std::string>("climate_file");
	std::string species_filenames = simulationParameters.get_val<std::string>("species_filenames");
	std::string species_path = simulationParameters.get_val<std::string>("species_path");
	// std::string init_filenamePattern = simulationParameters.get_val<std::string>("init_filenamePattern");
	// std::string init_path = simulationParameters.get_val<std::string>("init_path");
	std::string forestDataFile = simulationParameters.get_val<std::string>("forestDataFile");
	
	Species* sp = new Species(species_filenames, species_path, " = "); // Be extremely careful with the delimiter, especially white spaces
	// std::cout << *sp << std::endl;

	// --- CRASH TEST ZONE TO REMOVE

	// Landscape *land = new Landscape(climate_file);
	// std::cout << *land << std::endl;

	// Forest forest(land, sp, forestDataFile);
	// forest.print();

	// forest.fct1();
	// forest.fct2();

	// forest.spatialDynamics();

	// double a = 0.84;
	// double b = 2.87;
	// double c = 5;
	// double d = 7;
	
	// double arrayParams[2] = {c, d};
	// double (*params)[2] = &arrayParams;

	// double value2d = 0;

	// Dispersal myDisp(sp, &land, land[0]);

	// pt2Object = (void*) &myDisp;

	// alglib::autogkstate ss;
	// alglib::autogkreport reprep;
	// alglib::autogksmooth(a, b, ss);
	// alglib::autogkintegrate(ss, Dispersal::wrapper_To_Call_Kintegral, params);
	// alglib::autogkresults(ss, value2d, reprep);

	// std::cout << value2d << std::endl;
	
	// --- END CRASH TEST ZONE

	// try
	// {
	// 	Landscape land(climate_file);
		// std::cout << *(land[0]) << std::endl;
		// std::cout << *(land[1]) << std::endl;

		// std::cout << land << std::endl;

		// std::cout << land[0]->distance(*land[1]) << std::endl;
		
		// try
		// {
		// 	Population pop2(maxCohorts, sp, init_filenamePattern, land[0], 0);
		
		// 	pop2.sort(true);
		// 	std::ofstream os_init("init.txt", std::ofstream::out);
		// 	std::ofstream os_end("end.txt", std::ofstream::out);
		// 	os_init << "iteration iterationBirth density dbh" << std::endl;
		// 	os_init << pop2;
		// 	try
		// 	{
		// 		pop2.euler(n_t, t0, t_max);
		// 	}
		// 	catch (const std::exception& e)
		// 	{
		// 		std::cerr << e.what() << '\n';
		// 		exit(EXIT_FAILURE);
		// 	}
		// 	os_end << "iteration iterationBirth density dbh" << std::endl;
		// 	os_end << pop2;
		// }
		// catch(const std::exception& e)
		// {
		// 	std::cerr << e.what() << '\n';
		// 	exit(EXIT_FAILURE);
		// }
	// }
	// catch(const std::exception& e)
	// {
	// 	std::cerr << e.what() << '\n';
	// 	exit(EXIT_FAILURE);
	// }
	
	// auto finish = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> elapsed = finish - start;
	// std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
	return 0;
}
