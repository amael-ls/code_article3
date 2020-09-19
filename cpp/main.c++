
#include <stdexcept> // exceptions (bad_alloc, bad_cast, out_of_range, ...)
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> // for high_resolution_clock
#include <cmath>

#include "Dispersal.h++"
#include "Forest.h++"
#include "Params.h++"

// #include "integration.h" // You need to copy paste all the src files from ALGLIB and then compile. Example:
// g++ -I. -std=c++17 -o demo.out *.cpp main.c++ Species.c++ Params.c++ Landscape.c++ Environment.c++ Error_classes.c++
// Linking the library won't work
// g++ -Wall -std=c++17 -I ~/Documents/Cpp_libraries/alglib_v3-16/src main.c++ -o test

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
	
	/***** Read simulation parameters and build species *****/
	// --- Load params
	par::Params const simulationParameters(argv[1], " = "); // Be extremely careful with the delimiter, especially white spaces
	
	std::string climate_file = simulationParameters.get_val<std::string>("climate_file");
	std::string species_path = simulationParameters.get_val<std::string>("species_path");
	
	// --- Create species
	std::vector<std::string> species = simulationParameters.get_val<std::vector<std::string> >("species_filenames");
	std::vector<std::string>::const_iterator species_filenames_it = species.cbegin();
	std::vector<Species*> speciesList;

	double t0 = simulationParameters.get_val<double>("t0");
	double tmax = simulationParameters.get_val<double>("tmax");
	double nIter = simulationParameters.get_val<double>("nIter");
	double delta_t = (tmax - t0)/(nIter - 1);

	for (; species_filenames_it != species.cend(); ++species_filenames_it)
	{
		Species* sp = new Species(*species_filenames_it, species_path, " = ", delta_t); // Be extremely careful with the delimiter, especially white spaces
		speciesList.emplace_back(sp);
	}

	/***** Build forest and run simulation *****/
	try
	{
		Forest test(simulationParameters, speciesList, climate_file);
		test.dynamics();
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	// auto finish = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> elapsed = finish - start;
	// std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
	return 0;
}
