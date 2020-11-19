
#include <stdexcept> // exceptions (bad_alloc, bad_cast, out_of_range, ...)
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> // for high_resolution_clock
#include <cmath>

#include "Dispersal.h++"
#include "Forest.h++"
#include "Params.h++"

// Compilation on the super computer; the -ltbb flag might be required (not with gcc/10.2.0):
// g++ -I. -std=c++2a -Wall -o demo.out main.c++ Cohort.c++ Species.c++ Params.c++ Environment.c++ Error_classes.c++ Population.c++ Patch.c++ Forest.c++ Dispersal.c++ Distance.c++ *.cpp

void* pt2Object; // global variable which points to an arbitrary Dispersal object for kernel integration

// Compiler command: g++ -Wall -std=c++17 *.cpp -lstdc++fs -o test

int main(int argc, char *argv[])
{
	auto start = std::chrono::high_resolution_clock::now();
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

	for (; species_filenames_it != species.cend(); ++species_filenames_it)
	{
		Species* sp = new Species(*species_filenames_it, species_path, " = "); // Be extremely careful with the delimiter, especially white spaces
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

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
	return 0;
}
