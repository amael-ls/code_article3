
#include <stdexcept> // exceptions (bad_alloc, bad_cast, out_of_range, ...)
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> // for high_resolution_clock
#include <cmath>

#include "Population.h++"
#include "Landscape.h++"
#include "Params.h++"

// Compiler command: g++ -Wall -std=c++17 *.cpp -lstdc++fs -o test

// RK4
// RK4-5 adptatif

// CODER RK4, JE PENSE QU'EULER EST TROP INSTABLE
int main(int argc, char *argv[])
{
	auto start = std::chrono::high_resolution_clock::now();
	std::cout << "Program " << argv[0] << " is running" << std::endl;
	std::cout << "Number of arguments: " << argc - 1 << std::endl;
	
	if (argc != 2)
	{
		std::cout << "ERROR (from main): missing argument" << std::endl;
		std::cout << "prototype: ./" << argv[0] << " simulationParameters.txt" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	par::Params simulationParameters(argv[1], ": "); // Be extremely careful with the delimiter, especially white spaces
	std::cout << simulationParameters << std::endl;
	
	double maxCohorts = simulationParameters.get_val<double>("maxCohorts");
	double n_t = simulationParameters.get_val<double>("n_t");
	double t0 = simulationParameters.get_val<double>("t0");
	double t_max = simulationParameters.get_val<double>("t_max");
	std::string climate_file = simulationParameters.get_val<std::string>("climate_file");
	std::string species_filenames = simulationParameters.get_val<std::string>("species_filenames");
	std::string species_path = simulationParameters.get_val<std::string>("species_path");
	std::string init_file = simulationParameters.get_val<std::string>("init_file");
	
	Species* sp = new Species(species_filenames, species_path, " = "); // Be extremely careful with the delimiter, especially white spaces
	// std::cout << *sp << std::endl;

	try
	{
		Landscape land(climate_file);
		std::cout << *(land[0]) << std::endl;
		try
		{
			Population pop2(maxCohorts, sp, init_file, land[0], 0);
		
			pop2.sort(true);
			std::ofstream os_init("init.txt", std::ofstream::out);
			std::ofstream os_end("end.txt", std::ofstream::out);
			os_init << "iteration iterationBirth density dbh" << std::endl;
			os_init << pop2;
			try
			{
				pop2.euler(n_t, t0, t_max);
			}
			catch (const std::exception& e)
			{
				std::cerr << e.what() << '\n';
				exit(EXIT_FAILURE);
			}
			os_end << "iteration iterationBirth density dbh" << std::endl;
			os_end << pop2;
		}
		catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		exit(EXIT_FAILURE);
	}
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		exit(EXIT_FAILURE);
	}
	
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
	return 0;
}
