
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "Population.h++"
#include "Cohort.h++"

// CODER RK4, JE PENSE QU'EULER EST TROP INSTABLE

int main(int argc, char *argv[])
{
	// std::cout << "Hello world" << std::endl;

	double s_inf = 5;
	double maxCohorts = 10;
	double delta_s = s_inf/maxCohorts;

	std::cout << "Δs = " << delta_s << std::endl;

	// Initial condition
	std::vector<double> l(1);
	l[0] = 1 + M_PI;
	std::vector<double> m(1);
	m[0] = M_PI;

	// for (unsigned int i = 0; i < 10; ++i) // Inititial condition
	// {
	// 	l[i] = 1 + i*delta_s; // 1 + s_i
	// 	m[i] = i*delta_s; // s_i
	// }

	// Population pop2(2000, 1);
	Population pop2(maxCohorts, s_inf, l, m); // Population(maxCohorts, s_inf, lambda, mu);
	pop2.sort(false);
	pop2.printNonZero();
	std::cout << pop2 << std::endl;
	std::ofstream os_init("init.txt", std::ofstream::out);
	std::ofstream os_end("end.txt", std::ofstream::out);
	os_init << pop2;
	pop2.sort(true);
	pop2.euler2(500, 0, 1); // euler2(n_t, t0, t_max)
	pop2.sort(false);
	os_end << pop2;
	pop2.printNonZero();

	return 0;
}


// #include <iostream>
// #include <boost/array.hpp>
//
// #include <boost/numeric/odeint.hpp>
//
// using namespace std;
// using namespace boost::numeric::odeint;
//
// const double sigma = 10.0;
// const double R = 28.0;
// const double b = 8.0 / 3.0;
//
// typedef boost::array< double , 3 > state_type;
//
// void lorenz( const state_type &x , state_type &dxdt , double t )
// {
//     dxdt[0] = sigma * ( x[1] - x[0] );
//     dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
//     dxdt[2] = -b * x[2] + x[0] * x[1];
// }
//
// void write_lorenz( const state_type &x , const double t )
// {
//     cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
// }
//
// int main(int argc, char **argv)
// {
//     state_type x = { 10.0 , 1.0 , 1.0 }; // initial conditions
//     integrate( lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz );
// }
