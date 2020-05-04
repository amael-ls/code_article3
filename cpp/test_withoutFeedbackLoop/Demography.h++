
#ifndef DEMOGRAPHY_H
#define DEMOGRAPHY_H

#include <cmath>

/*---	I chose c = 2	---*/

namespace demography
{
	double v(double const t, double const s, double const s_star)
	{
		return 2.0/(1 + s);
	}

	double d(double const t, double const s, double const s_star)
	{
		return 0.4; // (1 + s_star)*std::exp(- s_star);
	}

	double dv_ds(double const t, double const s, double const s_star)
	{
		return -2.0/((1 + s)*(1 + s));
	}

	double dd_ds(double const t, double const s, double const s_star)
	{
		return 0;
	}
}

#endif
