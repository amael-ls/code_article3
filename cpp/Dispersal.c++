
#ifndef DISPERSAL_C
#define DISPERSAL_C

// My headers
#include "Dispersal.h++"

/****************************************/
/******        Constructors        ******/
/****************************************/
Dispersal::Dispersal(Species const* const sp, std::string const climateFilename):
	m_species(sp), m_dispersalProbaThreshold(sp->dispersalProbaThreshold), m_min_dispersalProba(sp->min_dispersalProba),
	m_dispersalDistThreshold(sp->dispersalDistThreshold), m_max_dispersalDist(sp->max_dispersalDist), m_totalIntegral(0)
{
	/**** Read landscape parameters from file climateFilename ****/
	par::Params climateParams(climateFilename.c_str(), "=");

	m_nRow_land = climateParams.get_val<unsigned int>("nRow");
	m_nCol_land = climateParams.get_val<unsigned int>("nCol");
	m_dim_land = m_nRow_land*m_nCol_land;

	m_deltaLon = climateParams.get_val<double>("deltaLon");
	m_deltaLat = climateParams.get_val<double>("deltaLat");
}

/*************************************************/
/******        Dispersal integration        ******/
/*************************************************/
void Dispersal::kernel(double r, double xminusa, double bminusx, double &y, void *ptr) const
{
	// double *param = (double *) ptr; // casting
	// double const z = param[0];

	// y = m_species->K(x, z); //y = x*x*cos(sqrt(x - 3)) + 2 - cos(z)*sin(z) + exp(-z*z/2.0);
	y = r * m_species->K(r); // Because integrate g(r) r drdθ // ! Using polar coords
}

void Dispersal::wrapper_r_integral(double x, double xminusa, double bminusx, double &y, void *ptr)
{
	// explicitly cast global variable <pt2Object> to a pointer to Dispersal
	// warning: <pt2Object> MUST point to an appropriate object!
	Dispersal* mySelf = (Dispersal*) pt2Object;

	// call member
	mySelf->kernel(x, xminusa, bminusx, y, ptr);
}

void Dispersal::theta_integral(double theta, double xminusa, double bminusx, double &y, void *ptr)
{
	// Access parameters from void *ptr
	// --- Theta bounds
	double **boundsInteg = (double **) ptr;
	double *bounds = boundsInteg[0];
	
	double x1 = bounds[0];
	double x2 = bounds[1];
	double y1 = bounds[2];
	double y2 = bounds[3];

	// --- Case and zone
	std::string **strInteg = (std::string **) ptr;
	std::string *param_str = strInteg[1];
	std::string ratioCase = param_str[0];
	std::string zone = param_str[1];

	// Declare and compute r1 and r2, the radius boundaries of r_integral
	double r1 = 0;
	double r2 = 0;

	// Integral bounds for r_integral
	if (ratioCase == "y2/x2 <= y1/x1") // --- Case y2/x2 <= y1/x1
	{
		if (zone == "zone1")
		{
			r1 = y1/std::sin(theta);
			r2 = x2/std::cos(theta);
		}
		else if (zone == "zone2")
		{
			r1 = y1/std::sin(theta);
			r2 = y2/std::sin(theta);
		}
		else if (zone == "zone3")
		{
			r1 = x1/std::cos(theta);
			r2 = y2/std::sin(theta);
		}
	}
	else if (ratioCase == "y1/x1 < y2/x2") // --- Case y1/x1 < y2/x2
	{
		if (zone == "zone1")
		{
			r1 = y1/std::sin(theta);
			r2 = x2/std::cos(theta);
		}
		else if (zone == "zone2")
		{
			r1 = x1/std::cos(theta);
			r2 = x2/std::cos(theta);
		}
		else if (zone == "zone3")
		{
			r1 = x1/std::cos(theta);
			r2 = y2/std::sin(theta);
		}
	}
	else if (ratioCase == "pathos_x1")
	{
		if (zone == "zone01_1")
		{
			r1 = y1/std::sin(theta);
			r2 = x2/std::cos(theta);
		}
		else // zone = zone01_2
		{
			r1 = y1/std::sin(theta);
			r2 = y2/std::sin(theta);
		}
	}
	else if (ratioCase == "pathos_y1")
	{
		if (zone == "zone10_1")
		{
			r1 = x1/std::cos(theta);
			r2 = x2/std::cos(theta);
		}
		else // zone = zone10_2
		{
			r1 = x2/std::cos(theta);
			r2 = y2/std::sin(theta);
		}
	}
	else // pathos_origin
	{
		r1 = 0;
		if (zone == "zone00_1")
			r2 = x2/std::cos(theta);
		else // zone = zone00_2
			r2 = y2/std::sin(theta);
	}

	alglib::autogkstate s;
	alglib::autogkreport rep;
	alglib::autogksmooth(r1, r2, s);
	alglib::autogkintegrate(s, Dispersal::wrapper_r_integral, &y);
	alglib::autogkresults(s, y, rep);
}

void Dispersal::wrapper_theta_integral(double theta, double xminusa, double bminusx, double &y, void *ptr)
{
	// explicitly cast global variable <pt2Object> to a pointer to Dispersal
	// warning: <pt2Object> MUST point to an appropriate object!
	Dispersal* mySelf = (Dispersal*) pt2Object;

	// call member
	mySelf->theta_integral(theta, xminusa, bminusx, y, ptr);
}

void landscapeIntegrals(Dispersal& disp)
{
	// Integral bounds
	// --- Cartesian
	double x1, x2;
	double y1, y2;

	// --- Polar (angle theta only)
	double theta1, theta2;

	// Dimension patch (rectangle x1,y1 ---> x2,y2, with x2 = x1 + deltaLon and y2 = y1 + deltaLat)
	double deltaLon(disp.m_deltaLon);
	double deltaLat(disp.m_deltaLat);
	double euclideanDistanceToZero;

	// Array and pointer (parameters to integrate in Cartesian a Polar function)
	double coordsPatch[4]; // {x1, x2, y1, y2};
	std::string caseZone[2]; // {"case", "zone"}

	void *params[2];
	params[0] = &coordsPatch;
	params[1] = &caseZone;

	// Values for integral
	pt2Object = (void*) &disp;

	double value2d(0), integPatch(0);
	disp.m_totalIntegral = 0;

	alglib::autogkstate ss;
	alglib::autogkreport reprep;

	for (unsigned int row = 0; row < disp.m_nRow_land; ++row) // latitude direction
	{
		for (unsigned int col = 0; col < disp.m_nCol_land; ++col) // longitude direction
		{
			// Reset integPatch for current patch
			integPatch = 0;

			// Coordinates and distance, row = y = latitude, and col = x = longitude
			x1 = col*deltaLon;
			y1 = row*deltaLat;
			Distance distanceToZero(0, 0, row, col, deltaLat, deltaLon);
			euclideanDistanceToZero = sqrt(x1*x1 + y1*y1);

			x2 = x1 + deltaLon;
			y2 = y1 + deltaLat;

			coordsPatch[0] = x1;
			coordsPatch[1] = x2;
			coordsPatch[2] = y1;
			coordsPatch[3] = y2;

			// There are three pathological cases with two zones, and two non pathological cases with three zones to differentiate
			if (y2/x2 <= y1/x1) // --- Non pathological case 1: y2/x2 <= y1/x1
			{
				caseZone[0] = "y2/x2 <= y1/x1";

				// --- First zone, case 2
				caseZone[1] = "zone1";

				theta1 = std::atan(y1/x2);
				theta2 = std::atan(y2/x2);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;

				// --- Second zone, case 2
				caseZone[1] = "zone2";

				theta1 = std::atan(y2/x2);
				theta2 = std::atan(y1/x1);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;

				// --- Third zone, case 2
				caseZone[1] = "zone3";

				theta1 = std::atan(y1/x1);
				theta2 = std::atan(y2/x1);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;
			}
			else if (y1/x1 < y2/x2) // --- Non pathological case 2: y1/x1 < y2/x2
			{
				caseZone[0] = "y1/x1 < y2/x2";

				// --- First zone, case 1
				caseZone[1] = "zone1";

				theta1 = std::atan(y1/x2);
				theta2 = std::atan(y1/x1);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;

				// --- Second zone, case 1
				caseZone[1] = "zone2";

				theta1 = std::atan(y1/x1);
				theta2 = std::atan(y2/x2);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;

				// --- Third zone, case 1
				caseZone[1] = "zone3";

				theta1 = std::atan(y2/x2);
				theta2 = std::atan(y2/x1);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;
			}
			else if ((x1 == 0) && (y1 != 0)) // Pathological case 1: x1 = 0
			{
				caseZone[0] = "pathos_x1";

				// --- First zone, pathology x1 = 0
				caseZone[1] = "zone01_1";
				theta1 = std::atan(y1/x2);
				theta2 = std::atan(y2/x2);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;

				// --- Second zone, pathology x1 = 0
				caseZone[1] = "zone01_2";
				theta1 = std::atan(y2/x2);
				theta2 = M_PI_2; // pi/2

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;
			}
			else if ((x1 != 0) && (y1 == 0)) // Pathological case2: y1 = 0, i.e., targetPatch = sourcePatch
			{
				caseZone[0] = "pathos_y1";

				// --- First zone, pathology y1 = 0
				caseZone[1] = "zone10_1";
				theta1 = 0;
				theta2 = std::atan(y2/x2);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;

				// --- Second zone, pathology y1 = 0
				caseZone[1] = "zone10_2";
				theta1 = std::atan(y2/x2);
				theta2 = std::atan(y2/x1);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;
			}
			else // Pathological case 3: (x1, y1) = (0, 0), i.e., targetPatch = sourcePatch
			{
				caseZone[0] = "pathos_origin";

				// --- First zone, pathology at origin
				caseZone[1] = "zone00_1";
				theta1 = 0;
				theta2 = std::atan(y2/x2);

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;

				// --- Second zone, pathology at origin
				caseZone[1] = "zone00_2";
				theta1 = std::atan(y2/x2);
				theta2 = M_PI_2; // pi/2

				alglib::autogksmooth(theta1, theta2, ss);
				alglib::autogkintegrate(ss, Dispersal::wrapper_theta_integral, params);
				alglib::autogkresults(ss, value2d, reprep);

				integPatch += value2d;
			}
			
			disp.m_totalIntegral += integPatch;

			if ((disp.m_dispersalProbaThreshold) &&(integPatch >= disp.m_min_dispersalProba)) // If using proba threshold
			{
				if (disp.m_map_distance_integral.find(distanceToZero) == disp.m_map_distance_integral.end())
					disp.m_map_distance_integral[distanceToZero] = integPatch;
			}

			if ((disp.m_dispersalDistThreshold) && (euclideanDistanceToZero <= disp.m_dispersalDistThreshold)) // If using distance threshold
			{
				if (disp.m_map_distance_integral.find(distanceToZero) == disp.m_map_distance_integral.end())
					disp.m_map_distance_integral[distanceToZero] = integPatch;
			}
		}
	}

	// Because we integrated only on a quarter of the plane, we need to multiply the total integral by 4
	disp.m_totalIntegral *= 4; // If the landscape is big enough, this should be close to 1

	// std::map<double, double>::const_iterator it = disp.m_map_distance_integral.cbegin();
	// for (; it != disp.m_map_distance_integral.cend(); ++it)
	// 	std::cout << it->first << "    " << it->second << std::endl;
	// std::cout << std::endl;

	if ((disp.m_totalIntegral > 1) || (disp.m_totalIntegral < 0))
		throw Except_Dispersal(disp.m_totalIntegral, disp.m_species->getName());
}

/************************************/
/******        Overload        ******/
/************************************/
std::ostream& operator<<(std::ostream& os, Dispersal const &dispersal)
{
	// Landscape
	os << std::endl << "Dimensions (row x col): " << dispersal.m_nRow_land << " x " << dispersal.m_nCol_land << std::endl;
	os << "Resolution (lon x lat): " << dispersal.m_deltaLon << " x " << dispersal.m_deltaLat << std::endl;
	
	os << "Integral on Γ: " << dispersal.m_totalIntegral << std::endl;

	os << "Dimension of the dispersal map: " << (dispersal.m_map_distance_integral).size();
	return os;
}

#endif
