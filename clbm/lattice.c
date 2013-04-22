#include "lattice.h"
#include <math.h>

#ifdef D2Q9
	double cx[9] = {0, 1.0, 0, -1.0, 0, 1.0, -1.0, -1.0, 1.0};
	double cy[9] = {0, 0, 1.0, 0, -1.0, 1.0, 1.0, -1.0 , -1.0};
	double weight[9] = {4.0/9.0,
			1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
			1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

	double cs =  0.577350269189625842081;
	double cs2 = 0.333333333333333333333;

	double feq(unsigned int k, double rho0, double ux0, double uy0)
	{
		double cu = 3*(cx[k]*ux0 + cy[k]*uy0);
		return rho0 * weight[k] *
				(1 + cu + 0.5*cu*cu - (3.0/2.0)*(ux0*ux0 + uy0*uy0));
	}
#endif
