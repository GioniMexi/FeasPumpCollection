/**
 * @file numbers.cpp
 * @brief Numerical utilities Source
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#include <cmath>

#include "utils/numbers.h"

namespace dominiqs {

double getMachineEps()
{
	double one = 1.0;
	double eps = 1.0;
	double res;
	do
	{
		res = eps;
		eps /= 2.0;
	}
	while( (one + eps) > one );
	return res;
}

bool float2frac(double r, double mindelta, double maxdelta, int maxden, int& n, int& d)
{
	double epsilon = std::min(-mindelta, maxdelta) / 2.0;
	double b = r;
	double a = std::floor(b + epsilon);
	double n0 = a;
	double n1 = 1.0;
	double d0 = 1.0;
	double d1 = 0.0;
	double delta = r - n0 / d0;
	double nx;
	double dx;

	while( delta < mindelta || delta > maxdelta )
	{
		b = 1.0 / (b - a);
		a = std::floor(b + epsilon);
		nx = n0;
		dx = d0;

		n0 = a * n0 + n1;
		d0 = a * d0 + d1;

		n1 = nx;
		d1 = dx;

		if( d0 > maxden ) return false;

		delta = r - n0/d0;
	}

	n = (int)n0;
	d = (int)d0;
	return true;
}

static const double simpledens[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 25.0, -1.0};

bool float2simplefrac(double r, double mindelta, double maxdelta, int maxden, int& n, int& d)
{
	// try the simple denominators
	for(int i = 0; simpledens[i] > 0.0; i++)
	{
		double num;
		double delta;
		// try powers of 10 of the current simple denominator
		double den = simpledens[i];
		while( den <= maxden )
		{
			num = std::floor(r * den);
			delta = r - num / den;
			if( mindelta <= delta && delta <= maxdelta )
			{
				n = (int)num;
				d = (int)den;
				// simplify
				int common = gcd(n,d);
				n /= common;
				d /= common;
				return true;
			}
			den *= 10.0;
		}
	}
	return false;
}

bool findFrac(double lb, double ub, int maxden, int& n, int& d)
{
	double center = 0.5*(lb+ub);
	double delta = 0.5*(ub-lb);
	return float2frac(center, -delta, +delta, maxden, n, d);
}

} // namespace dominiqs
