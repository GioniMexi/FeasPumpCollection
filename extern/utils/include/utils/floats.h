/**
 * @file floats.h
 * @brief Floating Point Utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef FLOATS_H
#define FLOATS_H

#include <cmath>
#include <algorithm>

namespace dominiqs {

const double defaultEPS = 1e-6;

inline bool equal(double x, double y, double eps = defaultEPS)
{
	return (fabs(x - y) <= eps);
}

inline bool relEqual(double x, double y, double eps = defaultEPS)
{
	if (x == y) return true;
	double tol = std::max(fabs(x), fabs(y));
	return (fabs(x - y) <= eps * (1 + tol));
}

inline bool different(double x, double y, double eps = defaultEPS)
{
	return fabs(x -y) > eps;
}

inline bool lessThan(double x, double y, double eps = defaultEPS)
{
	return (x - y) < -eps;
}

inline bool greaterThan(double x, double y, double eps = defaultEPS)
{
	return (x - y) > eps;
}

inline bool lessEqualThan(double x, double y, double eps = defaultEPS)
{
	return (x - y) <= eps;
}

inline bool greaterEqualThan(double x, double y, double eps = defaultEPS)
{
	return (x - y) >= -eps;
}

inline bool isNegative(double x, double eps = defaultEPS)
{
	return x < - eps;
}

inline bool isPositive(double x, double eps = defaultEPS)
{
	return x > eps;
}

inline bool isNull(double x, double eps = defaultEPS)
{
	return (fabs(x) <= eps);
}

inline bool isNotNull(double x, double eps = defaultEPS)
{
	return (fabs(x) > eps);
}

inline double floorEps(double x, double eps = defaultEPS)
{
	double tmp = x + eps;
	double r = static_cast<double>(static_cast<int>(tmp));
	return (r - (r > tmp));
}

inline double ceilEps(double x, double eps = defaultEPS)
{
	double tmp = x - eps;
	double r = static_cast<double>(static_cast<int>(tmp));
	return (r + (r < tmp));
}

inline double fractionalPart(double x, double eps = defaultEPS)
{
	double ret = x - floorEps(x, eps);
	return std::max(ret, 0.0);
}

inline double integralityViolation(double x, double eps = defaultEPS)
{
	double f = fractionalPart(x, eps);
	return std::min(f, (1.0 - f));
}

inline bool isInteger(double x, double eps = defaultEPS)
{
	return isNull(fractionalPart(x, eps), eps);
}

template<class int_type>
inline int_type double2int(double x)
{
	if (isNull(x)) return (int_type)0;
	return (x > 0) ? (int_type)(x+0.5) : (int_type)(x-0.5);
}

} // namespace dominiqs

#endif /* FLOATS_H */
