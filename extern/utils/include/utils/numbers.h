/**
 * @file numbers.h
 * @brief Numerical utilities Header
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef NUMBERS_H
#define NUMBERS_H

#include <algorithm>
#include <iostream>

#include "asserter.h"

namespace dominiqs {

/**
 * returns the smallest number eps > 0, for which 1.0 + eps > 1.0 on the current machine
 */

double getMachineEps();

/**
 * calculate the greatest common divisor of two integral types n and m
 */

template <typename IntType> IntType gcd(IntType n, IntType m)
{
	IntType zero(0);
	if (n < zero) n = -n;
	if (m < zero) m = -m;
		// check for zeros
	if (m == zero) return n;
	if (n == zero) return m;
	DOMINIQS_ASSERT( n > zero );
	DOMINIQS_ASSERT( m > zero );
	// binary Stein algorithm
	IntType const one(1);
	// the greatest divisors of m, n that are powers of 2
	IntType mOdds(zero);
	IntType nOdds(zero);
	IntType const mask(one); //< binary mask to check for % 2
	// n,m stripped of trailing zeros
	for( ; !(m & mask) ; m >>= 1, mOdds++ );
	for( ; !(n & mask) ; n >>= 1, nOdds++ );
	for( ; n != m ; )
	{
		// n == m is gcd
		if( n < m ) for( m -= n, m >>= 1; !(m & mask) ; m >>= 1 ); //< strips trailing zeros from m
		else for( n -= m, n >>= 1; !(n & mask) ; n >>= 1 );
	}
	return m << (std::min)( mOdds, nOdds );
}

/**
 * calculates the least common multiple of two integral types n and m
 */

template <typename IntType> IntType lcm(IntType n, IntType m)
{
	return n / gcd<IntType>(n, m) * m;
}

/**
 * Compute the best approximating fraction n/d to the value r
 * such that  mindelta <= r - n/d <= maxdelta and d <= maxden
 * @return true in case of success, false otherwise
 */

bool float2frac(double r, double mindelta, double maxdelta, int maxden, int& n, int& d);

/**
 * Compute a simple fraction n/d approximating the value r
 * such that  mindelta <= r - n/d <= maxdelta and d <= maxden
 * @return true in case of success, false otherwise
 */

bool float2simplefrac(double r, double mindelta, double maxdelta, int maxden, int& n, int& d);

/**
 * Find a fraction n/d in the (very small) range [lb, ub]
 * auch that d <= maxden
 * @return true if successful, false otherwise
 */

bool findFrac(double lb, double ub, int maxden, int& n, int& d);

/**
 * Calculate the binomial coefficient (n choose m)
 * @return the binomial coefficient
 */

template<typename IntType>
double choose(IntType m, IntType n)
{
	if (n > m) return 0;
	if (n == m) return 1;
	double r = 1;
	for (IntType d = 1; d <= n; d++)
	{
		r *= m--;
		r /= d;
	}
	return r;
}

} // namespace dominiqs

#endif /* NUMBERS_H */
