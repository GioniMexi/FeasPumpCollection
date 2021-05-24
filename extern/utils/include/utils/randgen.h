/**
 * @file randgen.h
 * @brief Random Number Generator
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2007
 */

#ifndef RANDGEN_H
#define RANDGEN_H

#include <cmath>
#include <ctime>
#include <random>

#include "asserter.h"

namespace dominiqs {

static const uint64_t DEFAULT_SEED = 20070512;

/**
 * Generate a seed for the generator based on user input
 * If user input is not null, use it as seed
 * Otherwise use the default seed.
 */

inline uint64_t generateSeed(uint64_t seed)
{
	if (seed) return seed;
	return DEFAULT_SEED;
}

/**
 * Basic Random numbers generator class
 * Generates a floating number in the range [0,1) with uniform distribution
 */

class RandGen
{
public:
	/**
	 * default constructor
	 * @param _seed (uses generateSeed!)
	 * @param min minimum value of the range
	 * @param max maximum value of the range
	 */
	RandGen(uint64_t _seed = 0, double min = 0.0, double max = 1.0) : engine{(unsigned int)generateSeed(_seed)}, rnd{min, max} {}
	/**
	 * Seed the generator
	 * @param _seed: if 0 use the default seed, otherwise the supplied seed
 	 */
	void setSeed(uint64_t _seed) { engine.seed(generateSeed(_seed)); }
	/**
	 * set uniform distribution range
	 * @param _min lower range value
	 * @param _max upper range value
	 */
	void setRange(double min, double max)
	{
		std::uniform_real_distribution<double>::param_type pt(min, max);
		rnd.param(pt);
	}
	/** @return a random floating point value in the specified range */
	double getFloat() { return rnd(engine); }
	/** @return a random integer value in the specified range */
	long getInteger() { return lrint(getFloat()); }
	/** call the generator a few times to discard first values */
	void warmUp() { engine.discard(WARMUP_TRIES); }
protected:
	std::default_random_engine engine;
	std::uniform_real_distribution<double> rnd;
	static const int WARMUP_TRIES = 1000;
};


/**
 * Random numbers generator class
 * Generates integers in the range [0,N) in a uniform distribution
 * This interface is needed for STL algorithms
 */

class STLRandGen : public RandGen
{
public:
	/**
	 * default constructor
	 * @param _seed (uses generateSeed!)
	 */
	STLRandGen(uint64_t _seed = 0) : RandGen(_seed) {}
	long operator()(long N)
	{
		long res = lrint(getFloat() * (N - 1));
		DOMINIQS_ASSERT( res >= 0 );
		DOMINIQS_ASSERT( res < N );
		return res;
	}
};

/**
 * Random alphanumeric character generator
 */
namespace details
{
	static const char* genChars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
}

class RandCharGen
{
public:
	/**
	 * default constructor
	 * @param _seed (uses generateSeed!)
	 */
	RandCharGen(uint64_t _seed = 0) : engine{(unsigned int)generateSeed(_seed)}, rnd{0, RANDGEN_NUM_CHARS - 1} {}
	/**
	 * Seed the generator
	 * @param _seed: if 0 use the default seed, otherwise the supplied seed
 	 */
	void setSeed(uint64_t _seed) { engine.seed(generateSeed(_seed)); }
	/** @return random character */
	char operator()()
	{
		return details::genChars[rnd(engine)];
	}
	/** call the generator a few times to discard first values */
	void warmUp() { engine.discard(WARMUP_TRIES); }
private:
	std::default_random_engine engine;
	std::uniform_int_distribution<int> rnd;
	static const int WARMUP_TRIES = 1000;
	static const int RANDGEN_NUM_CHARS = 62;
};

} // namespace dominiqs

#endif /* RANDGEN_H */
