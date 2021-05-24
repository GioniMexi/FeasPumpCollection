/**
 * @file transformers.h
 * @brief Solution Transformers (a.k.a. rounders) Header
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2008
 */

#ifndef TRANSFORMERS_H
#define TRANSFORMERS_H

#include <utils/randgen.h>

#include <propagator/domain.h>
#include <propagator/prop_engine.h>

#include "fp_interface.h"
#include "ranking.h"

/**
 * Rounding helpers
 */

inline void doRound(const double& in, double& out, const double& thr) { out = floor(in + thr); }

inline double getRoundingThreshold(bool isRandom, dominiqs::RandGen& gen)
{
	if (isRandom)
	{
		//double t = gen.getFloat();
		//if (t <= 0.5) return 2 * t * (1 - t);
		//else return 2 * t * (t - 1) + 1;
		return gen.getFloat();
	}
	else return 0.5;
}

/**
 * Standard frac -> int transformer: performs a simple rounding rounding
 */

class SimpleRounding : public dominiqs::SolutionTransformer
{
public:
	SimpleRounding();
	void readConfig();
	void init(MIPModelPtr model, bool ignoreGeneralInt = true);
	void ignoreGeneralIntegers(bool flag);
	void apply(const std::vector<double>& in, std::vector<double>& out);
protected:
	std::vector<int> binaries;
	std::vector<int> gintegers;
	std::vector<int> integers;
	dominiqs::RandGen roundGen;
	bool randomizedRounding;
	bool logDetails;
};

/**
 * Ranked rounding + constraint propagation
 */

class PropagatorRounding : public SimpleRounding
{
public:
	PropagatorRounding();
	~PropagatorRounding() { clear(); }
	void readConfig();
	void init(MIPModelPtr model, bool ignoreGeneralInt = true);
	void ignoreGeneralIntegers(bool flag);
	void apply(const std::vector<double>& in, std::vector<double>& out);
	void clear();
protected:
	// data
	DomainPtr domain;
	StatePtr state;
	PropagationEngine prop;
	std::map<int, PropagatorFactoryPtr> factories;
	RankerPtr ranker;
	bool filterConstraints;
};

#endif /* TRANSFORMERS_H */
