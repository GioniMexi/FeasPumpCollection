/**
 * @file ranking.h
 * @brief Variable Rankers Header
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2008
 */

#ifndef RANKING_H
#define RANKING_H

#include <vector>
#include <set>
#include <memory>

#include <utils/randgen.h>
#include <utils/singleton.h>
#include <utils/factory.h>
#include <propagator/domain.h>

// forward declarations
class Propagator;

/**
 * @brief Ranker Base Class
 */

class Ranker
{
public:
	virtual ~Ranker() {}
	virtual void readConfig() {}
	virtual void init(DomainPtr d, bool ignoreGeneralInt = true);
	virtual void ignoreGeneralIntegers(bool flag);
	virtual void setCurrentState(const std::vector<double>& x) = 0;
	virtual int next() = 0;
protected:
	DomainPtr domain;
	std::vector<int> binaries;
	std::vector<int> gintegers;
	std::vector<int> integers;
	int nCalled = 0;
};

typedef std::shared_ptr<Ranker> RankerPtr;

typedef dominiqs::SingletonHolder< dominiqs::Factory<Ranker, std::string> > RankerFactory;

/**
 * @brief Rank variables from left to right
 */

class LeftToRightRanker : public Ranker
{
public:
	void setCurrentState(const std::vector<double>& x);
	int next();
protected:
	// data
	unsigned int nextItr;
};

/**
 * @brief Rank variables in order of increasing fractionality
 * always preferring binary variables to general integer ones
 */

class FractionalityRanker : public Ranker
{
public:
	void readConfig();
	void ignoreGeneralIntegers(bool flag);
	void setCurrentState(const std::vector<double>& x);
	int next();
protected:
	// options
	bool reverse;
	double rankNoise;
	int noiseAfter;
	dominiqs::STLRandGen rnd;
	// data
	typedef std::pair<double, int> Score;
	std::vector<Score> scores;
	std::vector<int> perm;
	unsigned int nextItr;
};

/**
 * @brief Rank variables randomly
 */

class RandomRanker : public Ranker
{
public:
	void readConfig();
	void ignoreGeneralIntegers(bool flag);
	void setCurrentState(const std::vector<double>& x);
	int next();
protected:
	// data
	dominiqs::STLRandGen rnd;
	std::vector<int> perm;
	unsigned int nextItr;
};

#endif /* RANKING_H */
