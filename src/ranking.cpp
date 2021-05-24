/**
 * @file ranking.cpp
 * @brief Variable Rankers Source
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 */

#include <algorithm>
#include <iostream>

#include <utils/sorting.h>
#include <utils/fileconfig.h>
#include <utils/consolelog.h>

#include "feaspump/ranking.h"

using namespace dominiqs;


// macro type savers
#define READ_FROM_CONFIG( what, defValue ) what = gConfig().get("fp."#what, defValue)
#define LOG_ITEM(name, value) consoleLog("{} = {}", name, value)
#define LOG_CONFIG( what ) LOG_ITEM("fp."#what, what)


// Base Ranker Class

void Ranker::init(DomainPtr d, bool ignoreGeneralInt)
{
	DOMINIQS_ASSERT( d );
	nCalled = 0;
	domain = d;
	binaries.clear();
	gintegers.clear();
	unsigned int n = domain->size();
	for (unsigned int i = 0; i < n; i++)
	{
		if (!domain->isVarFixed(i))
		{
			if (domain->varType(i) == 'B') binaries.push_back(i);
			if (domain->varType(i) == 'I') gintegers.push_back(i);
		}
	}
	ignoreGeneralIntegers(ignoreGeneralInt);
}

void Ranker::ignoreGeneralIntegers(bool flag)
{
	if (flag) integers = binaries;
	else
	{
		// append general integers to binaries
		integers.resize(binaries.size() + gintegers.size());
		std::vector<int>::iterator itr = copy(binaries.begin(), binaries.end(), integers.begin());
		std::vector<int>::iterator end = copy(gintegers.begin(), gintegers.end(), itr);
		DOMINIQS_ASSERT( end == integers.end() );
	}
}

// LeftToRight

void LeftToRightRanker::setCurrentState(const std::vector<double>& x)
{
	nextItr = 0;
	nCalled++;
}

int LeftToRightRanker::next()
{
	while (nextItr < integers.size() && domain->isVarFixed(integers[nextItr])) nextItr++;
	if (nextItr < integers.size()) return integers[nextItr];
	return -1;
}

// FractionalityRanker

static const bool FRAC_RANKER_REVERSE_DEF = false;
static const double FRAC_RANKER_RANK_NOISE_DEF = 0.1;
static const int FRAC_RANKER_NOISE_AFTER_DEF = 10;

void FractionalityRanker::readConfig()
{
	READ_FROM_CONFIG( reverse, FRAC_RANKER_REVERSE_DEF );
	READ_FROM_CONFIG( rankNoise, FRAC_RANKER_RANK_NOISE_DEF );
	READ_FROM_CONFIG( noiseAfter, FRAC_RANKER_NOISE_AFTER_DEF );
	uint64_t seed = gConfig().get<uint64_t>("seed", 0);
	rnd.setSeed(seed);
	rnd.warmUp();
	consoleInfo("[config ranker]");
	LOG_CONFIG( reverse );
	LOG_CONFIG( rankNoise );
	LOG_CONFIG( noiseAfter );
}

void FractionalityRanker::ignoreGeneralIntegers(bool flag)
{
	Ranker::ignoreGeneralIntegers(flag);
	// reserve memory for scores and permutation
	scores.resize(integers.size());
	perm.resize(integers.size());
}

void FractionalityRanker::setCurrentState(const std::vector<double>& x)
{
	// calculate scores
	for (unsigned int i = 0; i < integers.size(); i++)
	{
		double s = integralityViolation(x[integers[i]]); // number between 0 and 0.5
		if (domain->varType(integers[i]) == 'B') s -= 10; // binaries score is always less than general integers
		scores[i] = Score(s, integers[i]);
		perm[i] = i;
	}
	// sort
	if (reverse) permShellSort(&scores[0], &perm[0], scores.size(), std::greater<Score>(), true);
	else permShellSort(&scores[0], &perm[0], scores.size(), std::less<Score>(), true);
	// perturbe sorting
	if (isNotNull(rankNoise) && (nCalled > noiseAfter))
	{
		int n = perm.size();
		int nSwaps = int(rankNoise * n);
		for (int i = 0; i < nSwaps; i++)
		{
			int fromIdx = rnd(n);
			int toIdx = rnd(n);
			std::swap(perm[fromIdx], perm[toIdx]);
		}
	}
	nextItr = 0;
	nCalled++;
}

int FractionalityRanker::next()
{
	while (nextItr < integers.size() && domain->isVarFixed(integers[perm[nextItr]])) nextItr++;
	if (nextItr < integers.size()) return integers[perm[nextItr]];
	return -1;
}

// RandomRanker

void RandomRanker::readConfig()
{
	uint64_t seed = gConfig().get<uint64_t>("seed", 0);
	rnd.setSeed(seed);
	rnd.warmUp();
}

void RandomRanker::ignoreGeneralIntegers(bool flag)
{
	Ranker::ignoreGeneralIntegers(flag);
	// reserve memory for permutation
	perm.resize(integers.size());
}

void RandomRanker::setCurrentState(const std::vector<double>& x)
{
	std::iota(perm.begin(), perm.end(), 0);
	std::random_shuffle(perm.begin(), perm.end(), rnd);
	nextItr = 0;
	nCalled++;
}

int RandomRanker::next()
{
	while (nextItr < integers.size() && domain->isVarFixed(integers[perm[nextItr]])) nextItr++;
	if (nextItr < integers.size()) return integers[perm[nextItr]];
	return -1;
}

// auto registration

class RANK_FACTORY_RECORDER
{
public:
	RANK_FACTORY_RECORDER()
	{
		std::cout << "Registering Var Rankers...";
		RankerFactory::getInstance().registerClass<LeftToRightRanker>("LR");
		RankerFactory::getInstance().registerClass<FractionalityRanker>("FRAC");
		RankerFactory::getInstance().registerClass<RandomRanker>("RND");
		std::cout << "done" << std::endl;
	}
};

RANK_FACTORY_RECORDER my_rank_factory_recorder;
