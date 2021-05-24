/**
 * @file metaheur_basics.h
 * @brief Basic buildin blocks for metaheuristics
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2018 Domenico Salvagnin
 */

#ifndef METAHEUR_BASICS_H
#define METAHEUR_BASICS_H

#include <stdexcept>
#include <limits>
#include <random>
#include "floats.h"
#include "timer.h"

namespace dominiqs {

class EmptyUserData {};

// Interface for metaheuristics:
// users should overload the following functions to provide support for their types
// In all functions below, UserData is an opaque template parameter that is needed
// to support customized behaviour for the same Solution type.
//
// /**
//  * Get the objective value (of type ObjType) of a given solution of type Solution
//  */
// ObjType getSolutionObjValue(const Solution& sol, const UserData& ud);
//
// /**
//  * Generate a random solution of type Solution, using a random engine of type RandomEngine
//  */
// void randomGenerateSolution(Solution& sol, RandomEngine& engine, const UserData& ud);
//
// /**
//  * Generate a given solution of type Solution, using a random engine of type RandomEngine, and integer parameter K
//  */
// void perturbeSolution(Solution& sol, RandomEngine& engine, int K, const UserData& ud);
//
// /**
//  * Explore the neighborhood of a solution of type Solution, using an observer of type Observer
//  */
// bool exploreNeighborhood(Solution& sol, Observer& rec, const UserData& ud);


/**
 * Interface for observers of solutions of type Solution and objective type ObjType
 */
template<typename Solution, typename ObjType = double>
class IObserver
{
public:
	virtual ~IObserver() {}
	virtual void solution(const Solution& sol) {}
	virtual void iteration(int iter, ObjType currentObj, ObjType bestObj) {}
};


/**
 * Termination Criteria class for localsearch-based metaheuristics
 * It can deal with:
 * - maximum number of iterations
 * - maximum number of iterations without improvement
 * - time limit
 * - abort flag (to support Ctrl-C)
 */
class TerminationCriteria
{
public:
	TerminationCriteria(StopWatch& sw = gStopWatch()) : stopwatch(sw) {}
	// params
	int iterLimit = std::numeric_limits<int>::max();
	int noImproveLimit = std::numeric_limits<int>::max();
	double timeLimit = std::numeric_limits<double>::max(); //< time limit in seconds
	int* terminate = nullptr;
	// clock
	StopWatch& stopwatch;
	// exec
	bool operator()(int currentItr, int noImproveCnt)
	{
		if (terminate && *terminate) return true;
		if (currentItr > iterLimit) return true;
		if (noImproveCnt > noImproveLimit) return true;
		if (stopwatch.getElapsed() > timeLimit) return true;
		return false;
	}
};


/* Interface for acceptance criteria for ILS of solutions of type Solution and objective type ObjType */
template<typename Solution, typename ObjType = double, typename UserData = EmptyUserData>
class IAcceptanceCriterion
{
public:
	IAcceptanceCriterion(const UserData& _ud = UserData()) : ud(_ud) {}
	virtual ~IAcceptanceCriterion() {}
	virtual void operator()(int iteration, Solution& oldSol, Solution& newSol) {}
protected:
	const UserData& ud;
};

template<typename Solution, typename ObjType = double, typename UserData = EmptyUserData>
class HillClimbingCriterion : public IAcceptanceCriterion<Solution,ObjType,UserData>
{
public:
	void operator()(int iteration, Solution& oldSol, Solution& newSol)
	{
		if (dominiqs::lessThan(getSolutionObjValue(newSol, this->ud), getSolutionObjValue(oldSol, this->ud))) oldSol = newSol;
	}
};

template<typename Solution, typename ObjType = double, typename UserData = EmptyUserData>
class RandomWalkCriterion : public IAcceptanceCriterion<Solution,ObjType,UserData>
{
public:
	void operator()(int iteration, Solution& oldSol, Solution& newSol)
	{
		oldSol = newSol;
	}
};

template<typename Solution, typename RandomEngine = std::default_random_engine, typename ObjType = double, typename UserData = EmptyUserData>
class AnnealingCriterion : public IAcceptanceCriterion<Solution,ObjType,UserData>
{
public:
	AnnealingCriterion(double T, RandomEngine& e, const UserData& _ud = UserData())
		: IAcceptanceCriterion<Solution,ObjType,UserData>(_ud), temperature(T), engine(e), rnd(0.0, 1.0) {}
	void operator()(int iteration, Solution& oldSol, Solution& newSol)
	{
		ObjType energyOld = getSolutionObjValue(oldSol, this->ud);
		ObjType energyNew = getSolutionObjValue(newSol, this->ud);
		if (iteration && ((iteration % interval) == 0)) temperature *= ratio;
		if (dominiqs::lessThan(energyNew, energyOld)) oldSol = newSol;
		else
		{
			double prob = std::exp(ObjType(energyOld - energyNew) / temperature);
			if (rnd(engine) <= prob) oldSol = newSol;
		}
	}
private:
	double temperature;
	double ratio = 0.99;
	int interval = 100;
	RandomEngine& engine;
	std::uniform_real_distribution<double> rnd;
};

} // namespace dominiqs

#endif /* end of include guard: METAHEUR_BASICS_H */
