/**
 * @file local_search.h
 * @brief Local search based metaheuristics: LS, ILS and RRLS
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * @author Gioni Mexi <gionimexi at gmail dot com>
 * Copyright 2018 Domenico Salvagnin
 */

#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include "metaheur_basics.h"

namespace dominiqs {

/**
 * Local search for solutions of type Solution and observers of type Observer
 */
template<typename Solution, typename Observer = IObserver<Solution,double>, typename UserData = EmptyUserData>
void localSearch(Solution& sol, Observer& observer, const UserData& ud = UserData())
{
	observer.solution(sol);
	while(exploreNeighborhood(sol, observer, ud)) {}
}

/**
 * ILS for solutions of type Solution and observers of type Observer,
 * using a random engine of type RandomEngine, acceptance criterion of type AcceptanceCriterion
 * and objective value of type ObjValue
 */
template<typename Solution,
	typename AcceptanceCriterion = IAcceptanceCriterion<Solution>,
	typename RandomeEngine = std::default_random_engine,
	typename ObjValue = double,
	typename Observer = IObserver<Solution,ObjValue>,
	typename UserData = EmptyUserData>
void iteratedLocalSearch(Solution& sol,
	AcceptanceCriterion& accept,
	RandomeEngine& rng,
	Observer& observer,
	TerminationCriteria& toStop,
	int k_min, int k_max = -1,
	const UserData& ud = UserData())
{
	if (k_max == -1) k_max = k_min;
	int k = k_min;

	int noImproveCount = 0;
	int iter = 0;

	// generate a random starting point and perform local search
	Solution currSol(sol);
	randomGenerateSolution(currSol, rng, ud);
	localSearch(currSol, observer, ud);
	// keep best
	ObjValue bestObj = getSolutionObjValue(currSol, ud);
	sol = currSol;

	// main loop
	while(true)
	{
		// iterated local search main step
		if (toStop(iter, noImproveCount)) break;
		Solution candidate(currSol);
		perturbeSolution(candidate, rng, k, ud);
		localSearch(candidate, observer, ud);
		accept(iter, currSol, candidate);

		// keep best
		ObjValue currentObj = getSolutionObjValue(currSol, ud);
		if (currentObj < bestObj - 0.5)
		{
			bestObj = currentObj;
			sol = currSol;
			noImproveCount = 0;
			k = k_min;
		}
		else
		{
			noImproveCount++;
			k++;
			if (k > k_max) k = k_min;
		}

		// optional restart
		if (noImproveCount > toStop.noImproveLimit)
		{
			randomGenerateSolution(currSol, rng, ud);
			noImproveCount = 0;
		}

		observer.iteration(iter, currentObj, bestObj);
		iter++;
	}
}

/**
 * RRLS for solutions of type Solution and solution recorders of type SolutionRecorder,
 * using a random engine of type RandomEngine and objective value of type ObjValue
 */
template<typename Solution,
	typename RandomeEngine = std::default_random_engine,
	typename ObjValue = double,
	typename Observer = IObserver<Solution,ObjValue>,
	typename UserData = EmptyUserData>
void randomRestartLocalSearch(Solution& sol,
	RandomeEngine& rng,
	Observer& observer,
	TerminationCriteria& toStop,
	const UserData& ud = UserData())
{
	ObjValue bestObj = std::numeric_limits<ObjValue>::max();
	int iter = 0;
	int noImproveCount = 0;

	while(true)
	{
		if (toStop(iter, noImproveCount)) break;

		Solution currSol(sol);
		// generate a random starting point
		randomGenerateSolution(currSol, rng, ud);
		// local search
		localSearch(currSol, observer, ud);
		ObjValue currentObj = getSolutionObjValue(currSol, ud);
		if (currentObj < bestObj - 0.5)
		{
			bestObj = currentObj;
			sol = currSol;
			noImproveCount = 0;
		}
		else noImproveCount++;

		observer.iteration(iter, currentObj, bestObj);
		iter++;
	}
}

} // namespace dominiqs

#endif /* end of include guard: LOCAL_SEARCH_H */
