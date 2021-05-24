/**
 * @file feaspump.h
 * @brief Feasibility Pump algorithm Header
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2008
 */

#ifndef FEASPUMP_H
#define FEASPUMP_H

#include <list>
#include <set>

#include <utils/randgen.h>
#include <utils/it_display.h>
#include <utils/timer.h>

#include "fp_interface.h"

namespace dominiqs {

/**
 * @brief Basic Feasibility Pump Scheme
 * Accepts custom rounders
 */

class FeasibilityPump
{
public:
	FeasibilityPump();
	// config
	void readConfig();
	/** init algorithm
	 * @param env: cplex environment
	 * @param lp: problem object (this is modified by the algorithm: you may want to pass a copy!)
	 * @param ctype: if the problem object is an LP, you can provide variable type info with this vector
	 */
	void init(MIPModelPtr model, const std::vector<char>& ctype = std::vector<char>());
	/** pump
	 * @param xStart: starting fractional solution (will solve LP if empty)
	 * @param pFeas: primal feasiblity status of supplied vector
	 */
	bool pump(const std::vector<double>& xStart = std::vector<double>(), bool pFeas = false);
	// get solution info
	bool foundSolution() const;
	void getSolution(std::vector<double>& x) const;
	double getSolutionValue(const std::vector<double>& x) const;
	int getIterations() const;
	// reset
	void reset();
private:
	// FP options
	double timeLimit;
	double timeMult;
	double lpIterMult;
	int stageIterLimit;
	int iterLimit;
	int avgFlips;
	double integralityEps;
	uint64_t seed;
	double alpha;
	double alphaFactor;
	double alphaDist;
	bool doStage3;
	bool walksatPerturbe;
	bool randomizeLP;
	bool penaltyObj;
	// LP options
	char firstOptMethod;
	char reOptMethod;
	// FP data
	MIPModelPtr model;
	double objOffset;
	SolutionTransformerPtr frac2int; /**< rounder */
	std::vector<double> frac_x; /**< fractional x^* */
	int primalFeas; /**< is current fractional x^* primal feasible? */
	std::vector<double> integer_x; /**< integer x^~ */
	typedef std::pair<double, std::vector<double>> AlphaVector;
	std::list<AlphaVector> lastIntegerX; /**< integer x cache */
	RandGen rnd;
	std::vector<double> closestPoint; /**< point closest to feasibility */
	double closestDist;
	// problem data
	bool isPureInteger; /**< is it a mixed linear program? */
	bool isBinary; /**< are there general integers? */
	std::vector<double> obj; /**< original objective function */
	double objNorm;
	std::vector<double> lb; /**< lower bounds */
	std::vector<double> ub; /**< upper bounds */
	std::vector<char> xType; /**< column types */
	std::vector<bool> fixed; /**< fixed status in the original formulation (usually due to presolve) */
	std::vector<int> binaries; /**< list of binary vars indexes */
	std::vector<int> gintegers; /**< list of general integer vars indexes */
	std::vector<int> integers; /**< list of non continuous vars indexes (binaries + gintegers) */
	std::vector<ConstraintPtr> rows; /**< constraints of the model */
	IterationDisplay display;
	// solution
	bool hasIncumbent;
	std::vector<double> incumbent; /**< current incumbent */
	double primalBound;
	// stats
	int firstPerturbation;
	int pertCnt;
	int restartCnt;
	int walksatCnt;
	int nitr; /**< pumping iterations */
	int lastRestart;
	int flipsInRestart;
	int maxFlipsInRestart;
	StopWatch chrono;
	StopWatch lpWatch;
	StopWatch roundWatch;
	double rootTime;
	int rootLpIter;
	// helpers
	void solveInitialLP();
	void perturbe(std::vector<double>& x, bool ignoreGeneralIntegers);
	void restart(std::vector<double>& x, bool ignoreGeneralIntegers);
	bool pumpLoop(double& runningAlpha, int stage);
	bool stage3();
	void foundIncumbent(const std::vector<double>& x, double objval);
	bool isInCache(double a, const std::vector<double>& x, bool ignoreGeneralIntegers);
	void infeasibleSupport(const std::vector<double>& x, std::set<int>& supp, bool ignoreGeneralIntegers);
};

} // namespace dominiqs

#endif /* FEASPUMP_H */
