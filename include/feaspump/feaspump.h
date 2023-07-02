/**
 * @file feaspump.h
 * @brief Feasibility Pump algorithm Header
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * @author Gioni Mexi <gionimexi at gmail dot com>
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
	int stage1IterLimit;
	int stage2IterLimit;
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

	// new parameters
 	int itr1NoImpr; 	// max iterations without 10% improvent in stage 1
	int itr2NoImpr;		// max iterations without 10% improvent in stage 2
	bool expObj;
	bool logisObj;
	bool analcenterFP;

	bool newScaleC;	// new scaling for c in the obj fp
	bool newScaleDelta;	// new scaling for Delta in the obj fp

	int numIntegersObj;	// number of integer points used in obj
	bool aggregateInts; // accumulate integer points before creating objective (not good)
	int numFracsObj;    // number of fractional points used (not good)
	double fracScaleFactor; // scaling factor for multiple fractional points in objective
	double intScaleFactor; // scaling factor for multiple integer points in objective
	bool lessViolatedIntegers; // pick reference points with smallest violation of constraints
	bool lessDistanceIntegers; // pick reference points with smallest distance to projection LP solution
	bool bestObjIntegers; // pick reference points with best objective
	bool harmonicWeights; // harmonic weights for multiple integers
	bool exponDecayWeights; // exponential decay weights for multiple integers -> with factor intScaleFactor, e.g. 0.5

	bool normalMIPStage3;  // just solves the normal mip in stage 3 (for comparison purposes)
	bool newStage3;  // use more integers in stage 3
	bool rensStage3; // apply rens instead of stage 3
	bool rensClosestDistStage3; // apply rens on relaxation of closest integer
	bool multirensStage3; // rens with multiple reference points
	bool stage3lessViolatedIntegers; // pick reference points with smallest violation of constraints in stage 3
	bool stage3bestObjIntegers; // pick reference points with best objective in stage 3
	bool stage3harmonicWeights; // harmonic weights for multiple integers in stage 3
	int stage3IntegersObj;	// number of integer points used in the stage 3 objective
	double stage3ScaleFactor; // scaling factor for multiple integer points in stage 3 objective
	double stage3Time; 

	double pdlpTol; // tolerance for pdlp
	double pdlpTolDecreaseFactor; // decrease factor for pdlp tolerance
	bool pdlpWarmStart; // true if pdlp should use warm start

	// LP options
	char firstOptMethod;
	char reOptMethod;
	// FP data
	MIPModelPtr model;
	double objOffset;
	SolutionTransformerPtr frac2int; /**< rounder */
	std::vector<double> frac_x; /**< fractional x^* */
	std::vector<double> ac_x; /**< analytic center */
	int primalFeas; /**< is current fractional x^* primal feasible? */
	std::vector<double> integer_x; /**< integer x^~ */
	typedef std::pair<double, std::vector<double>> AlphaVector;
	std::list<AlphaVector> lastIntegerX; /**< integer x cache */

	typedef std::pair<int, std::vector<double>> NumberVector;
	std::list<NumberVector> lastFracX; /**< fractional x cache */

	typedef std::pair<double, std::vector<double>> DistVector;
	std::list<DistVector> multipleIntegerX; /**< integer reference points */
	std::list<DistVector> closestIntegerXs1; /**< closest integer x cache stage 1*/
	std::list<DistVector> closestIntegerXs2; /**< closest integer x cache stage 2*/


	RandGen rnd;
	std::vector<double> closestPoint; /**< point closest to feasibility */
	std::vector<double> closestFrac; /**< projection of closest point to feasibility */

	double closestDist;
	// problem data
	bool isPureInteger; /**< is it a mixed linear program? */
	bool isBinary; /**< are there general integers? */
	std::vector<double> obj; /**< original objective function */
	double objNorm;
	std::vector<double> lb; /**< lower bounds */
	std::vector<double> ub; /**< upper bounds */
	std::vector<double> newlb; /**< lower bounds */
	std::vector<double> newub; /**< upper bounds */


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
	double acTime;
	int rootLpIter;
	// helpers
	void solveInitialLP();
	void perturbe(std::vector<double>& x, bool ignoreGeneralIntegers);
	void restart(std::vector<double>& x, bool ignoreGeneralIntegers);
	bool pumpLoop(double& runningAlpha, int stage, double& dualBound, double& timeModel);
	bool stage3();
	void foundIncumbent(const std::vector<double>& x, double objval);
	bool isInCache(double a, const std::vector<double>& x, bool ignoreGeneralIntegers);
	void infeasibleSupport(const std::vector<double>& x, std::set<int>& supp, bool ignoreGeneralIntegers);


	// added function
	void aggregateFracs(std::vector<double>& aggr_frac_x, std::vector<double>& scaleVector);
	void computeAC( std::vector<double>& x);
	void integerFromAC(std::vector<double>& x, double& bestgamma, double step);
	bool isComponentSame(int idx);
	bool isComponentSameS3(int idx);

};

} // namespace dominiqs

#endif /* FEASPUMP_H */
