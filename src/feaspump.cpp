/**
 * @file feaspump.cpp
 * @brief Feasibility Pump algorithm Source
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 */

#include <map>
#include <signal.h>
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include <utils/asserter.h>
#include <utils/floats.h>
#include <utils/maths.h>
#include <utils/fileconfig.h>
#include <utils/consolelog.h>
#include <fmt/format.h>

#include "feaspump/feaspump.h"

using namespace dominiqs;

// macro type savers
#define READ_FROM_CONFIG( what, defValue ) what = gConfig().get("fp."#what, defValue)
#define LOG_ITEM(name, value) consoleLog("{} = {}", name, value)
#define LOG_CONFIG( what ) LOG_ITEM("fp."#what, what)


static bool isSolutionInteger(const std::vector<int>& integers, const std::vector<double>& x, double eps)
{
	for (int j: integers) if (!isInteger(x[j], eps)) return false;
	return true;
}

static int solutionNumFractional(const std::vector<int>& integers, const std::vector<double>& x, double eps)
{
	int tot = 0;
	for (int j: integers) if (!isInteger(x[j], eps)) tot++;
	return tot;
}

static double solutionsDistance(const std::vector<int>& integers, const std::vector<double>& x1, const std::vector<double>& x2)
{
	double tot = 0.0;
	for (int j: integers) tot += fabs(x1[j] - x2[j]);
	return tot;
}

static bool areSolutionsEqual(const std::vector<int>& integers, const std::vector<double>& x1, const std::vector<double>& x2, double eps)
{
	for (int j: integers) if (different(x1[j], x2[j], eps)) return false;
	return true;
}

static bool isSolutionFeasible(const std::vector<double>& x, const std::vector<ConstraintPtr>& rows)
{
	for (auto c: rows) {
		if (!c->satisfiedBy(&x[0]))  return false;
	}
	return true;
}


namespace dominiqs {

static const double DEF_TIME_LIMIT = 3600.0;
static const double DEF_TIME_MULT = 100.0;
static const double DEF_LPITER_MULT = -1.0;
static const int DEF_STAGE_ITER_LIMIT = 50;
static const int DEF_ITER_LIMIT = 100;
static const int DEF_AVG_FLIPS = 20;
static const double DEF_INTEGRALITY_EPS = 1e-6;
static const uint64_t DEF_SEED = 0;
static const double DEF_ALPHA = 0.0;
static const double DEF_ALPHA_FACTOR = 0.9;
static const double DEF_ALPHA_DIST = 0.005;
static const bool DEF_DO_STAGE_3 = false;
static const bool DEF_WALKSAT_PERTURBE = true;
static const bool DEF_RANDOMIZE_LP = false;
static const bool DEF_PENALTYOBJ = false;
static const char DEF_FIRST_OPT_METHOD = 'S';
static const char DEF_REOPT_METHOD = 'S';


static const double GEOM_FACTOR = 0.85;
static const double BIGM = 1e9;
static const double BIGBIGM = 1e15;
static const double INFBOUND = 1e20;


FeasibilityPump::FeasibilityPump() :
	timeLimit(DEF_TIME_LIMIT), timeMult(DEF_TIME_MULT), lpIterMult(DEF_LPITER_MULT),
	stageIterLimit(DEF_STAGE_ITER_LIMIT), iterLimit(DEF_ITER_LIMIT),
	avgFlips(DEF_AVG_FLIPS), integralityEps(DEF_INTEGRALITY_EPS), seed(DEF_SEED),
	alpha(DEF_ALPHA), alphaFactor(DEF_ALPHA_FACTOR), alphaDist(DEF_ALPHA_DIST),
	doStage3(DEF_DO_STAGE_3), walksatPerturbe(DEF_WALKSAT_PERTURBE),
	randomizeLP(DEF_RANDOMIZE_LP), penaltyObj(DEF_PENALTYOBJ),
	firstOptMethod(DEF_FIRST_OPT_METHOD), reOptMethod(DEF_REOPT_METHOD),
	objOffset(0.0), hasIncumbent(false)
{
}


void FeasibilityPump::readConfig()
{
	std::string frac2intName = gConfig().get("fp.frac2int", std::string("propround"));
	frac2int = SolutionTransformerPtr(TransformersFactory::getInstance().create(frac2intName));
	DOMINIQS_ASSERT( frac2int );
	// optimization methods
	std::string firstMethod = gConfig().get("fp.firstOptMethod", std::string("default"));
	if (firstMethod == "default") firstOptMethod = 'S';
	else if (firstMethod == "primal") firstOptMethod = 'P';
	else if (firstMethod == "dual") firstOptMethod = 'D';
	else if (firstMethod == "barrier") firstOptMethod = 'B';
	else throw std::runtime_error(std::string("Unknown optimization method: ") + firstMethod);
	std::string reMethod = gConfig().get("fp.reOptMethod", std::string("default"));
	if (reMethod == "default") reOptMethod = 'S';
	else if (reMethod == "primal") reOptMethod = 'P';
	else if (reMethod == "dual") reOptMethod = 'D';
	else if (reMethod == "barrier") reOptMethod = 'B';
	else throw std::runtime_error(std::string("Unknown optimization method: ") + reMethod);
	//other options
	READ_FROM_CONFIG( timeLimit, DEF_TIME_LIMIT );
	READ_FROM_CONFIG( timeMult, DEF_TIME_MULT );
	READ_FROM_CONFIG( lpIterMult, DEF_LPITER_MULT );
	READ_FROM_CONFIG( stageIterLimit, DEF_STAGE_ITER_LIMIT );
	READ_FROM_CONFIG( iterLimit, DEF_ITER_LIMIT );
	READ_FROM_CONFIG( avgFlips, DEF_AVG_FLIPS );
	READ_FROM_CONFIG( integralityEps, DEF_INTEGRALITY_EPS );
	seed = gConfig().get<uint64_t>("seed", DEF_SEED);
	READ_FROM_CONFIG( alpha, DEF_ALPHA );
	READ_FROM_CONFIG( alphaFactor, DEF_ALPHA_FACTOR );
	READ_FROM_CONFIG( alphaDist, DEF_ALPHA_DIST );
	READ_FROM_CONFIG( doStage3, DEF_DO_STAGE_3 );
	READ_FROM_CONFIG( walksatPerturbe, DEF_WALKSAT_PERTURBE );
	READ_FROM_CONFIG( randomizeLP, DEF_RANDOMIZE_LP );
	READ_FROM_CONFIG( penaltyObj, DEF_PENALTYOBJ );
	// display options
	display.headerInterval = gConfig().get("headerInterval", 10);
	display.iterationInterval = gConfig().get("iterationInterval", 1);
	// log config
	consoleInfo("[config fp]");
	LOG_ITEM("fp.frac2int", frac2intName);
	LOG_ITEM("fp.firstOptMethod", firstMethod);
	LOG_ITEM("fp.reOptMethod", reMethod);
	LOG_CONFIG( timeLimit );
	LOG_CONFIG( timeMult );
	LOG_CONFIG( lpIterMult );
	LOG_CONFIG( iterLimit );
	LOG_CONFIG( stageIterLimit );
	LOG_CONFIG( avgFlips );
	LOG_CONFIG( integralityEps );
	LOG_CONFIG( seed );
	LOG_CONFIG( alpha );
	LOG_CONFIG( alphaFactor );
	LOG_CONFIG( alphaDist );
	LOG_CONFIG( doStage3 );
	LOG_CONFIG( walksatPerturbe );
	LOG_CONFIG( randomizeLP );
	LOG_CONFIG( penaltyObj );
	rnd.setSeed(seed);
	rnd.warmUp();
	frac2int->readConfig();
}


bool FeasibilityPump::foundSolution() const
{
	return hasIncumbent;
}


void FeasibilityPump::getSolution(std::vector<double>& x) const
{
	DOMINIQS_ASSERT( hasIncumbent );
	x.resize(incumbent.size());
	copy(incumbent.begin(), incumbent.end(), x.begin());
}


double FeasibilityPump::getSolutionValue(const std::vector<double>& x) const
{
	return std::inner_product(x.begin(), x.end(), obj.begin(), 0.0) + objOffset;
}


int FeasibilityPump::getIterations() const
{
	return nitr;
}

void FeasibilityPump::reset()
{
	nitr = 0;
	firstPerturbation = 0;
	pertCnt = 0;
	restartCnt = 0;
	walksatCnt = 0;
	lastRestart = 0;
	flipsInRestart = 0;
	maxFlipsInRestart = 0;
	lastIntegerX.clear();
	chrono.reset();
	lpWatch.reset();
	roundWatch.reset();
	fixed.clear();
	binaries.clear();
	gintegers.clear();
	integers.clear();
	rows.clear();
	isPureInteger = false;
	isBinary = false;
	objOffset = 0.0;
	model = MIPModelPtr();
	closestPoint.clear();
	closestDist = INFBOUND;
	hasIncumbent = false;
}

void FeasibilityPump::init(MIPModelPtr _model, const std::vector<char>& ctype)
{
	DOMINIQS_ASSERT( _model );
	DOMINIQS_ASSERT( frac2int );
	// INIT
	consoleInfo("[fpInit]");
	reset();
	model = _model;
	int n = model->ncols();
	// if user passed column type information, used it
	if (ctype.size())
	{
		DOMINIQS_ASSERT(n == (int)ctype.size());
		for (int j = 0; j < n; j++)  model->ctype(j, ctype[j]);
	}
	frac2int->init(model, true);
	frac_x.resize(n, 0);
	integer_x.resize(n, 0);
	obj.resize(n, 0);
	lb.resize(n, 0);
	ub.resize(n, 0);
	xType.resize(n);
	closestPoint.clear();
	closestDist = INFBOUND;
	model->objcoefs(&obj[0]);
	objNorm = sqrt(dotProduct(&obj[0], &obj[0], n));
	model->lbs(&lb[0]);
	model->ubs(&ub[0]);
	model->ctypes(&xType[0]);
	for (int i = 0; i < n; i++)
	{
		// fixed variables are not considered integer variables, of any kind
		// this is correct if the rounding function does not alter their value!
		// this is trivially true for simple rounding, but some care must be
		// taken for more elaborate strategies!!!
		if ( equal(lb[i], ub[i], integralityEps) ) fixed.push_back(i);
		else if ( xType[i] != 'C' )
		{
			integers.push_back(i);
			if (xType[i] == 'B') binaries.push_back(i);
			else                 gintegers.push_back(i);
		}
	}
	// extract the rows
	int m = model->nrows();
	rows.resize(m);
	for (int i = 0; i < m; i++)
	{
		ConstraintPtr c = std::make_shared<Constraint>();
		model->row(i, c->row, c->sense, c->rhs, c->range);
		rows[i] = c;
	}

	isBinary = (gintegers.size() == 0);
	isPureInteger = (fixed.size() + integers.size() == (unsigned int)n);
	consoleLog("#cols = {} #bins = {} #integers = {}", n, binaries.size(), gintegers.size());
	consoleLog("fixedCnt = {} isBinary = {} isPureInteger = {}",
				fixed.size(), isBinary, isPureInteger);
	objOffset = model->objOffset();
	model->switchToLP();
}


bool FeasibilityPump::pump(const std::vector<double>& xStart, bool pFeas)
{
	DOMINIQS_ASSERT( model );
	DOMINIQS_ASSERT( frac2int );
	chrono.start();
	int n = model->ncols();
	primalFeas = false;
	ObjSense origObjSense = model->objSense();
	double dualBound = -static_cast<int>(origObjSense) * INFBOUND;
	primalBound = static_cast<int>(origObjSense) * INFBOUND;

	// setup iteration display
	display.addColumn("iter", 0, 6);
	display.addColumn("stage", 1, 6);
	display.addColumn("alpha", 2, 10, 4);
	display.addColumn("origObj", 3, 20);
	display.addColumn("projObj", 6, 15, 4);
	display.addColumn("dist", 8, 15, 4);
	display.addColumn("#frac", 9, 8);
	display.addColumn("P", 11, 3);
	display.addColumn("#flips", 12, 8);
	display.addColumn("lpiter", 13, 8);
	display.addColumn("time", 15, 10);

	// Ctrl-C handling
	model->handleCtrlC(true);

	// find first fractional solution (or use user supplied one)
	// this sets up frac_x
	if (xStart.size())
	{
		DOMINIQS_ASSERT( (int)xStart.size() == n );
		frac_x = xStart;
		primalFeas = pFeas;
	}
	else
	{
		solveInitialLP();
		if (primalFeas)  dualBound = getSolutionValue(frac_x);
	}
	consoleDebug(DebugLevel::Verbose, "startNumFrac = {}", solutionNumFractional(integers, frac_x, integralityEps));
	consoleLog("");


	// setup for changing objective function
	model->objSense(ObjSense::MIN); //< change obj sense: minimize distance
	model->objOffset(0.0); //< get rid of offset

	double runningAlpha = alpha;
	// If there is no objective, there is no need to use the objective FP
	if (isNull(objNorm))
	{
		objNorm = 1.0;
		runningAlpha = 0.0;
	}
	if (runningAlpha == 0.0)  display.setVisible("alpha", false);


	consoleInfo("[pump]");
	display.printHeader(std::cout);

	// stage 1
	bool found = pumpLoop(runningAlpha, 1);

	// stage 2 specific setup
	maxFlipsInRestart = std::max(int(gintegers.size() / 10.0), 10);
	consoleDebug(DebugLevel::Verbose, "maxFlipsInRestart = {}", maxFlipsInRestart);
	if (closestPoint.empty())
	{
		// no iterations were made in stage 1
		// we can still use frac_x
	}
	else
	{
		// use closest point as starting vector for the next pumping loop
		frac_x = closestPoint;
		primalFeas = false;
	}
	closestDist = INFBOUND;

	// stage 2
	// can skip stage 2 only if we have found a stage-1 solution and there are no general integers
	if (!found || gintegers.size())  found = pumpLoop(runningAlpha, 2);

	consoleLog("");

	// stage 3
	if (!found && doStage3)  found = stage3();

	if (found)  foundIncumbent(frac_x, getSolutionValue(frac_x));

	model->handleCtrlC(false);
	chrono.stop();

	// restore objective stuff
	model->objSense(origObjSense);
	model->objOffset(objOffset);

	model = MIPModelPtr();

	consoleLog("");
	consoleInfo("[results]");
	LOG_ITEM("primalBound", primalBound);
	LOG_ITEM("dualBound", dualBound);
	LOG_ITEM("objSense", origObjSense);
	LOG_ITEM("numSols", (int)found);
	LOG_ITEM("totalLpTime", lpWatch.getTotal());
	LOG_ITEM("totalRoundingTime", roundWatch.getTotal());
	LOG_ITEM("iterations", nitr);
	LOG_ITEM("rootTime", rootTime);
	LOG_ITEM("time", chrono.getTotal());
	LOG_ITEM("firstPerturbation", firstPerturbation);
	LOG_ITEM("perturbationCnt", pertCnt);
	LOG_ITEM("restartCnt", restartCnt);
	LOG_ITEM("walksatCnt", walksatCnt);
	return found;
}

// feasiblity pump helpers

void FeasibilityPump::solveInitialLP()
{
	consoleInfo("[initialSolve]");
	model->logging(true);
	double timeLeft = std::max(timeLimit - chrono.getElapsed(), 0.0);
	model->dblParam(DblParam::TimeLimit, timeLeft);
	model->lpopt(firstOptMethod);
	rootTime = chrono.getElapsed();
	int simplexIt = model->intAttr(IntAttr::SimplexIterations);
	int barrierIt = model->intAttr(IntAttr::BarrierIterations);
	rootLpIter = std::max(simplexIt, barrierIt);
	model->sol(&frac_x[0]);
	primalFeas = model->isPrimalFeas();
	model->logging(false);
	double dualBound = getSolutionValue(frac_x);
	consoleLog("Initial LP: lpiter={} barit={} time={:.4f} pfeas={} dualbound={:.2f}",
				simplexIt, barrierIt, rootTime, primalFeas, dualBound);
}


void FeasibilityPump::perturbe(std::vector<double>& x, bool ignoreGeneralIntegers)
{
	pertCnt++;

	std::multimap< double, int, std::greater<double> > toOrder;
	int nflips = avgFlips * (rnd.getFloat() + 0.5);
	int flipsDone = 0;

	// add fractional variables
	double sigma;
	if (ignoreGeneralIntegers)
	{
		for (int j: binaries)
		{
			sigma = fabs(x[j] - frac_x[j]);
			if (sigma > integralityEps) toOrder.emplace(sigma, j);
		}
	}
	else
	{
		for (int j: integers)
		{
			sigma = fabs(x[j] - frac_x[j]);
			if (sigma > integralityEps) toOrder.emplace(sigma, j);
		}
	}

	// compute number of flips to do
	int nFracFlips = std::min(nflips, (int) toOrder.size()); //number of vars to be flipped from fractional

	// add variables from walksat if needed
	if (walksatPerturbe && (nFracFlips < nflips))
	{
		int nneeded = nflips - nFracFlips;
		std::set<int> supp;

		// compute support of infeasible constraints
		infeasibleSupport(x, supp, ignoreGeneralIntegers);

		// remove variables already in toOrder
		std::multimap<double, int>::const_iterator itr = toOrder.begin();
		std::multimap<double, int>::const_iterator end = toOrder.end();
		while (itr != end)
		{
			int j = itr->second;
			supp.erase(j);
			++itr;
		}

		// pick randomly some elements from xsupp to add to toOrder
		std::vector<int> xsupp(supp.begin(), supp.end());
		std::shuffle(xsupp.begin(), xsupp.end(), std::default_random_engine(seed));
		std::vector<int>::const_iterator xitr = xsupp.begin();
		std::vector<int>::const_iterator xend = xsupp.end();
		while ((xitr != xend) && (nneeded > 0))
		{
			int j = *xitr++;
			toOrder.emplace(0.0, j);
			nneeded--;
		}

		walksatCnt++;
	}

	// do flips
	std::multimap<double, int>::const_iterator itr = toOrder.begin();
	std::multimap<double, int>::const_iterator end = toOrder.end();
	while ((itr != end) && (flipsDone < nflips))
	{
		int toFlip = itr->second;

		if (equal(x[toFlip], lb[toFlip], integralityEps)) { x[toFlip] += 1.0; ++flipsDone; }
		else if (equal(x[toFlip], ub[toFlip], integralityEps)) { x[toFlip] -= 1.0; ++flipsDone; }
		else
		{
			// this can happen only with general integer variables
			if (lessThan(x[toFlip], frac_x[toFlip], integralityEps)) { x[toFlip] += 1.0; ++flipsDone; }
			if (greaterThan(x[toFlip], frac_x[toFlip], integralityEps)) { x[toFlip] -= 1.0; ++flipsDone; }
		}
		++itr;
	}
	DOMINIQS_ASSERT( flipsDone );
	display.set("P", " *");
	display.set("#flips", flipsDone);
}


void FeasibilityPump::restart(std::vector<double>& x, bool ignoreGeneralIntegers)
{
	restartCnt++;
	// get previous solution
	DOMINIQS_ASSERT( lastIntegerX.size() );
	const std::vector<double>& previousSol = lastIntegerX.begin()->second;
	// perturbe
	double sigma;
	double r;
	// perturbe binaries
	int changed = 0;
	unsigned int size = binaries.size();
	for (unsigned int i = 0; i < size; i++)
	{
		int j = binaries[i];
		r = rnd.getFloat() - 0.47;
		if ( r > 0 && equal(x[j], previousSol[j], integralityEps) )
		{
			sigma = fabs(x[j] - frac_x[j]);
			if ( (sigma + r) > 0.5 )
			{
				x[j] = isNull(x[j], integralityEps) ? 1.0 : 0.0;
				++changed;
			}
		}
	}

	if (!ignoreGeneralIntegers && gintegers.size())
	{
		if (lastRestart < nitr)
		{
			while (lastRestart++ < nitr) flipsInRestart *= GEOM_FACTOR; // geometric decrease
		}
		flipsInRestart = std::min(flipsInRestart + 2 * avgFlips + 1, maxFlipsInRestart); // linear increase
		DOMINIQS_ASSERT( flipsInRestart );
		// perturbe general integers
		double newValue;
		for (int i = 0; i < flipsInRestart; i++)
		{
			int randIdx = int(rnd.getFloat() * (gintegers.size() - 1));
			DOMINIQS_ASSERT( randIdx >= 0 && randIdx < (int)gintegers.size() );
			int j = gintegers[randIdx];
			double xlb = lb[i];
			double xub = ub[i];
			r = rnd.getFloat();
			if( (xub - xlb) < BIGBIGM ) newValue = floor(xlb + (1 + xub - xlb) * r);
			else if( (x[j] - xlb) < BIGM ) newValue = xlb + (2 * BIGM - 1) * r;
			else if( (xub - x[j]) < BIGM ) newValue = xub - (2 * BIGM - 1) * r;
			else newValue = x[j] + (2 * BIGM - 1) * r - BIGM;
			newValue = std::max(std::min(newValue, xub), xlb);
			if( different(newValue, x[j], integralityEps) )
			{
				x[j] = newValue;
				++changed;
			}
		}
		DOMINIQS_ASSERT( changed );
	}
	else
	{
		// there are no general integers (or we are in stage1), so we only flipped binary variables.
		// This however can lead to an infinite loop, for example in the case in which we flip all of them:
		// then two consecutive integer solutions have nothing in common, and we end up flipping nothing...
		// As a safety net, we do a big random shake if the changed counter is still zero at this point.
		if (!changed)
		{
			for (unsigned int i = 0; i < size; i++)
			{
				int j = binaries[i];
				r = rnd.getFloat();
				if ( r > 0.5 )
				{
					x[j] = isNull(x[j], integralityEps) ? 1.0 : 0.0;
					++changed;
				}
			}
			DOMINIQS_ASSERT( changed );
		}
	}
	display.set("P", "**");
	display.set("#flips", changed);
}


bool FeasibilityPump::pumpLoop(double& runningAlpha, int stage)
{
	// setup
	lastIntegerX.clear();
	int n = model->ncols();
	std::vector<double> distObj(n, 0);
	std::vector<int> colIndices(n);
	std::iota(colIndices.begin(), colIndices.end(), 0);
	int oldIterCnt = nitr;
	bool ignoreGenerals = (stage == 1) ? true : false;
	const auto& intSubset = (stage == 1) ? binaries : integers;
	frac2int->ignoreGeneralIntegers(ignoreGenerals);
	double pumpTimeLimit = (timeMult > 0.0) ? timeMult*rootTime : std::numeric_limits<double>::max();
	std::vector<std::string> xNames;
	model->colNames(xNames);
	int lpIterLimit = -1;
	if (lpIterMult > 0.0)
	{
		lpIterLimit = int(rootLpIter * lpIterMult);
		lpIterLimit = std::max(lpIterLimit, 10);
	}

	while (!model->aborted()
		&& ((nitr - oldIterCnt) < stageIterLimit)
		&& (nitr < iterLimit))
	{
		// check if frac_x is feasible (w.r.t. the integer variables in this stage)
		bool found = (primalFeas && isSolutionInteger(intSubset, frac_x, integralityEps));
		if (found)
		{
			// update closest point
			closestDist = 0.0;
			closestPoint = frac_x;
			// then break this stage
			break;
		}

		// global timelimit check
		double timeLeft = std::max(std::min(timeLimit, pumpTimeLimit) - chrono.getElapsed(), 0.0);
		if (timeLeft <= 0.0)  break;

		// display logger
		nitr++;
		display.resetIteration();
		if (display.needHeader(nitr)) display.printHeader(std::cout);

		// frac -> int
		roundWatch.start();
		frac2int->apply(frac_x, integer_x);
		roundWatch.stop();
		consoleDebug(DebugLevel::Verbose, "roundingTime = {}", roundWatch.getPartial());

		// cycle detection and antistalling actions
		// is it the same of the last one? If yes perturbe
		if (lastIntegerX.size()
			&& areSolutionsEqual(intSubset, integer_x, (*(lastIntegerX.begin())).second, integralityEps)
			&& equal(runningAlpha, (*(lastIntegerX.begin())).first, alphaDist))
		{
			if (!pertCnt)  firstPerturbation = nitr;
			perturbe(integer_x, ignoreGenerals);
		}
		// do a restart until we are able to insert it in the cache
		for (int rtry = 0; rtry < 10; rtry++)
		{
			if (isInCache(runningAlpha, integer_x, ignoreGenerals)) restart(integer_x, ignoreGenerals);
			else break;
		}
		lastIntegerX.push_front(AlphaVector(runningAlpha, integer_x));

		// int -> frac
		lpWatch.start();

		double thisAlpha = runningAlpha;
		// if the distance function is not the pure distance one
		// then we might not realize the current integer_x is feasible.
		// so we explictly check for feasibility and, if so,
		// temporarily set the running alpha to zero.
		if (isSolutionFeasible(integer_x, rows))  thisAlpha = 0.0;

		// setup distance objective
		int addedVars = 0;
		int addedConstrs = 0;
		std::fill(distObj.begin(), distObj.end(), 0.0);
		for (int j: binaries)
		{
			double fracj = std::min(1.0, std::max(frac_x[j], 0.0)); //< clip fractional value to [0,1]
			if (isNull(integer_x[j], integralityEps))
			{
				double distCoef = 1.0;
				if (penaltyObj)  distCoef = 1.0 / std::max(1.0 - fracj, 1e-6);
				distObj[j] = distCoef;
			}
			else
			{
				double distCoef = -1.0;
				if (penaltyObj)  distCoef = -1.0 / std::max(fracj, 1e-6);
				distObj[j] = distCoef;
			}
		}
		if (stage > 1)
		{
			for (int j: gintegers)
			{
				// TODO: penalty objective for general integers?
				if (equal(integer_x[j], lb[j], integralityEps))
				{
					distObj[j] = 1.0;
				}
				else if (equal(integer_x[j], ub[j], integralityEps))
				{
					distObj[j] = -1.0;
				}
				else
				{
					// add auxiliary variable
					std::string deltaName = xNames[j] + "_delta";
					model->addEmptyCol(deltaName, 'C', 0.0, INFBOUND, 0.0);
					int auxIdx = model->ncols() - 1;
					colIndices.push_back(auxIdx);
					distObj.push_back(1.0);
					addedVars++;
					// add constraints
					SparseVector vec;
					vec.push(j, 1.0);
					vec.push(auxIdx, -1.0);
					model->addRow(xNames[j] + "_d1", vec.idx(), vec.coef(), 2, 'L', integer_x[j]);
					addedConstrs++;
					vec.coef()[1] = 1.0;
					model->addRow(xNames[j] + "_d2", vec.idx(), vec.coef(), 2, 'G', integer_x[j]);
					addedConstrs++;
				}
			}
			consoleDebug(DebugLevel::Verbose, "addedVars={} addedConstrs={}", addedVars, addedConstrs);
			DOMINIQS_ASSERT( distObj.size() == (unsigned int)(n + addedVars) );
			DOMINIQS_ASSERT( distObj.size() == colIndices.size() );
		}

		// randomize distance coefficients
		if (randomizeLP)
		{
			for (int j = 0; j < (n + addedVars); j++)
			{
				double randMult = rnd.getFloat() * 0.1 + 0.9; //< random float in [0.9,1.0)
				distObj[j] *= randMult;
			}
		}

		// objective FP
		if (thisAlpha > 0.0)
		{
			// compute distance norm and scale distance objective by (1-thisAlpha)
			double distNorm = 0.0;
			for (int j = 0; j < (n + addedVars); j++)
			{
				distNorm += (distObj[j]*distObj[j]);
				distObj[j] *= (1.0 - thisAlpha);
			}
			distNorm = sqrt(distNorm) * (1.0 - thisAlpha);

			// add objective with proper weight
			double weight = (thisAlpha * distNorm) / objNorm;
			accumulate(&distObj[0], &obj[0], n, weight);
		}

		// set objective
		model->objcoefs(colIndices.size(), &colIndices[0], &distObj[0]);

		// solve LP
		if (lpIterLimit > 0)  model->intParam(IntParam::IterLimit, lpIterLimit);
		model->dblParam(DblParam::TimeLimit, timeLeft);
		model->lpopt(reOptMethod);
		lpWatch.stop();

		// get solution
		model->sol(&frac_x[0], 0, n-1);
		primalFeas = model->isPrimalFeas();
		consoleDebug(DebugLevel::VeryVerbose, "Iteration {}: time={} pFeas={} lpiter={}",
				lpWatch.getPartial(), primalFeas, model->intAttr(IntAttr::SimplexIterations));
		double projObj = model->objval();

		// cleanup added vars and constraints
		if (addedConstrs)
		{
			int begin = model->nrows() - addedConstrs;
			model->delRows(begin, begin + addedConstrs - 1);
		}
		if (addedVars)
		{
			int begin = model->ncols() - addedVars;
			model->delCols(begin, begin + addedVars - 1);
		}
		colIndices.resize(n);
		distObj.resize(n);
		DOMINIQS_ASSERT( model->ncols() == n );

		// get some statistics
		double origObj = dotProduct(&obj[0], &frac_x[0], n) + objOffset;
		double dist = solutionsDistance(intSubset, frac_x, integer_x);
		int numFrac = solutionNumFractional(intSubset, frac_x, integralityEps);

		// save integer_x as best point if distance decreased
		if (dist < closestDist)
		{
			closestDist = dist;
			closestPoint = integer_x;
		}

		// display log
		if (display.needPrint(nitr))
		{
			int simplexIt = model->intAttr(IntAttr::SimplexIterations);
			int barrierIt = model->intAttr(IntAttr::BarrierIterations);
			display.set("stage", stage);
			display.set("iter", nitr);
			display.set("alpha", runningAlpha);
			display.set("origObj", origObj);
			display.set("time", chrono.getElapsed());
			display.set("dist", dist);
			display.set("#frac", numFrac);
			display.set("projObj", projObj);
			display.set("lpiter", std::max(simplexIt, barrierIt));
			display.printIteration(std::cout);
		}

		// update running alpha
		runningAlpha *= alphaFactor;
		if (runningAlpha <= 1e-4)  runningAlpha = 0.0;
	}

	return (primalFeas && isSolutionInteger(integers, frac_x, integralityEps));
}


bool FeasibilityPump::stage3()
{
	if (model->aborted()) return false;
	if (closestPoint.empty()) return false;
	consoleInfo("[stage3]");
	double elapsedTime = chrono.getElapsed();
	double remainingTime = timeLimit - elapsedTime;
	double s3TimeLimit = std::max(std::min(remainingTime, elapsedTime), 1.0);
	if (lessThan(remainingTime, 0.1)) return false;

	consoleLog("Starting stage3 from point with distance={} [timeLimit={}]", closestDist, s3TimeLimit);

	int n = model->ncols();
	bool found = false;
	std::vector<std::string> xNames;
	model->colNames(xNames);

	// restore type information
	std::vector<char> ctype(n, 'C');
	for (int j: binaries) ctype[j] = 'B';
	for (int j: gintegers) ctype[j] = 'I';
	for (int j = 0; j < n; j++)  model->ctype(j, ctype[j]);

	// load best point and generate objective
	integer_x = closestPoint;
	int addedVars = 0;
	int addedConstrs = 0;
	std::vector<double> distObj(n, 0.0);
	std::vector<int> colIndices(n);
	std::iota(colIndices.begin(), colIndices.end(), 0);
	for (int j: binaries) distObj[j] = (isNull(integer_x[j], integralityEps) ? 1.0 : 1.0);
	for (int j: gintegers)
	{
		if (equal(integer_x[j], lb[j], integralityEps)) distObj[j] = 1.0;
		else if (equal(integer_x[j], ub[j], integralityEps)) distObj[j] = -1.0;
		else
		{
			// add auxiliary variable
			std::string deltaName = xNames[j] + "_delta";
			model->addEmptyCol(deltaName, 'C', 0.0, INFBOUND, 0.0);
			int auxIdx = model->ncols() - 1;
			colIndices.push_back(auxIdx);
			distObj.push_back(1.0);
			addedVars++;
			// add constraints
			SparseVector vec;
			vec.push(j, 1.0);
			vec.push(auxIdx, -1.0);
			model->addRow(xNames[j] + "_d1", vec.idx(), vec.coef(), 2, 'L', integer_x[j]);
			addedConstrs++;
			vec.coef()[1] = 1.0;
			model->addRow(xNames[j] + "_d2", vec.idx(), vec.coef(), 2, 'G', integer_x[j]);
			addedConstrs++;
		}
	}
	consoleDebug(DebugLevel::Normal, "addedVars={} addedConstrs={}", addedVars, addedConstrs);
	DOMINIQS_ASSERT( distObj.size() == (unsigned int)(n + addedVars) );
	DOMINIQS_ASSERT( distObj.size() == colIndices.size() );
	model->objcoefs(colIndices.size(), &colIndices[0], &distObj[0]);

	model->logging(true);
	model->intParam(IntParam::SolutionLimit, 1);
	model->dblParam(DblParam::TimeLimit, timeLimit);
	model->mipopt();
	primalFeas = model->isPrimalFeas();
	model->logging(false);
	if (primalFeas)
	{
		model->sol(&frac_x[0], 0, n-1);
		DOMINIQS_ASSERT( isSolutionInteger(integers, frac_x, integralityEps) );
		found = true;
	}

	// cleanup added vars and constraints and restore obj
	if (addedConstrs)
	{
		int begin = model->nrows() - addedConstrs;
		model->delRows(begin, begin + addedConstrs - 1);
	}
	if (addedVars)
	{
		int begin = model->ncols() - addedVars;
		model->delCols(begin, begin + addedVars - 1);
	}

	// restore original objective function
	model->objcoefs(model->ncols(), &colIndices[0], &obj[0]);

	return found;
}

void FeasibilityPump::foundIncumbent(const std::vector<double>& x, double objval)
{
	incumbent = x;
	primalBound = objval;
	hasIncumbent = true;
	frac2int->newIncumbent(incumbent, primalBound);
	lastIntegerX.clear();
}

bool FeasibilityPump::isInCache(double a, const std::vector<double>& x, bool ignoreGeneralIntegers)
{
	bool found = false;
	std::list< AlphaVector >::const_iterator itr = lastIntegerX.begin();
	std::list< AlphaVector >::const_iterator end = lastIntegerX.end();
	if (ignoreGeneralIntegers)
	{
		while ((itr != end) && !found)
		{
			if ((fabs(a - itr->first) < alphaDist) && areSolutionsEqual(binaries, itr->second, integer_x, integralityEps)) found = true;
			++itr;
		}
	}
	else
	{
		while ((itr != end) && !found)
		{
			if ((fabs(a - itr->first) < alphaDist) && areSolutionsEqual(integers, itr->second, integer_x, integralityEps)) found = true;
			++itr;
		}
	}
	return found;
}

void FeasibilityPump::infeasibleSupport(const std::vector<double>& x, std::set<int>& supp, bool ignoreGeneralIntegers)
{
	if (supp.empty())
	{
		// find set of infeasible constraints
		std::set<int> infeas;
		int m = model->nrows();
		for (int i = 0; i < m; i++) {
			ConstraintPtr c = rows[i];
			if (!c->satisfiedBy(&x[0])) infeas.insert(i);
		}
		// find support of infeasible constraints
		supp.clear();
		for (auto i: infeas)
		{
			ConstraintPtr c = rows[i];
			const int* idx = c->row.idx();
			for (int k = 0; k < c->row.size(); k++)
			{
				int j = idx[k];
				if ((xType[j] == 'B') || ((xType[j] == 'I') && !ignoreGeneralIntegers)) supp.insert(j);
			}
		}
	}
}

} // namespace dominiqs
