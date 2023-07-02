/**
 * @file feaspump.cpp
 * @brief Feasibility Pump algorithm Source
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * @author Gioni Mexi <gionimexi at gmail dot com>
 * 2023
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
static double solutionsLinfDistance(const std::vector<int>& integers, const std::vector<double>& x1, const std::vector<double>& x2)
{
	double dist = 0.0;
	for (int j: integers) dist = std::max(dist, fabs(x1[j] - x2[j]));
	return dist;
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

static std::vector<double> exponDecayVector(const double factor,const int numPoints)
{
		// compute the sum of all elements -> used to make a convex combination
		double sumAll = 0.0;

		// resVector contains the barycentric coordinates of the n points
		std::vector<double> resVector;

		double currElement = 1.0;
		resVector.push_back(currElement);
		sumAll += currElement;
		for (int i = 1; i < numPoints; i++)
		{
			currElement = currElement * factor;
			resVector.push_back(currElement);
			sumAll += currElement;

		}
		// make the sum of all vector elements convex
		for (int i = 0; i < numPoints; i++)
		{
			resVector[i] = resVector[i] / sumAll;
		}

		DOMINIQS_ASSERT( resVector.size() == numPoints);
		return resVector;
}

static std::vector<double> harmonicVector(const int numPoints)
{
		// compute the sum of all elements -> used to make a convex combination
		double sumAll = 0.0;

		// resVector contains the barycentric coordinates of the n points
		std::vector<double> resVector;

		double currElement = 1.0;		
		resVector.push_back(currElement);
		sumAll += currElement;
		for (int i = 1; i < numPoints; i++)
		{
			currElement = 1.0 / (i+1);
			resVector.push_back(currElement);
			sumAll += currElement;

		}
		// make the sum of all vector elements convex
		for (int i = 0; i < numPoints; i++)
		{
			resVector[i] = resVector[i] / sumAll;
		}

		DOMINIQS_ASSERT( resVector.size() == numPoints);
		return resVector;
}


namespace dominiqs {

static const double DEF_TIME_LIMIT = 3600.0;
static const double DEF_TIME_MULT = 100.0;
static const double DEF_LPITER_MULT = -1.0;
static const int DEF_STAGE_1_ITER_LIMIT = 50;
static const int DEF_STAGE_2_ITER_LIMIT = 50;
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
static const bool DEF_EXPOBJ = false;
static const bool DEF_LOGISOBJ = false;
static const bool DEF_ACFP = false;

static const char DEF_FIRST_OPT_METHOD = 'S';
static const char DEF_REOPT_METHOD = 'S';


static const double GEOM_FACTOR = 0.85;
static const double BIGM = 1e9;
static const double BIGBIGM = 1e15;
static const double INFBOUND = 1e20;

// using more than one integer in the objective possible
static const int DEF_INTEGERS_OBJ = 1;
static const double DEF_INTEGER_SCALE_OBJ = 1.0;
static const int DEF_LESS_VIOLATED_INT = 0;
static const int DEF_LESS_DISTANCE_INT = 0;
static const int DEF_BEST_OBJ_INT = 0;
static const int DEF_HARMONIC_WEIGHTS = 0;
static const int DEF_EXP_DECAY_WEIGHTS = 0;

// using more than one integer in the stage 3 objective is possible
static const int DEF_INTEGERS_OBJ_STAGE_3 = 1;
static const double DEF_INTEGER_SCALE_STAGE_3 = 1.0;
static const int DEF_LESS_VIOLATED_INT_STAGE_3 = 0;
static const int DEF_BEST_OBJ_INT_STAGE_3 = 0;
static const int DEF_HARMONIC_WEIGHTS_STAGE_3 = 0;

// use more integers in stage3 objective
static const bool DEF_NEW_STAGE_3 = false;
static const bool DEF_RENS_STAGE_3 = false;
static const bool DEF_MULTIRENS_STAGE_3 = false;
static const bool DEF_NORMAL_MIP_STAGE_3 = false;
static const bool DEF_RENS_CLOSEST_DIST_STAGE_3 = false;


// number of fractional numbers to be considered in objective (not working)
static const int DEF_FRACS_OBJ = 1;

// barycenter coordinates either for the aggregated fractional point or the objective when using multiple integers
static const double DEF_FRAC_SCALE_FACTOR = 1.0;

// TODO: custom barycenter coordinates. For example if we want to use the first and third integers with
// coordinates 0.75 and 0.25 the coordinate vector should be defined as e.g. {0.75, 0 , 0.25}

// alternative objective scale
static const bool DEF_NEW_OBJ_SCALE = false;
static const bool DEF_SCALE_C_AND_DELTA = false;

// for the KK parameter in stage 1 and 2
static const int ITR1_NO_IMPR = 70;
static const int ITR2_NO_IMPR = 600;

//accumulate integers
static const bool DEF_AGGREGATE_INT = false;

static const double DEF_PDLP_TOLERANCE = 1.0e-6;
static const double DEF_PDLP_TOLERANCE_DECREASE = 1.0;

static const bool DEF_PDLP_WARMSTART = false;

FeasibilityPump::FeasibilityPump() :
	timeLimit(DEF_TIME_LIMIT), timeMult(DEF_TIME_MULT), lpIterMult(DEF_LPITER_MULT),
	stage1IterLimit(DEF_STAGE_1_ITER_LIMIT),stage2IterLimit(DEF_STAGE_2_ITER_LIMIT),
	iterLimit(DEF_ITER_LIMIT), avgFlips(DEF_AVG_FLIPS), integralityEps(DEF_INTEGRALITY_EPS),
	seed(DEF_SEED), alpha(DEF_ALPHA), alphaFactor(DEF_ALPHA_FACTOR), alphaDist(DEF_ALPHA_DIST),
	doStage3(DEF_DO_STAGE_3), walksatPerturbe(DEF_WALKSAT_PERTURBE),
	randomizeLP(DEF_RANDOMIZE_LP), penaltyObj(DEF_PENALTYOBJ),
	firstOptMethod(DEF_FIRST_OPT_METHOD), reOptMethod(DEF_REOPT_METHOD),
	objOffset(0.0), hasIncumbent(false), numIntegersObj(DEF_INTEGERS_OBJ),numFracsObj(DEF_FRACS_OBJ), 
	newScaleC(DEF_NEW_OBJ_SCALE),newScaleDelta(DEF_SCALE_C_AND_DELTA),itr1NoImpr(ITR1_NO_IMPR),
	itr2NoImpr(ITR2_NO_IMPR), fracScaleFactor(DEF_FRAC_SCALE_FACTOR), intScaleFactor(DEF_INTEGER_SCALE_OBJ),
	aggregateInts(DEF_AGGREGATE_INT),newStage3(DEF_NEW_STAGE_3), stage3ScaleFactor(DEF_INTEGER_SCALE_STAGE_3),
	stage3IntegersObj(DEF_INTEGERS_OBJ_STAGE_3),stage3lessViolatedIntegers(DEF_LESS_VIOLATED_INT_STAGE_3),
	stage3bestObjIntegers(DEF_BEST_OBJ_INT_STAGE_3), stage3harmonicWeights(DEF_HARMONIC_WEIGHTS_STAGE_3),
	expObj(DEF_EXPOBJ), logisObj(DEF_LOGISOBJ), analcenterFP(DEF_ACFP),
	lessViolatedIntegers(DEF_LESS_VIOLATED_INT), bestObjIntegers(DEF_BEST_OBJ_INT), lessDistanceIntegers(DEF_LESS_DISTANCE_INT),
	harmonicWeights(DEF_HARMONIC_WEIGHTS), exponDecayWeights(DEF_EXP_DECAY_WEIGHTS),
	rensStage3(DEF_RENS_STAGE_3),multirensStage3(DEF_MULTIRENS_STAGE_3), normalMIPStage3(DEF_NORMAL_MIP_STAGE_3),
	rensClosestDistStage3(DEF_RENS_CLOSEST_DIST_STAGE_3), pdlpTol(DEF_PDLP_TOLERANCE),
	pdlpTolDecreaseFactor(DEF_PDLP_TOLERANCE_DECREASE), pdlpWarmStart(DEF_PDLP_WARMSTART)
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
	// READ_FROM_CONFIG( timeMult, DEF_TIME_MULT );
	timeMult = 0.0;
	
	READ_FROM_CONFIG( lpIterMult, DEF_LPITER_MULT );
	READ_FROM_CONFIG( stage1IterLimit, DEF_STAGE_1_ITER_LIMIT );
	READ_FROM_CONFIG( stage2IterLimit, DEF_STAGE_2_ITER_LIMIT );

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

	// new options
	READ_FROM_CONFIG( numIntegersObj, DEF_INTEGERS_OBJ );
	READ_FROM_CONFIG( numFracsObj, DEF_FRACS_OBJ );
	READ_FROM_CONFIG( lessViolatedIntegers, DEF_LESS_VIOLATED_INT );
	READ_FROM_CONFIG( lessDistanceIntegers, DEF_LESS_DISTANCE_INT );
	READ_FROM_CONFIG( bestObjIntegers, DEF_BEST_OBJ_INT );
	READ_FROM_CONFIG( harmonicWeights, DEF_HARMONIC_WEIGHTS );
	READ_FROM_CONFIG( exponDecayWeights, DEF_EXP_DECAY_WEIGHTS );

	READ_FROM_CONFIG( fracScaleFactor, DEF_FRAC_SCALE_FACTOR );
	READ_FROM_CONFIG( intScaleFactor, DEF_INTEGER_SCALE_OBJ );
	READ_FROM_CONFIG( newScaleC, DEF_NEW_OBJ_SCALE );
	READ_FROM_CONFIG( newScaleDelta, DEF_SCALE_C_AND_DELTA );
	READ_FROM_CONFIG( aggregateInts, DEF_AGGREGATE_INT );
	READ_FROM_CONFIG( expObj, DEF_EXPOBJ );
	READ_FROM_CONFIG( logisObj, DEF_LOGISOBJ );
	READ_FROM_CONFIG( analcenterFP, DEF_ACFP );

	READ_FROM_CONFIG( newStage3, DEF_NEW_STAGE_3 );
	READ_FROM_CONFIG( rensStage3, DEF_RENS_STAGE_3 );
	READ_FROM_CONFIG( multirensStage3, DEF_MULTIRENS_STAGE_3 );
	READ_FROM_CONFIG( normalMIPStage3, DEF_NORMAL_MIP_STAGE_3 );
	READ_FROM_CONFIG( rensClosestDistStage3, DEF_RENS_CLOSEST_DIST_STAGE_3 );

	READ_FROM_CONFIG( stage3ScaleFactor, DEF_INTEGER_SCALE_STAGE_3 );
	READ_FROM_CONFIG( stage3IntegersObj, DEF_INTEGERS_OBJ_STAGE_3 );
	READ_FROM_CONFIG( stage3lessViolatedIntegers, DEF_LESS_VIOLATED_INT_STAGE_3 );
	READ_FROM_CONFIG( stage3bestObjIntegers, DEF_BEST_OBJ_INT_STAGE_3 );
	READ_FROM_CONFIG( stage3harmonicWeights, DEF_HARMONIC_WEIGHTS_STAGE_3 );

	READ_FROM_CONFIG( pdlpTol, DEF_PDLP_TOLERANCE );
	READ_FROM_CONFIG( pdlpTolDecreaseFactor, DEF_PDLP_TOLERANCE_DECREASE );
	READ_FROM_CONFIG( pdlpWarmStart, DEF_PDLP_WARMSTART );


	//till here
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
	LOG_CONFIG( stage1IterLimit );
	LOG_CONFIG( stage2IterLimit );
	LOG_CONFIG( avgFlips );
	LOG_CONFIG( integralityEps );
	LOG_CONFIG( seed );

	itr1NoImpr =  stage1IterLimit;
	itr2NoImpr =  stage2IterLimit;


	// added
	LOG_CONFIG( numIntegersObj );
	LOG_CONFIG( lessViolatedIntegers );
	LOG_CONFIG( lessDistanceIntegers );
	LOG_CONFIG( bestObjIntegers );
	LOG_CONFIG( harmonicWeights );
	LOG_CONFIG( exponDecayWeights );
	LOG_CONFIG( intScaleFactor );
	LOG_CONFIG( numFracsObj );
	LOG_CONFIG( fracScaleFactor );
	LOG_CONFIG( newScaleC );
	LOG_CONFIG( newScaleDelta );
	LOG_CONFIG( aggregateInts );
	LOG_CONFIG( newStage3 );
	LOG_CONFIG( rensStage3 );
	LOG_CONFIG( multirensStage3 );
	LOG_CONFIG( normalMIPStage3 );
	LOG_CONFIG( rensClosestDistStage3 );

	LOG_CONFIG( stage3ScaleFactor );
	LOG_CONFIG( stage3IntegersObj );
	LOG_CONFIG( stage3lessViolatedIntegers );
	LOG_CONFIG( stage3bestObjIntegers );
	LOG_CONFIG( stage3harmonicWeights );

	LOG_CONFIG( expObj );
	LOG_CONFIG( logisObj );
	LOG_CONFIG( analcenterFP );

	LOG_CONFIG( pdlpTol );
	LOG_CONFIG( pdlpTolDecreaseFactor );
	LOG_CONFIG( pdlpWarmStart );

	// till here

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

	lastFracX.clear();
	closestIntegerXs1.clear();
	closestIntegerXs2.clear();
	multipleIntegerX.clear();

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
	closestFrac.clear();
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
	closestFrac.clear();
	closestDist = INFBOUND;
	model->objcoefs(&obj[0]);

	objNorm = sqrt(dotProduct(&obj[0], &obj[0], n));
	LOG_ITEM("objNorm", objNorm);

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
	newlb.resize(n,100000000000.0); /**< new lower bounds */
	newub.resize(n,-100000000000.0); /**< new upper bounds */

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
	display.addColumn("PDLP status", 20, 15);
	display.addColumn("PDLP feas", 21, 15);
	display.addColumn("PDLP iter", 25, 15);

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
	if (rensStage3 || multirensStage3)
	{
		for (int j: integers)
		{	
			if (equal(frac_x[j], floor(frac_x[j]), integralityEps))
			{
				newlb[j] = std::min(newlb[j], floor(frac_x[j]));
				newub[j] = std::max(newub[j], floor(frac_x[j]));
			}
			else 
			{
				newlb[j] = std::min(newlb[j], floor(frac_x[j]));
				newub[j] = std::max(newub[j], ceil(frac_x[j]));
			}
		}
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

	int stage = 1;
	// stage 1
	double timeModelS1 = 0.0;
	double timeModelS2 = 0.0;
	bool found = pumpLoop(runningAlpha, stage, dualBound, timeModelS1);

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

	// stage 2
	// can skip stage 2 only if we have found a stage-1 solution and there are no general integers
	if (!found || gintegers.size())
	{
		closestDist = INFBOUND;
		stage = 2;
		found = pumpLoop(runningAlpha, stage, dualBound, timeModelS2);
	}
	consoleLog("");

	// stage 3
	if (!found && doStage3)
	{
		stage = 3;
		found = stage3();
	}

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
	LOG_ITEM("stage", stage);
	LOG_ITEM("objSense", origObjSense);
	LOG_ITEM("numSols", (int)found);
	LOG_ITEM("totalLpTime", lpWatch.getTotal());
	LOG_ITEM("totalRoundingTime", roundWatch.getTotal());
	LOG_ITEM("iterations", nitr);
	LOG_ITEM("rootTime", rootTime);
	if (analcenterFP) 	LOG_ITEM("acTime", acTime);
	LOG_ITEM("time", chrono.getTotal());
	LOG_ITEM("timeModel", timeModelS1 + timeModelS2);
	LOG_ITEM("firstPerturbation", firstPerturbation);
	LOG_ITEM("perturbationCnt", pertCnt);
	LOG_ITEM("restartCnt", restartCnt);
	LOG_ITEM("walksatCnt", walksatCnt);
	if (stage == 3) LOG_ITEM("stage3Time", stage3Time);
	return found;
}

// feasiblity pump helpers

void FeasibilityPump::solveInitialLP()
{

	double timeLeft;
	if (analcenterFP)
	{
		consoleInfo("[computeAC]");
		int n = model->ncols();
		ac_x.resize(n, 0);
		timeLeft = std::max(timeLimit - chrono.getElapsed(), 0.0);
		model->dblParam(DblParam::TimeLimit, timeLeft);
		computeAC(ac_x);
		acTime = chrono.getElapsed();
	}

	consoleInfo("[initialSolve]");
	double elapsedBefore = chrono.getElapsed();
	timeLeft = std::max(timeLimit - chrono.getElapsed(), 0.0);
	model->dblParam(DblParam::TimeLimit, timeLeft);
	double time_model = model->lpopt(firstOptMethod, false, true);
	rootTime = chrono.getElapsed() - elapsedBefore;
	int simplexIt = model->intAttr(IntAttr::SimplexIterations);
	int barrierIt = model->intAttr(IntAttr::BarrierIterations);
	int pdlpIt = model->intAttr(IntAttr::PDLPIterations);

	rootLpIter = std::max(simplexIt, barrierIt);
	model->sol(&frac_x[0]);
	primalFeas = (model->isPrimalFeas() && isSolutionFeasible(frac_x, rows));
	double dualBound = -static_cast<int>(model->objSense()) * INFBOUND;
	if (primalFeas) dualBound = getSolutionValue(frac_x);
	consoleLog("Initial LP: lpiter={} barit={} pdlp={} time={:.4f} pfeas={} dualbound={:.2f}",
				simplexIt, barrierIt, pdlpIt, rootTime, primalFeas, dualBound);
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


bool FeasibilityPump::pumpLoop(double& runningAlpha, int stage, double& dualBound, double& timeModel)
{
	// setup
	lastIntegerX.clear();
	multipleIntegerX.clear();
	lastFracX.clear();
	int n = model->ncols();
	std::vector<double> distObj(n, 0);
	std::vector<int> colIndices(n);
	std::iota(colIndices.begin(), colIndices.end(), 0);
	int oldIterCnt = nitr;
	int stageIterLimit = (stage == 1) ? stage1IterLimit : stage2IterLimit;
	bool ignoreGenerals = (stage == 1) ? true : false;
	const auto& intSubset = (stage == 1) ? binaries : integers;
	frac2int->ignoreGeneralIntegers(ignoreGenerals);
	double pumpTimeLimit = (timeMult > 0.0) ? timeMult*rootTime : std::numeric_limits<double>::max();
	std::vector<std::string> xNames;
	model->colNames(xNames);
	int lpIterLimit = -1;

	int iterationsNoImpr = 0;
	int maxItrNoImpr = (stage == 1) ? itr1NoImpr : itr2NoImpr;
	bool usedOrigFpNoRestart = false;  // when using multiple reference points, if a cycle occurs use original FP instead of restart

	if (lpIterMult > 0.0)
	{
		lpIterLimit = int(rootLpIter * lpIterMult);
		lpIterLimit = std::max(lpIterLimit, 10);
	}
	bool lpfeasible = isSolutionFeasible(frac_x, rows);

	// Disable PdLP warm start in stage 2 because of numerical issues
	if (stage == 2)
		model->intParam(IntParam::PdlpWarmStart, 0);
	
	timeModel = 0;

	while (!model->aborted()
		&& ((nitr - oldIterCnt) < stageIterLimit)
		&& (nitr < iterLimit))
	{
		lpfeasible = isSolutionFeasible(frac_x, rows);
		bool applyPdlpRestart = (!lpfeasible && isSolutionInteger(intSubset, frac_x, integralityEps));
		// ToDo check the value of primFeas
		// check if frac_x is feasible (w.r.t. the integer variables in this stage)
		bool found = (lpfeasible && isSolutionInteger(intSubset, frac_x, integralityEps));
		if (found || (stage == 1 && binaries.size() == 0 ))
		{
			// update closest point
			closestDist = 0.0;
			closestPoint = frac_x;
			closestFrac = frac_x;
			// then break this stage
			return true;
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

		// If we want to round a convex combination of fractional points instead of the last one
		if (numFracsObj != 1 ) // aggregate fracs and then round
		{
			std::vector<double> frac2round; // convex combination of fractional points
			frac2round.resize(n, 0);
			lastFracX.push_front(NumberVector(nitr,frac_x));
			int numPoints = std::min( numFracsObj, (int) lastFracX.size());
			//scaleVec contains the barycentric coordinates of the new point
			std::vector<double> scaleVec;
			scaleVec = exponDecayVector(fracScaleFactor,numPoints);
			// get new fractional point
			aggregateFracs( frac2round, scaleVec);
			// get rounding
			frac2int->apply(frac2round,integer_x);		
		}
		// get rounding by using AC  
		else if (analcenterFP) 
		{
			double bestgamma;
			integerFromAC(integer_x, bestgamma, 0.05);
		}
		// round the last fractional point
		else 
		{	
			frac2int->apply(frac_x, integer_x);
		}
		roundWatch.stop();
		consoleDebug(DebugLevel::Verbose, "roundingTime = {}", roundWatch.getPartial());

		// cycle detection and antistalling actions
		// is it the same of the last one? If yes perturbe
		if (lastIntegerX.size()
			&& areSolutionsEqual(intSubset, integer_x, (*(lastIntegerX.begin())).second, integralityEps)
			&& equal(runningAlpha, (*(lastIntegerX.begin())).first, alphaDist))
		{
			if (!pertCnt)  firstPerturbation = nitr;
			// directly apply restart if PDLP stalls at an integer but not primal feasible solution
			if (!applyPdlpRestart) perturbe(integer_x, ignoreGenerals);
			else restart(integer_x, ignoreGenerals);
		}
		

		if (isInCache(runningAlpha, integer_x, ignoreGenerals))
		{

			if ((!usedOrigFpNoRestart) && (numFracsObj != 1))
			{
				usedOrigFpNoRestart = true;
				integer_x.resize(n, 0);
				frac2int->apply(frac_x, integer_x);

			}

			else if ((!usedOrigFpNoRestart) && (numIntegersObj != 1))
			{
				usedOrigFpNoRestart = true;
			}
			
			else
			{
				usedOrigFpNoRestart = false;
				// do a restart until we are able to insert it in the cache
				for (int rtry = 0; rtry < 10; rtry++)
				{
					if (isInCache(runningAlpha, integer_x, ignoreGenerals)) restart(integer_x, ignoreGenerals);
					else break;
				}
				lastIntegerX.push_front(AlphaVector(runningAlpha, integer_x));
			}

		}
		else
		{
			usedOrigFpNoRestart = false;
			lastIntegerX.push_front(AlphaVector(runningAlpha, integer_x));
		}
		
		// int -> frac
		lpWatch.start();

		double thisAlpha = runningAlpha;
		

		// points to be used in current obj 		
		int currNumPointsInObj = ((usedOrigFpNoRestart)) ? 1 : std::min((int) lastIntegerX.size(), numIntegersObj);

		// if the distance function is not the pure distance one
		// then we might not realize the current integer_x is feasible.
		// so we explictly check for feasibility and, if so,
		// temporarily set the running alpha to zero.
		if (isSolutionFeasible(integer_x, rows))  
		{
			thisAlpha = 0.0;
			currNumPointsInObj = 1;
		}

		// setup distance objective
		int addedVars = 0;
		int addedConstrs = 0;
		std::fill(distObj.begin(), distObj.end(), 0.0);

		// if no score for the integers it is 0 for all
		double scoreInt = 0.0;
		if ((bestObjIntegers) || (lessViolatedIntegers))
		{
			scoreInt = 0.0;
			multipleIntegerX.sort();

		}
		multipleIntegerX.push_front(DistVector(scoreInt, integer_x));
		if (multipleIntegerX.size() > numIntegersObj) multipleIntegerX.resize(numIntegersObj);

		// iterator over previous integer points
		std::list<AlphaVector>::const_iterator itrInt = multipleIntegerX.begin();
		std::list< AlphaVector >::const_iterator itrIntEnd = multipleIntegerX.end();
		
		// int currNumPointsInObj = std::min((int) lastIntegerX.size(), numIntegersObj);

		// scales-weights of the integers in the objective
		std::vector<double> scaleIntVec;

		if (harmonicWeights) scaleIntVec = harmonicVector(currNumPointsInObj);
		else scaleIntVec = exponDecayVector(intScaleFactor,currNumPointsInObj);
		// If aggregateInts=False we use a convex combination of distance functions of previous integers as abjoctive in the projection LP. 
		if (!aggregateInts)
		{
			for(int iter = 0; iter < currNumPointsInObj; iter++)
			{
				for (int j: binaries)
				{
					double distCoef = scaleIntVec[iter];
					
					if (isNull(itrInt->second[j], integralityEps))
					{
						if (expObj) distCoef = scaleIntVec[iter] * 0.5 * exp(-0.5 * frac_x[j]);
						if (logisObj) distCoef = scaleIntVec[iter] * 0.1 * exp(-0.1 * frac_x[j]) / pow(1 + exp(-0.1 * frac_x[j]), 2);
						distObj[j] += distCoef;
					}
					else
					{
						if (expObj) distCoef = scaleIntVec[iter] * 0.5 * exp(-0.5 * ( 1 - frac_x[j] ));
						if (logisObj) distCoef = scaleIntVec[iter] * 0.1 * exp(-0.1 * ( 1 - frac_x[j])) / pow(1 + exp(-0.1 * (1 - frac_x[j])) , 2);
						distObj[j] -= distCoef;
					}

				}

				itrInt++;
			}


			
			
			if (stage > 1)
			{
				itrInt = multipleIntegerX.begin();
				for(int iter = 0 ;iter < currNumPointsInObj ; iter++)
				{


					for (int j: gintegers)
					{		
						double distCoef = scaleIntVec[iter];

						if (equal(itrInt->second[j], lb[j], integralityEps))
						{									
							if (expObj) distCoef = scaleIntVec[iter] * 0.5 * exp(-0.5 * (frac_x[j] - lb[j]));
							if (logisObj) distCoef = scaleIntVec[iter] * 0.1 * exp(-0.1 * (frac_x[j] - lb[j])) / pow(1 + exp(-0.1 * (frac_x[j] - lb[j])), 2);
							distObj[j] += distCoef;

						}
						else if (equal(itrInt->second[j], ub[j], integralityEps))
						{
							if (expObj) distCoef = scaleIntVec[iter] * 0.5 * exp(-0.5 * (ub[j] - frac_x[j]));
							if (logisObj) distCoef = scaleIntVec[iter] * 0.1 * exp(-0.1 * (ub[j] - frac_x[j])) / pow(1 + exp(-0.1 * (ub[j] - frac_x[j])) , 2);
							distObj[j] -= distCoef;							
						}
						
						else
						{	
							bool isSame = isComponentSame(j);

							if (isSame)
							{
								if (iter == 0)
								{
									std::string deltaName = xNames[j] + "_delta_" + std::to_string(iter) ;
									model->addEmptyCol(deltaName, 'C', 0.0, INFBOUND, 0.0);
									int auxIdx = model->ncols() - 1;
									colIndices.push_back(auxIdx);
									if (expObj) distCoef = scaleIntVec[iter] * 0.5 * exp(-0.5 * fabs(frac_x[j] - itrInt->second[j]));
									if (logisObj) distCoef = scaleIntVec[iter] * 0.1 * exp(-0.1 * fabs(frac_x[j] - itrInt->second[j])) / pow(1 + exp(-0.1 * fabs(frac_x[j] - itrInt->second[j])), 2);
									
									distObj.push_back(1.0);
									addedVars++;
									// add constraints
									SparseVector vec;
									vec.push(j, 1.0);
									vec.push(auxIdx, -1.0);
									model->addRow(xNames[j] + "_d1p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'L', itrInt->second[j]);
									addedConstrs++;
									vec.coef()[1] = 1.0;
									model->addRow(xNames[j] + "_d2p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'G', itrInt->second[j]);
									addedConstrs++;

								}
								
							}


							else
							{							
								// add auxiliary variable
								std::string deltaName = xNames[j] + "_delta_" + std::to_string(iter) ;
								model->addEmptyCol(deltaName, 'C', 0.0, INFBOUND, 0.0);
								int auxIdx = model->ncols() - 1;
								colIndices.push_back(auxIdx);
								if (expObj) distCoef = scaleIntVec[iter] * 0.5 * exp(-0.5 * fabs(frac_x[j] - itrInt->second[j]));
								if (logisObj) distCoef = scaleIntVec[iter] * 0.1 * exp(-0.1 * fabs(frac_x[j] - itrInt->second[j])) / pow(1 + exp(-0.1 * fabs(frac_x[j] - itrInt->second[j])), 2);
								
								distObj.push_back(distCoef);
								addedVars++;
								// add constraints
								SparseVector vec;
								vec.push(j, 1.0);
								vec.push(auxIdx, -1.0);
								model->addRow(xNames[j] + "_d1p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'L', itrInt->second[j]);
								addedConstrs++;
								vec.coef()[1] = 1.0;
								model->addRow(xNames[j] + "_d2p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'G', itrInt->second[j]);
								addedConstrs++;
							}						
						}					
					}

					itrInt++;
				}

				consoleDebug(DebugLevel::Verbose, "addedVars={} addedConstrs={}", addedVars, addedConstrs);
				DOMINIQS_ASSERT( distObj.size() == (unsigned int)(n + addedVars) );
				DOMINIQS_ASSERT( distObj.size() == colIndices.size() );
			}
		}
		// If aggregateInts=True we use a convex combination of integers as reference in the projection step. 
		// The resulting point is not an integer, hence it should be rounded.
		else
		{
			
			std::vector<double> aggr_integers;
			aggr_integers.resize(n, 0);
			std::vector<double> new_integer;
			new_integer.resize(n, 0);

			if (isSolutionFeasible(integer_x, rows))  
			{
				new_integer = integer_x;
				thisAlpha = 0.0;
			}
			else
			{
				itrInt = multipleIntegerX.begin();
				for(int iter = 0 ;iter < currNumPointsInObj ; iter++)
				{ 
					accumulate(&aggr_integers[0], &itrInt->second[0], n, scaleIntVec[iter]);
					itrInt++;
				}
				frac2int->apply(aggr_integers, new_integer);
				// the new point may be already in cache! 
				bool sameAsLast = true;
				bool inCache = false;
				std::list< AlphaVector >::const_iterator itr = lastIntegerX.begin();
				std::list< AlphaVector >::const_iterator end = lastIntegerX.end();
				if (stage == 1)
				{
					if (!areSolutionsEqual(binaries, itr->second, new_integer, integralityEps)) 
					{
						sameAsLast = false;
					}
				}
				else
				{
					if (!areSolutionsEqual(integers, itr->second, new_integer, integralityEps)) 
					{
						sameAsLast = false;
					}
				}
				itr ++; // to ignore first integer
				if (!sameAsLast)
				{
					if (stage == 1)
					{
						while ((itr != end) && !inCache)
						{
							if ((fabs(thisAlpha - itr->first) < alphaDist) && areSolutionsEqual(binaries, itr->second, new_integer, integralityEps)) inCache = true;
							++itr;
						}
					}
					else
					{
						while ((itr != end) && !inCache)
						{
							if ((fabs(thisAlpha - itr->first) < alphaDist) && areSolutionsEqual(integers, itr->second, new_integer, integralityEps)) inCache = true;
							++itr;
						}
					}
					if (inCache) 
					{	
						new_integer = integer_x;
					}
					else if (isSolutionFeasible(new_integer, rows))  
					{
						thisAlpha = 0.0;
					}

				}

			}
			

			for (int j: binaries)
			{
				if (isNull(new_integer[j], integralityEps))
				{
					double distCoef = 1.0;
					distObj[j] = distCoef;
				}
				else
				{
					double distCoef = -1.0;
					distObj[j] = distCoef;
				}
			}


			if (stage > 1)
			{
				for (int j: gintegers)
				{
					if (equal(new_integer[j], lb[j], integralityEps))
					{
						distObj[j] = 1.0;
					}
					else if (equal(new_integer[j], ub[j], integralityEps))
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
						model->addRow(xNames[j] + "_d1", vec.idx(), vec.coef(), 2, 'L', new_integer[j]);
						addedConstrs++;
						vec.coef()[1] = 1.0;
						model->addRow(xNames[j] + "_d2", vec.idx(), vec.coef(), 2, 'G', new_integer[j]);
						addedConstrs++;
					}
				}

				consoleDebug(DebugLevel::Verbose, "addedVars={} addedConstrs={}", addedVars, addedConstrs);
				DOMINIQS_ASSERT( distObj.size() == (unsigned int)(n + addedVars) );
				DOMINIQS_ASSERT( distObj.size() == colIndices.size() );
			}
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

			double distScale = 0.0;
			double objScale = 0.0;

			if (newScaleC)
			{
				// scale by the objective value of the current LP solution
				objScale = dotProduct(&obj[0], &frac_x[0], n) + objOffset;
				objScale = fabs(objScale);
			}

			// else scale as normal by the l2 norm of the objective vector
			else objScale = objNorm;

			if (newScaleDelta)
			{	
				itrInt = multipleIntegerX.begin();
				for (int iter = 0; iter < currNumPointsInObj; iter++)
				{
					distScale += scaleIntVec[iter] * solutionsDistance(intSubset, frac_x, itrInt->second);
					itrInt++;
				}
			}

			else
			{	
				// distScale = sqrt(intSubset.size());
				distScale = sqrt(dotProduct(&distObj[0], &distObj[0], (int) distObj.size()));
			}

			// scale distance objective by (1-thisAlpha)
			for (int j = 0; j < (n + addedVars); j++)
			{
				distObj[j] *= (1.0 - thisAlpha);
			}	

			if (isNull(objScale))
			{
				objScale = 1.0;
			}

			if (isNull(distScale))
			{
				distScale = 1.0;
			}
			// add objective with proper weight
			double weight = (thisAlpha * distScale) / objScale;
			accumulate(&distObj[0], &obj[0], n, weight);
		
		}
		// set objective
		model->objcoefs(colIndices.size(), &colIndices[0], &distObj[0]);

		// solve LP
		if (lpIterLimit > 0)  model->intParam(IntParam::IterLimit, lpIterLimit);
		model->dblParam(DblParam::TimeLimit, timeLeft);
		timeModel += model->lpopt(reOptMethod, false, false);
		lpWatch.stop();

		// if no solution is found because of timeLeft return false
		primalFeas = model->isPrimalFeas();
		if (primalFeas)
		{
			// get solution
			model->sol(&frac_x[0], 0, n-1);
			// try to tighten the tolerance for pdlp
			while (pdlpTol > 1e-6 && !isSolutionFeasible(frac_x, rows) && isSolutionInteger(intSubset, frac_x, integralityEps))
			{
				pdlpTol *= 0.1;
				/* decrease the tolerance for pdlp */
				timeModel += model->lpopt(reOptMethod, true, false);
				if(model->isPrimalFeas())
					model->sol(&frac_x[0], 0, n-1);
				else
					break;
			}
		}
		if (!primalFeas) 
		{
			// if lpopt was interrupted because of the timelimit postsolve the model manually and exit the current stage
			model->postsolve();
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
			return false;
		}
		consoleDebug(DebugLevel::VeryVerbose, "Iteration {}: time={} pFeas={} lpiter={} pdlp={}",
				lpWatch.getPartial(), primalFeas, model->intAttr(IntAttr::SimplexIterations), model->intAttr(IntAttr::PDLPIterations));
		double projObj = model->objval();

		dominiqs::StopWatch chronoDelRows;
		chronoDelRows.start();

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
		
		chronoDelRows.stop();
		timeModel += chronoDelRows.getElapsed();

		colIndices.resize(n);
		distObj.resize(n);
		DOMINIQS_ASSERT( model->ncols() == n );

		// get some statistics
		double origObj = dotProduct(&obj[0], &frac_x[0], n) + objOffset;
		double dist = solutionsDistance(intSubset, frac_x, integer_x);
		int numFrac = solutionNumFractional(intSubset, frac_x, integralityEps);
		
		if (lessViolatedIntegers)
		{
			std::list<AlphaVector>::iterator itrIntegers = multipleIntegerX.begin();
			std::vector<double> integer_afterProj;
			integer_afterProj = frac_x;
			for (int j: integers) integer_afterProj[j] = integer_x[j];
			double intViolation = 0.0;
			for (auto c: rows) 
			{	
				intViolation += std::max(0.0, c->violation(&integer_afterProj[0]));
			}
			itrIntegers->first = intViolation;
			if (bestObjIntegers) 
			{
				double objOfInt = dotProduct(&obj[0], &integer_afterProj[0], n) + objOffset;
				if ( objOfInt > 0 ) itrIntegers->first = intViolation * objOfInt;
				else if ( objOfInt == 0 ) itrIntegers->first = intViolation;
				else itrIntegers->first = (1.0 / (1.0 + intViolation)) * objOfInt;
			}
			else itrIntegers->first = intViolation;

		}
		if (lessDistanceIntegers)
		{
			std::list<AlphaVector>::iterator itrIntegers = multipleIntegerX.begin();
			if (bestObjIntegers) itrIntegers->first = dist*(dotProduct(&obj[0], &integer_x[0], n)+ objOffset);
			else itrIntegers->first = dist;
		}

		if (multirensStage3)
		{
			for (int j: integers)
			{	
				if (equal(frac_x[j], floor(frac_x[j]), integralityEps))
				{
					newlb[j] = std::min(newlb[j], floor(frac_x[j]));
					newub[j] = std::max(newub[j], floor(frac_x[j]));
				}
				else 
				{
					newlb[j] = std::min(newlb[j], floor(frac_x[j]));
					newub[j] = std::max(newub[j], ceil(frac_x[j]));
				}
			}
		}
		// get best integer cache for stage 3 based on score
		if (newStage3)
		{	
			double scoreIntS3 = 0.0;
			// if ((stage3bestObjIntegers) && (stage3lessViolatedIntegers))
			// {
			// 	// how much does the integer violate the constraints?
			// 	for (auto c: rows) 
			// 	{	
			// 		scoreIntS3 += std::max(0.0, c->violation(&integer_x[0]));
			// 	}
			// 	scoreIntS3 *= (dotProduct(&obj[0], &integer_x[0], n) + objOffset);

			// }
			if (stage3lessViolatedIntegers)
			{
				std::vector<double> integer_afterProj;
				integer_afterProj = frac_x;
				for (int j: integers) integer_afterProj[j] = integer_x[j];
				// how much does the integer violate the constraints?
				for (auto c: rows) 
				{	
					scoreIntS3 += std::max(0.0, c->violation(&integer_afterProj[0]));
				}
				if (stage3bestObjIntegers)
				{
					double objOfInt = dotProduct(&obj[0], &integer_afterProj[0], n) + objOffset;
					if ( objOfInt > 0 ) scoreIntS3 *= objOfInt;
					else if ( objOfInt < 0 ) scoreIntS3 = (1.0 / (1.0 + scoreIntS3)) * objOfInt;
				}
			}
			else
			{
				scoreIntS3 = dist;
				if (stage3bestObjIntegers)
				{

					std::vector<double> integer_afterProj;
					integer_afterProj = frac_x;
					for (int j: integers) integer_afterProj[j] = integer_x[j];
					double currentObjValue = getSolutionValue(integer_afterProj);
					scoreIntS3 *= (1 + abs(currentObjValue - dualBound));
				}

			}
			std::list<DistVector> *p;

			if (stage == 1) p = &closestIntegerXs1;
			else p = &closestIntegerXs2;

			std::list< DistVector >::const_iterator itr_cl = p->begin();
			std::list< DistVector >::const_iterator end_cl = p->end();
			
			bool pointAdded = false;
			while ((itr_cl != end_cl))
			{	

				if (lessEqualThan(scoreIntS3, itr_cl->first))
				{	
					pointAdded = true;
					p->insert(itr_cl,DistVector(scoreIntS3, integer_x));
					break;
				}
				itr_cl++;
			}
			if ((p->size() < stage3IntegersObj) && (!pointAdded)) 
			{
				p->push_back(DistVector(scoreIntS3, integer_x));
			}	
			if (p->size() > stage3IntegersObj) p->resize(stage3IntegersObj);
		}



		// save integer_x as best point if distance decreased and count number of itr without 10% improvement
		if (dist < closestDist)
		{
			// if 10% improvenent set iterationsNoImpr = 0
			if( dist / closestDist < 0.9) iterationsNoImpr = 0;
			closestDist = dist;
			closestPoint = integer_x;
			closestFrac = frac_x;
		}
		
		else iterationsNoImpr ++;


		//break loop if too many iterations without 10% improvement
		if( iterationsNoImpr > maxItrNoImpr )
		{
			std::cout << "Too many iterations without 10% improvement" << std::endl;
			break;
		}
		
		// ToDo this check is done twice per iteration. Remove one.
		lpfeasible = isSolutionFeasible(frac_x, rows);
		if (lpfeasible) 
		{	
			// For PDLP this is just the best known LP solution value and not the dual bound
			if (dualBound == -static_cast<int>(model->objSense()) * INFBOUND) dualBound = origObj;
			if (model->objSense() == ObjSense::MIN) dualBound = std::min(dualBound, origObj);
			else dualBound = std::max(dualBound, origObj);
		}
		// display log
		if (display.needPrint(nitr))
		{
			std::string reason;
			model->terminationReason(reason);
			int simplexIt = model->intAttr(IntAttr::SimplexIterations);
			int barrierIt = model->intAttr(IntAttr::BarrierIterations);
			int pdlpIt = model->intAttr(IntAttr::PDLPIterations);
			display.set("stage", stage);
			display.set("iter", nitr);
			display.set("alpha", runningAlpha);
			display.set("origObj", origObj);
			display.set("time", chrono.getElapsed());
			display.set("dist", dist);
			display.set("#frac", numFrac);
			display.set("projObj", projObj);
			display.set("lpiter", std::max(simplexIt, barrierIt));
			display.set("PDLP status", reason);
			display.set("PDLP feas", lpfeasible);
			display.set("PDLP iter", pdlpIt);
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
	// double s3TimeLimit = std::max(std::min(remainingTime, elapsedTime), 1.0);
	double s3TimeLimit = std::max(remainingTime, 1.0);
	if (lessThan(remainingTime, 0.1)) return false;

	model->switchToMIP();
	int n = model->ncols();
	bool found = false;
	std::vector<std::string> xNames;
	model->colNames(xNames);

	// restore type information
	std::vector<char> ctype(n, 'C');
	for (int j: binaries) ctype[j] = 'B';
	for (int j: gintegers) ctype[j] = 'I';
	for (int j = 0; j < n; j++)  model->ctype(j, ctype[j]);

	int addedVars = 0;
	int addedConstrs = 0;
	std::vector<double> distObj(n, 0.0);
	std::vector<int> colIndices(n);
	std::iota(colIndices.begin(), colIndices.end(), 0);

	if ( rensClosestDistStage3 )
	{
		for (int j: integers)
		{	
			if (equal(closestFrac[j], floor(closestFrac[j]), integralityEps))
			{
				newlb[j] = std::min(newlb[j], floor(closestFrac[j]));
				newub[j] = std::max(newub[j], floor(closestFrac[j]));
			}
			else 
			{
				newlb[j] = std::min(newlb[j], floor(closestFrac[j]));
				newub[j] = std::max(newub[j], ceil(closestFrac[j]));
			}
		}
	}
	if ( rensStage3 || multirensStage3 || rensClosestDistStage3 )
	{
		int changedBoundsBins = 0;
		int changedBoundsInts = 0;
		int FixedInts = 0;

		for (int j: binaries)
		{	
			if ( (newlb[j] == newub[j]) )
			{	
				model->lb(j,newlb[j]);
				model->ub(j,newub[j]);
				changedBoundsBins +=1;
			}	
		}

		for (int j: gintegers)
		{	
			if ( (newub[j] - newlb[j]) < (ub[j] - lb[j]) )
			{		
					model->lb(j,newlb[j]);
					model->ub(j,newub[j]);
					changedBoundsInts +=1;
					if ( newub[j] == newlb[j] ) FixedInts += 1;
			}
		}
		std::cout << "Out of " << binaries.size() << " binary variables " << changedBoundsBins << " changed bounds" << std::endl;
		std::cout << "Out of " << gintegers.size() << " integer variables " << changedBoundsInts << " changed bounds" << std::endl;
		std::cout << "Out of " << gintegers.size() << " integer variables " << FixedInts << " were fixed" << std::endl;


	}


	if (rensClosestDistStage3)
	{
		consoleLog("rens closest point stage3");
	}
	else if (normalMIPStage3)
	{
		consoleLog("solve normal MIP stage3");
	}
	else if (rensStage3)
	{
		consoleLog("RENS stage3");

	}
	else if (multirensStage3)
	{
		consoleLog("Multi-RENS stage3");
	}

	else if (!newStage3)
	{
		// load best point and generate objective
		integer_x = closestPoint;

		consoleLog("Starting stage3 from point with distance={} [timeLimit={}]", closestDist, s3TimeLimit);

		for (int j: binaries) distObj[j] = (isNull(integer_x[j], integralityEps) ? 1.0 : -1.0);
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
	}

	if (newStage3)
	{	

		std::list<DistVector> *p;

		if (closestIntegerXs2.size() >= 1) p = &closestIntegerXs2;
		else p = &closestIntegerXs1;


		std::list<DistVector>::const_iterator itrInt = p->begin();
		std::list< DistVector >::const_iterator itrInt_end = p->end();

		// scales-weights of the integers in the objective
		std::vector<double> scaleIntVec;


		int currNumPointsInObj = (int) p->size();

		consoleLog("Starting stage3 from {} closest points to LP [timeLimit={}] ", currNumPointsInObj, s3TimeLimit);

		if (stage3harmonicWeights) scaleIntVec = harmonicVector(currNumPointsInObj);
		else scaleIntVec = exponDecayVector(stage3ScaleFactor,currNumPointsInObj);

		for(int iter = 0; iter < currNumPointsInObj; iter++)
		{
			for (int j: binaries)
			{
				if (isNull(itrInt->second[j], integralityEps))
				{
					distObj[j] += scaleIntVec[iter];
				}
				else
				{
					distObj[j] += -scaleIntVec[iter];
				}

			}

			itrInt++;
		}

		itrInt = p->begin();
		for(int iter = 0 ;iter < currNumPointsInObj ; iter++)
		{
			for (int j: gintegers)
			{

				// TODO: penalty objective for general integers?
				if (equal(itrInt->second[j], lb[j], integralityEps))
				{
					distObj[j] += scaleIntVec[iter];
				}
				else if (equal(itrInt->second[j], ub[j], integralityEps))
				{
					distObj[j] += -scaleIntVec[iter];
				}
				
				else
				{
					bool isSame = isComponentSameS3(j);
					if (isSame)
					{	
						if (iter == 0)
						{
							// add auxiliary variable
							std::string deltaName = xNames[j] + "_deltaS3_" + std::to_string(iter) ;
							model->addEmptyCol(deltaName, 'C', 0.0, INFBOUND, 0.0);
							int auxIdx = model->ncols() - 1;
							colIndices.push_back(auxIdx);
							distObj.push_back(1.0);
							addedVars++;
							// add constraints
							SparseVector vec;
							vec.push(j, 1.0);
							vec.push(auxIdx, -1.0);
							model->addRow(xNames[j] + "_d1p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'L', itrInt->second[j]);
							addedConstrs++;
							vec.coef()[1] = 1.0;
							model->addRow(xNames[j] + "_d2p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'G', itrInt->second[j]);
							addedConstrs++;
						}
					}
					else
					{
						// add auxiliary variable
						std::string deltaName = xNames[j] + "_deltaS3_" + std::to_string(iter) ;
						model->addEmptyCol(deltaName, 'C', 0.0, INFBOUND, 0.0);
						int auxIdx = model->ncols() - 1;
						colIndices.push_back(auxIdx);
						distObj.push_back(scaleIntVec[iter]);
						addedVars++;
						// add constraints
						SparseVector vec;
						vec.push(j, 1.0);
						vec.push(auxIdx, -1.0);
						model->addRow(xNames[j] + "_d1p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'L', itrInt->second[j]);
						addedConstrs++;
						vec.coef()[1] = 1.0;
						model->addRow(xNames[j] + "_d2p" + std::to_string(iter) , vec.idx(), vec.coef(), 2, 'G', itrInt->second[j]);
						addedConstrs++;
					}
					
				}
			}
			itrInt++;
		}
	}



	LOG_ITEM("addedVars", addedVars);
	LOG_ITEM("addedConstrs", addedConstrs);

	consoleDebug(DebugLevel::Normal, "addedVars={} addedConstrs={}", addedVars, addedConstrs);
	DOMINIQS_ASSERT( distObj.size() == (unsigned int)(n + addedVars) );
	DOMINIQS_ASSERT( distObj.size() == colIndices.size() );

	if (normalMIPStage3) model->objcoefs(colIndices.size(), &colIndices[0], &obj[0]);
	else if (newStage3 && (multirensStage3 || rensClosestDistStage3)) model->objcoefs(colIndices.size(), &colIndices[0], &distObj[0]);
	else if ( rensStage3 || multirensStage3 || rensClosestDistStage3) model->objcoefs(colIndices.size(), &colIndices[0], &obj[0]);
	else model->objcoefs(colIndices.size(), &colIndices[0], &distObj[0]);

	// model->writeModel(std::to_string(binaries.size())+".lp", "lp");
	model->logging(true);
	model->intParam(IntParam::SolutionLimit, 1);
	model->intParam(IntParam::NodeLimit, 500);
	model->dblParam(DblParam::TimeLimit, s3TimeLimit);
	elapsedTime = chrono.getElapsed();
	model->mipopt();
	stage3Time = chrono.getElapsed() - elapsedTime;
	primalFeas = model->isPrimalFeas();

	model->logging(false);

	// postsolve the problem since it may still be presolved!
	model->postsolve();

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
		// ToDo check value of nrows() in all stages of feaspump
		for (int i = 0; i < rows.size(); i++) {
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
void FeasibilityPump::aggregateFracs( std::vector<double>& aggr_frac_x, std::vector<double>& scaleVector)
{
	// iterator over fractional point
	std::list< NumberVector >::const_iterator itr = lastFracX.begin();
	std::list< NumberVector >::const_iterator end = lastFracX.end();
	// number of points used to compute aggr_frac_x
	int numPoints = scaleVector.size();

	for (int i = 0; i < aggr_frac_x.size(); i++)
	{
		int curr = 0;
		while (curr < numPoints)
		{
			aggr_frac_x[i] += (itr->second[i])*scaleVector[curr];
			curr ++;
			itr ++;
		}
		itr = lastFracX.begin();
	}
}

void FeasibilityPump::computeAC(std::vector<double>& ac_x)
{
	int n = model->ncols();
	std::vector<double> OrigObj(n, 0);
	std::vector<double> EmptyObj(n, 0);
	std::vector<int> colIndices(n);
	std::iota(colIndices.begin(), colIndices.end(), 0);
	// get original objective
	model->objcoefs(&OrigObj[0]);
	// remove original objective
	model->objcoefs(colIndices.size(), &colIndices[0], &EmptyObj[0]);
	// compute AC
	model->lpopt('A', false, false);
	model->sol(&ac_x[0]);
	DOMINIQS_ASSERT(model->isPrimalFeas());
	// add original objective again
	model->objcoefs(colIndices.size(), &colIndices[0], &OrigObj[0]);
}

void FeasibilityPump::integerFromAC(std::vector<double>& x, double& bestgamma, double step)
{		
		// currentFrac obtained by (1.0 - gamma) * frac_x + gamma * ac_x
		int n = model->ncols();
		std::vector<double> currentFrac(n, 0);
		std::vector<double> currentInt(n, 0);

		double smallestViolation = INFBOUND;
		double intViolation = 0.0;

		for (double gamma = 0.0; gamma < (1+integralityEps); gamma = gamma + step)
	    {	
			for (int i = 0; i < n; i++) currentFrac[i] = (1.0 - gamma) * frac_x[i] + gamma * ac_x[i];
			frac2int->apply(currentFrac, currentInt);
			if (isSolutionFeasible(currentInt, rows)) 
			{	
				for (int i = 0; i < n; i++) x[i] = currentInt[i]; 
				bestgamma = gamma;
				std::cout << "AC-FP found a solution with gamma: "<< gamma << std::endl;
				return;
			}

			else
			{	
				// compute total constraint violation for the current integer point
				intViolation = 0.0;
				for (auto c: rows) 
				{	
					intViolation = std::max(intViolation, c->violation(&currentInt[0]));
				}
				// choose integer with the smallest total constraint violation
				if (intViolation < smallestViolation)
				{
					smallestViolation = intViolation;
					bestgamma = gamma;
					for (int i = 0; i < n; i++) x[i] = currentInt[i]; 
				}
			}
		}
		std::cout << bestgamma << std::endl;

}

bool FeasibilityPump::isComponentSame(int idx)
{
	bool isSame = true;
	// iterator over reference integer points
	std::list<AlphaVector>::const_iterator itrInt = multipleIntegerX.begin();
	std::list< AlphaVector >::const_iterator itrIntEnd = multipleIntegerX.end();
	double value =  itrInt->second[idx];
	while ((itrInt != itrIntEnd) && isSame)
	{
		if ( itrInt->second[idx] != value ) isSame = false;
		++itrInt;

	}
	return isSame;
}

bool FeasibilityPump::isComponentSameS3(int idx)
{
	bool isSame = true;
	std::list<DistVector> *p;

	if (closestIntegerXs2.size() >= 1) p = &closestIntegerXs2;
	else p = &closestIntegerXs1;


	std::list<DistVector>::const_iterator itrInt = p->begin();
	std::list< DistVector >::const_iterator itrIntEnd = p->end();
	double value =  itrInt->second[idx];
	while ((itrInt != itrIntEnd) && isSame)
	{
		if ( itrInt->second[idx] != value ) isSame = false;
		++itrInt;

	}
	return isSame;
}

} // namespace dominiqs
