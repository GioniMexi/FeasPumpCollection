/**
 * @file main.cpp
 * @brief Main App
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 */

#include <iostream>
#include <iomanip>
#include <fstream>

#include <utils/args_parser.h>
#include <utils/fileconfig.h>
#include <utils/path.h>
#include <utils/floats.h>
#include <utils/timer.h>
#include <utils/randgen.h>
#include <utils/consolelog.h>
#include <utils/path.h>

#include "feaspump/feaspump.h"
#include "feaspump/version.h"
#ifdef HAS_CPLEX
#include "feaspump/cpxmodel.h"
#endif
#ifdef HAS_XPRESS
#include "feaspump/xprsmodel.h"
#endif
#include <fmt/format.h>


// macro type savers
#define LOG_ITEM(name, value) consoleLog("{} = {}", name, value)

using namespace dominiqs;

static const uint64_t DEF_SEED = 0;

int main (int argc, char const *argv[])
{
	// config/options
	ArgsParser args;
	args.parse(argc, argv);
	if (args.input.size() < 1)
	{
		consoleError("usage: feaspump prob_file");
		return -1;
	}
	mergeConfig(args, gConfig());
	std::string runName = gConfig().get("runName", std::string("default"));
	std::string testset = gConfig().get("testset", std::string("unknown"));
	std::string solver = gConfig().get("solver", std::string("cpx"));
	bool mipPresolve = gConfig().get("mipPresolve", true);
	int numThreads = gConfig().get("numThreads", 0);
	bool printSol = gConfig().get("printSol", false);
	double timeLimit = gConfig().get("fp.timeLimit", 1e+75);
	std::string probName = getProbName(Path(args.input[0]).getBasename());
	// logger
	consoleInfo("Timestamp: {}", currentDateTime());
	consoleInfo("[config]");
	LOG_ITEM("probName", probName);
	LOG_ITEM("testset", testset);
	LOG_ITEM("solver", solver);
	LOG_ITEM("runName", runName);
	LOG_ITEM("presolve", mipPresolve);
	LOG_ITEM("numThreads", numThreads);
	LOG_ITEM("gitHash", FP_GIT_HASH);
	LOG_ITEM("fpVersion", FP_VERSION);
	LOG_ITEM("printSol", printSol);
	// seed
	uint64_t seed = gConfig().get<uint64_t>("seed", DEF_SEED);
	LOG_ITEM("seed", seed);
	seed = generateSeed(seed);
	gConfig().set<uint64_t>("seed", seed);

	MIPModelPtr model;
#ifdef HAS_CPLEX
	if (solver == "cpx")  model = MIPModelPtr(new CPXModel());
#else
	if (solver == "cpx")  throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif
#ifdef HAS_XPRESS
	if (solver == "xprs")  model = MIPModelPtr(new XPRSModel());
#else
	if (solver == "xprs")  throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif

	if (!model)  throw std::runtime_error("No solver available for FP");

	DOMINIQS_ASSERT(model);
	double integralityEps = model->dblParam(DblParam::IntegralityTolerance);
	gConfig().set("fp.integralityEps", integralityEps);
	model->logging(false);
	try
	{
		model->readModel(args.input[0]);
		consoleLog("originalProblem: #rows={} #cols={} #nnz={}",
					model->nrows(), model->ncols(), model->nnz());

		// presolve
		MIPModelPtr premodel;
		bool hasPresolve = false;

		if (mipPresolve)
		{
			model->dblParam(DblParam::TimeLimit, timeLimit);
			model->presolve();
			premodel = model->presolvedModel();
			if (!premodel)
			{
				// presolve made no reduction: just clone the original model
				premodel = model->clone();
			}
			else
			{
				hasPresolve = true;
#ifdef DUMPLP
				premodel->writeModel("presolved.sav.gz");
#endif //< DUMPLP
				consoleLog("presolvedProblem: #rows={} #cols={} #nnz={}",
						premodel->nrows(), premodel->ncols(), premodel->nnz());
			}
		}
		else
		{
			// presolve disabled: just clone the original model
			premodel = model->clone();
		}
		DOMINIQS_ASSERT( premodel );

		// feaspump
		FeasibilityPump solver;
		solver.readConfig();
		gStopWatch().start();
		solver.init(premodel);
		solver.pump();
		std::vector<double> x;
		if (solver.foundSolution())
		{
			// uncrush solution
			std::vector<double> preX;
			solver.getSolution(preX);
			if (hasPresolve)
			{
				DOMINIQS_ASSERT( (int)preX.size() == premodel->ncols() );
				x = model->postsolveSolution(preX);
				model->postsolve();
			}
			else x = preX;

			// compute objective in original space
			int n = model->ncols();
			std::vector<double> obj(n);
			model->objcoefs(&obj[0]);
			double objValue = model->objOffset();
			objValue += dotProduct(&obj[0], &x[0], n);

			if (printSol)
			{
				// print solution
				std::cout << "Solution:" << std::endl;
				std::cout << "=obj= " << std::setprecision(15) << objValue << std::endl;
				std::vector<std::string> xNames;
				model->colNames(xNames);
				DOMINIQS_ASSERT( xNames.size() == x.size() );
				for (unsigned int i = 0; i < x.size(); i++)
				{
					if (isNotNull(x[i], integralityEps))  std::cout << xNames[i] << " " << std::setprecision(15) << x[i] << std::endl;
				}
			}

			// check solution for feasibility
			int m = model->nrows();
			std::vector<std::string> rNames;
			model->rowNames(rNames);
			for (int i = 0; i < m; i++)
			{
				ConstraintPtr c = std::make_shared<Constraint>();
				model->row(i, c->row, c->sense, c->rhs, c->range);
				if (c->sense == 'N')  continue;
				if (!c->satisfiedBy(&x[0]))  throw std::runtime_error(fmt::format("Constraint {} violated by {}", rNames[i], c->violation(&x[0])));
			}

		}
		solver.reset();
		gStopWatch().stop();
	}
	catch(std::exception& e)
	{
		consoleError(e.what());
	}
	return 0;
}
