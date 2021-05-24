/**
 * @file transformers.cpp
 * @brief Solution Transformers (a.k.a. rounders) Source
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 */

#include <map>
#include <list>
#include <iterator>
#include <set>
#include <signal.h>
#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>

#include <utils/floats.h>
#include <utils/fileconfig.h>
#include <utils/consolelog.h>

#include "feaspump/transformers.h"

using namespace dominiqs;


// macro type savers
#define READ_FROM_CONFIG( what, defValue ) what = gConfig().get("fp."#what, defValue)
#define LOG_ITEM(name, value) consoleLog("{} = {}", name, value)
#define LOG_CONFIG( what ) LOG_ITEM("fp."#what, what)

// Rounders

static bool DEF_RANDOMIZED_ROUNDING = true;
static bool DEF_LOG_DETAILS = false;
static uint64_t DEF_SEED = 0;

SimpleRounding::SimpleRounding() : randomizedRounding(DEF_RANDOMIZED_ROUNDING), logDetails(DEF_LOG_DETAILS)
{
}

void SimpleRounding::readConfig()
{
	READ_FROM_CONFIG( randomizedRounding, DEF_RANDOMIZED_ROUNDING );
	READ_FROM_CONFIG( logDetails, DEF_LOG_DETAILS );
	consoleInfo("[config rounder]");
	LOG_CONFIG( randomizedRounding );
	LOG_CONFIG( logDetails );
	uint64_t seed = gConfig().get<uint64_t>("seed", DEF_SEED);
	roundGen.setSeed(seed);
	roundGen.warmUp();
}

void SimpleRounding::init(MIPModelPtr model, bool ignoreGeneralInt)
{
	binaries.clear();
	gintegers.clear();
	integers.clear();
	int ncols = model->ncols();
	std::vector<double> xLb(ncols);
	std::vector<double> xUb(ncols);
	std::vector<char> xType(ncols);
	model->lbs(&xLb[0]);
	model->ubs(&xUb[0]);
	model->ctypes(&xType[0]);
	for (int j = 0; j < ncols; j++)
	{
		if (different(xLb[j], xUb[j]))
		{
			if (xType[j] == 'B') binaries.push_back(j);
			if (xType[j] == 'I') gintegers.push_back(j);
		}
	}
	ignoreGeneralIntegers(ignoreGeneralInt);
}

void SimpleRounding::ignoreGeneralIntegers(bool flag)
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

void SimpleRounding::apply(const std::vector<double>& in, std::vector<double>& out)
{
	copy(in.begin(), in.end(), out.begin());
	int rDn = 0;
	int rUp = 0;
	double t = getRoundingThreshold(randomizedRounding, roundGen);
	for (int j: integers)
	{
		doRound(in[j], out[j], t);
		if (lessThan(out[j], in[j])) rDn++;
		if (greaterThan(out[j], in[j])) rUp++;
	}
	consoleDebug(DebugLevel::VeryVerbose, "rounding: thr={} #down={} #up={}", t, rDn, rUp);
}

PropagatorRounding::PropagatorRounding() {}

void PropagatorRounding::readConfig()
{
	SimpleRounding::readConfig();
	std::string rankerName = gConfig().get("fp.ranker", std::string("FRAC"));
	filterConstraints = gConfig().get("fp.filterConstraints", true);
	consoleInfo("[config rounder]");
	LOG_ITEM("fp.ranker", rankerName);
	LOG_ITEM("fp.filterConstraints", filterConstraints);
	ranker = RankerPtr(RankerFactory::getInstance().create(rankerName));
	ranker->readConfig();
}

void PropagatorRounding::init(MIPModelPtr model, bool ignoreGeneralInt)
{
	SimpleRounding::init(model, ignoreGeneralInt);
	domain = std::make_shared<Domain>();
	// add vars to domain
	int ncols = model->ncols();
	std::vector<double> xLb(ncols);
	std::vector<double> xUb(ncols);
	std::vector<char> xType(ncols);
	std::vector<std::string> xNames;
	model->lbs(&xLb[0]);
	model->ubs(&xUb[0]);
	model->ctypes(&xType[0]);
	model->colNames(xNames);
	for (int j = 0; j < ncols; j++) domain->pushVar(xNames[j], xType[j], xLb[j], xUb[j]);
	// connect domain to engine and ranker
	prop.setDomain(domain);
	ranker->init(domain, ignoreGeneralInt);
	// generate propagators with analyzers
	std::list<std::string> fNames;
	PropagatorFactories::getInstance().getIDs(std::back_insert_iterator< std::list<std::string> >(fNames));
	for (std::string name: fNames)
	{
		PropagatorFactoryPtr fact(PropagatorFactories::getInstance().create(name));
		factories[fact->getPriority()] = fact;
	}

	int filteredOut = 0;
	for (int i = 0; i < model->nrows(); i++)
	{
		std::map<int, PropagatorFactoryPtr>::iterator itr = factories.begin();
		std::map<int, PropagatorFactoryPtr>::iterator end = factories.end();
		ConstraintPtr c = std::make_shared<Constraint>();
		model->row(i, c->row, c->sense, c->rhs, c->range);
		// ignore nonbinding constraints
		if (c->sense == 'N')  continue;
		// constraint filter
		if (filterConstraints)
		{
			const int* idx = c->row.idx();
			const double* coef = c->row.coef();
			unsigned int size = c->row.size();
			bool allCont = true;
			double largest = std::numeric_limits<double>::min();
			double smallest = std::numeric_limits<double>::max();
			for (unsigned int k = 0; k < size; k++)
			{
				if (!domain->isVarFixed(idx[k]) && (domain->varType(idx[k]) != 'C'))
				{
					allCont = false;
					break;
				}
				double tmp = fabs(coef[k]);
				largest = std::max(largest, tmp);
				smallest = std::min(smallest, tmp);
			}
			double dynamism = (largest / smallest);
			if ((allCont && greaterThan(dynamism, 10.0)) || greaterThan(dynamism, 1000.0))
			{
				filteredOut++;
				continue;
			}
		}
		// try analyzers
		while (itr != end)
		{
			PropagatorPtr p = itr->second->analyze(*(domain.get()), c.get());
			if (p)
			{
				prop.pushPropagator(p);
				break;
			}
			itr++;
		}
	}
	// log prop stats
	consoleInfo("[propagator stats]");
	for (const auto& kv: factories)
	{
		consoleLog("{}: {}", kv.second->getName(), kv.second->created());
	}
	consoleLog("#filtered out: {}\n", filteredOut);

	// no initial propagation
	// prop.propagate();
	state = prop.getStateMgr();
	state->dump();
}

void PropagatorRounding::ignoreGeneralIntegers(bool flag)
{
	SimpleRounding::ignoreGeneralIntegers(flag);
	ranker->ignoreGeneralIntegers(flag);
}

void PropagatorRounding::apply(const std::vector<double>& in, std::vector<double>& out)
{
	copy(in.begin(), in.end(), out.begin());
	state->restore();
	double t = getRoundingThreshold(randomizedRounding, roundGen);
	ranker->setCurrentState(in);
	// main loop
	int next;
	while ((next = ranker->next()) >= 0)
	{
 		// standard rounding
		if (domain->varType(next) == 'B') doRound(in[next], out[next], t);
		else
		{
			// general integer variable: here we take into account also the tightened domain
			// in order to have a more clever rounding!
			if (lessEqualThan(in[next], domain->varLb(next))) out[next] = domain->varLb(next);
			else if (greaterEqualThan(in[next], domain->varUb(next))) out[next] = domain->varUb(next);
			else doRound(in[next], out[next], t);
		}
		// propagate
		prop.propagate(next, out[next]);
		DOMINIQS_ASSERT( domain->isVarFixed(next) );
		// update with fixings
		for (int j: prop.getLastFixed()) out[j] = domain->varLb(j);
	}
}

void PropagatorRounding::clear()
{
	// clear
	// delete state;
	prop.clear();
	factories.clear();
}

// auto registration
// please register your class here with an appropriate name

class TR_FACTORY_RECORDER
{
public:
	TR_FACTORY_RECORDER()
	{
		std::cout << "Registering SolutionTransformers...";
		TransformersFactory::getInstance().registerClass<SimpleRounding>("std");
		TransformersFactory::getInstance().registerClass<PropagatorRounding>("propround");
		std::cout << "done" << std::endl;
	}
};

TR_FACTORY_RECORDER my_tr_factory_recorder;
