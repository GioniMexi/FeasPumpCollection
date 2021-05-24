/**
 * @file prop_engine.h
 * @brief Constraint Propagation Engine
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008-2012
 */

#ifndef PROP_ENGINE_H
#define PROP_ENGINE_H

#include <vector>
#include <list>
#include <map>
#include <deque>
#include <string>
#include <iosfwd>

#include <utils/singleton.h>
#include <utils/factory.h>
#include <utils/maths.h>
#include <utils/asserter.h>

#include "advisors.h"
#include "history.h"
#include "propagator.h"

/**
 * @brief Decision Class
 * Holds a variable assignment
 */

class Decision
{
public:
	Decision(int vid = -1, double val = 0.0) : var(vid), value(val) {}
	int var;
	double value;
};

/**
 * Propagation Engine
 * Coordinates propagators actions, advisors and store the var domains
 */

class PropagationEngine
{
public:
	PropagationEngine() : stopPropagationIfFailed(false), domain(0), hasFailed(false) {}
	virtual ~PropagationEngine();
	inline DomainPtr getDomain() { return domain; }
	void setDomain(DomainPtr d);
	virtual void pushPropagator(PropagatorPtr prop);
	virtual bool propagate();
	virtual bool propagate(int var, double value);
	virtual bool propagate(const std::vector<int>& vars, const std::vector<double>& values);
	const std::vector<int>& getLastFixed() const { return lastFixed; }
	bool failed() const { return hasFailed; }
	// state handler
	StatePtr getStateMgr();
	// remove everything (advisors, propagators...)
	virtual void clear();
	// options
	bool stopPropagationIfFailed;
protected:
	// signal handlers
	void fixedBinUp(int j);
	void fixedBinDown(int j);
	void tightenedLb(int j, double newValue, double oldValue);
	void tightenedUb(int j, double newValue, double oldValue);
protected:
	friend class PropagationEngineState;

	DomainPtr domain;
	std::vector<int> vPropLbCount;
	std::vector<int> vPropUbCount;
	std::vector< std::vector<AdvisorPtr> > advisors;

	std::vector<PropagatorPtr> propagators;
	typedef std::deque<int> Queue;
	Queue queue;

	std::vector<Decision> decisions;
	std::vector<int> lastFixed;
	bool hasFailed;
	// helper
	PropagatorPtr top();
   void loop();
};

#endif /* PROP_ENGINE_H */
