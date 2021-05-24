/**
 * @file varbound_propagator.h
 * @brief Propagators for variable bound constraints
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2012
 */

#ifndef VARBOUND_PROPAGATOR_H
#define VARBOUND_PROPAGATOR_H

#include "propagator.h"

/**
 * @brief Propagator for Variable Lower Bound constraint
 *
 * Propagates a constraint of the form x + c y >= lb, with
 * x non-binary and y non-continuous
 */

class VarLowerBoundProp : public Propagator
{
public:
	VarLowerBoundProp(Domain& d, const std::string& _name, int _xIdx, int _yIdx, double _yCoef, double _lb);
	void createAdvisors(std::vector<AdvisorPtr>& advisors);
	void propagate();
	StatePtr getStateMgr();
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	friend class VarLowerBoundPropAdvisor;
	friend class VarLowerBoundPropState;
	// data
	int xIdx;
	int yIdx;
	double yCoef;
	double lb;
	double pendingLb;
};

/**
 * @brief Propagator for Variable Upper Bound constraint
 *
 * Propagates a constraint of the form x + c y <= ub, with
 * x non-binary and y non-continuous
 */

class VarUpperBoundProp : public Propagator
{
public:
	VarUpperBoundProp(Domain& d, const std::string& _name, int _xIdx, int _yIdx, double _yCoef, double _ub);
	void createAdvisors(std::vector<AdvisorPtr>& advisors);
	void propagate();
	StatePtr getStateMgr();
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	friend class VarUpperBoundPropAdvisor;
	friend class VarUpperBoundPropState;
	// data
	int xIdx;
	int yIdx;
	double yCoef;
	double ub;
	double pendingUb;
};

/**
 * Factory for both
 */

class VarBoundFactory : public PropagatorFactory
{
public:
	PropagatorFactoryPtr clone() const;
	int getPriority() const;
	const char* getName() const;
	PropagatorPtr analyze(Domain& d, dominiqs::Constraint* c);
};

#endif /* VARBOUND_PROPAGATOR_H */
