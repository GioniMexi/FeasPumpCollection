/**
 * @file logic_propagator.h
 * @brief Propagators for logic constraints
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008-2012
 */

#ifndef LOGIC_PROPAGATOR_H
#define LOGIC_PROPAGATOR_H

#include "propagator.h"
#include "advisors.h"


/**
 * Propagator for Implies constraint: x -> y
 */

class ImpliesProp : public Propagator
{
public:
	ImpliesProp(Domain& d, const std::string& _name, int xIdx, int yIdx);
	void createAdvisors(std::vector<AdvisorPtr>& advisors);
	void propagate();
	StatePtr getStateMgr();
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	friend class ImpliesAntecedentAdvisor;
	friend class ImpliesConsequentAdvisor;
	friend class ImpliesPropState;
	int antecedent;
	int consequent;
	int anteIdx;
	int consIdx;
};


/**
 * Propagator for Equivalence constraint: x <-> y
 */

class EquivProp : public Propagator
{
public:
	EquivProp(Domain& d, const std::string& _name, int xIdx, int yIdx);
	void createAdvisors(std::vector<AdvisorPtr>& advisors);
	void propagate();
	StatePtr getStateMgr();
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	void updateState();
	friend class EquivAdvisor;
	friend class EquivPropState;
	int firstIdx;
	int secondIdx;
};


/**
 * Factory for both
 */

class LogicFactory : public PropagatorFactory
{
public:
	std::shared_ptr<PropagatorFactory> clone() const;
	int getPriority() const;
	const char* getName() const;
	PropagatorPtr analyze(Domain& d, dominiqs::Constraint* c);
};


#endif /* LOGIC_PROPAGATOR_H */
