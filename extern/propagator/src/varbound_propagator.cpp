/**
 * @file varbound_propagator.cpp
 * @brief Propagators for variable bound constraints
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2012
 */

#include <iostream>

#include <utils/floats.h>
#include <fmt/format.h>

#include "propagator/varbound_propagator.h"
#include "propagator/advisors.h"

using namespace dominiqs;

class VarLowerBoundPropAdvisor : public AdvisorI
{
public:
	VarLowerBoundPropAdvisor(VarLowerBoundProp& p, int j) : AdvisorI(p, j) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarLowerBoundProp& p = getMyProp<VarLowerBoundProp>();
		if (isPositive(p.yCoef)) return;
		p.pendingLb = std::max(p.pendingLb, p.lb - p.yCoef);
		p.dirty = true;
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarLowerBoundProp& p = getMyProp<VarLowerBoundProp>();
		if (isNegative(p.yCoef)) return;
		p.pendingLb = std::max(p.pendingLb, p.lb);
		p.dirty = true;
	}
	// events for other variables
	void tightenLb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarLowerBoundProp& p = getMyProp<VarLowerBoundProp>();
		if (isPositive(p.yCoef)) return;
		p.pendingLb = std::max(p.pendingLb, p.lb - p.yCoef * p.domain.varLb(p.yIdx));
		p.dirty = propagate;
	}
	void tightenUb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarLowerBoundProp& p = getMyProp<VarLowerBoundProp>();
		if (isNegative(p.yCoef)) return;
		p.pendingLb = std::max(p.pendingLb, p.lb - p.yCoef * p.domain.varUb(p.yIdx));
		p.dirty = propagate;
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({}, vlb)", prop.getName());
	}
};

class VarLowerBoundPropState : public State
{
public:
	VarLowerBoundPropState(VarLowerBoundProp& p) : prop(p) {}
	std::shared_ptr<VarLowerBoundPropState> clone() const
	{
		return std::make_shared<VarLowerBoundPropState>(*this);
	}
	void dump()
	{
		state = prop.state;
		pendingLb = prop.pendingLb;
	}
	void restore()
	{
		prop.dirty = false;
		prop.state = state;
		prop.pendingLb = pendingLb;
	}
protected:
	VarLowerBoundProp& prop;
	PropagatorState state;
	double pendingLb;
};

VarLowerBoundProp::VarLowerBoundProp(Domain& d, const std::string& _name, int _xIdx, int _yIdx, double _yCoef, double _lb)
	: Propagator(d), xIdx(_xIdx), yIdx(_yIdx), yCoef(_yCoef), lb(_lb)
{
	name = _name;
	pendingLb = domain.varLb(xIdx);
	if (isNull(yCoef))
	{
		// numerically unstable case: mark it as entailed
		state = CSTATE_ENTAILED;
	}
}

void VarLowerBoundProp::createAdvisors(std::vector<AdvisorPtr>& advisors)
{
	advisors.push_back(std::make_shared<VarLowerBoundPropAdvisor>(*this, yIdx));
}

void VarLowerBoundProp::propagate()
{
	dirty = false;
	if (state == CSTATE_UNKNOWN)
	{
		if (domain.varType(xIdx) == 'I') pendingLb = ceilEps(pendingLb);
		if (greaterThan(pendingLb, domain.varLb(xIdx)))
		{
			domain.tightenLb(xIdx, pendingLb);
		}
	}
}

StatePtr VarLowerBoundProp::getStateMgr()
{
	return std::make_shared<VarLowerBoundPropState>(*this);
}

std::ostream& VarLowerBoundProp::print(std::ostream& out) const
{
	return out << fmt::format("VarLowerBoundProp({}, {}, {} + {} {} >= {})",
								name, PropagatorStateName[state], domain.varName(xIdx), yCoef, domain.varName(yIdx), lb);
}

class VarUpperBoundPropAdvisor : public AdvisorI
{
public:
	VarUpperBoundPropAdvisor(VarUpperBoundProp& p, int j) : AdvisorI(p, j) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarUpperBoundProp& p = getMyProp<VarUpperBoundProp>();
		if (isNegative(p.yCoef)) return;
		p.pendingUb = std::min(p.pendingUb, p.ub - p.yCoef);
		p.dirty = true;
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarUpperBoundProp& p = getMyProp<VarUpperBoundProp>();
		if (isPositive(p.yCoef)) return;
		p.pendingUb = std::min(p.pendingUb, p.ub);
		p.dirty = true;
	}
	// events for other variables
	void tightenLb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarUpperBoundProp& p = getMyProp<VarUpperBoundProp>();
		if (isNegative(p.yCoef)) return;
		p.pendingUb = std::min(p.pendingUb, p.ub - p.yCoef * p.domain.varLb(p.yIdx));
		p.dirty = propagate;
	}
	void tightenUb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		VarUpperBoundProp& p = getMyProp<VarUpperBoundProp>();
		if (isPositive(p.yCoef)) return;
		p.pendingUb = std::min(p.pendingUb, p.ub - p.yCoef * p.domain.varUb(p.yIdx));
		p.dirty = propagate;
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({}, vub)", prop.getName());
	}
};

class VarUpperBoundPropState : public State
{
public:
	VarUpperBoundPropState(VarUpperBoundProp& p) : prop(p) {}
	std::shared_ptr<VarUpperBoundPropState> clone() const
	{
		return std::make_shared<VarUpperBoundPropState>(*this);
	}
	void dump()
	{
		state = prop.state;
		pendingUb = prop.pendingUb;
	}
	void restore()
	{
		prop.dirty = false;
		prop.state = state;
		prop.pendingUb = pendingUb;
	}
protected:
	VarUpperBoundProp& prop;
	PropagatorState state;
	double pendingUb;
};

VarUpperBoundProp::VarUpperBoundProp(Domain& d, const std::string& _name, int _xIdx, int _yIdx, double _yCoef, double _ub)
	: Propagator(d), xIdx(_xIdx), yIdx(_yIdx), yCoef(_yCoef), ub(_ub)
{
	name = _name;
	pendingUb = domain.varUb(xIdx);
	if (isNull(yCoef))
	{
		// numerically unstable case: mark it as entailed
		state = CSTATE_ENTAILED;
	}
}

void VarUpperBoundProp::createAdvisors(std::vector<AdvisorPtr>& advisors)
{
	advisors.push_back(std::make_shared<VarUpperBoundPropAdvisor>(*this, yIdx));
}

void VarUpperBoundProp::propagate()
{
	dirty = false;
	if (state == CSTATE_UNKNOWN)
	{
		if (domain.varType(xIdx) == 'I') pendingUb = floorEps(pendingUb);
		if (lessThan(pendingUb, domain.varUb(xIdx)))
		{
			if (lessThan(pendingUb, domain.varLb(xIdx)))
			{
				state = CSTATE_INFEAS;
				return;
			}
			domain.tightenUb(xIdx, pendingUb);
		}
	}
}

StatePtr VarUpperBoundProp::getStateMgr()
{
	return std::make_shared<VarUpperBoundPropState>(*this);
}

std::ostream& VarUpperBoundProp::print(std::ostream& out) const
{
	return out << fmt::format("VarUpperBoundProp({}, {}, {} {} {} <= {})",
								name, PropagatorStateName[state], domain.varName(xIdx), yCoef, domain.varName(yIdx), ub);
}

PropagatorFactoryPtr VarBoundFactory::clone() const
{
	return std::make_shared<VarBoundFactory>();
}

int VarBoundFactory::getPriority() const
{
	return 20;
}

const char* VarBoundFactory::getName() const
{
	return "varbound";
}

PropagatorPtr VarBoundFactory::analyze(Domain& d, Constraint* c)
{
	if ((c->row.size() != 2) || (c->sense == 'E') || (c->sense == 'R')) return nullptr;
	DOMINIQS_ASSERT( (c->sense == 'L') || (c->sense == 'G') );
	const int* idx = c->row.idx();
	const double* coef = c->row.coef();
	// both binaries or both continuous are not allowed
	int numBin = (d.varType(idx[0]) == 'B') + (d.varType(idx[1]) == 'B');
	int numCont = (d.varType(idx[0]) == 'C') + (d.varType(idx[1]) == 'C');
	if ((numBin == 2) || (numCont == 2)) return nullptr;
	// map variables to x and y
	int xIdx = -1;
	double xCoef = 0.0;
	int yIdx = -1;
	double yCoef = 0.0;
	if (numBin)
	{
		if (d.varType(idx[0]) == 'B')
		{
			xIdx = idx[1];
			xCoef = coef[1];
			yIdx = idx[0];
			yCoef = coef[0];
		}
		else
		{
			xIdx = idx[0];
			xCoef = coef[0];
			yIdx = idx[1];
			yCoef = coef[1];
		}
	}
	else
	{
		if (d.varType(idx[0]) == 'I')
		{
			xIdx = idx[1];
			xCoef = coef[1];
			yIdx = idx[0];
			yCoef = coef[0];
		}
		else
		{
			xIdx = idx[0];
			xCoef = coef[0];
			yIdx = idx[1];
			yCoef = coef[1];
		}
	}
	// normalize coefficients
	char sense = c->sense;
	if (isNegative(xCoef))
	{
		// flip sense
		if (sense == 'L') sense = 'G';
		else sense = 'L';
	}
	xCoef = 1.0;
	yCoef /= xCoef;
	double rhs = c->rhs / xCoef;
   numCreated++;
	if (sense == 'L') return std::make_shared<VarUpperBoundProp>(d, c->name, xIdx, yIdx, yCoef, rhs);
	return std::make_shared<VarLowerBoundProp>(d, c->name, xIdx, yIdx, yCoef, rhs);
}

// auto registration

class VARBOUND_FACTORY_RECORDER
{
public:
	VARBOUND_FACTORY_RECORDER()
	{
		std::cout << "Registering VarBoundFactory...";
		PropagatorFactories::getInstance().registerClass<VarBoundFactory>("varbound");
		std::cout << "done" << std::endl;
	}
};

VARBOUND_FACTORY_RECORDER my_varbound_factory_recorder;
