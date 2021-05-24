/**
 * @file linear_propagator.cpp
 * @brief Linear propagators
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008-2012
 */

#include <iostream>

#include <fmt/format.h>
#include <utils/floats.h>

#include "propagator/linear_propagator.h"
#include "propagator/advisors.h"

using namespace dominiqs;

static const int LINEAR_DEFAULT_PRIORITY = 10000;
static const int KNAPSACK_DEFAULT_PRIORITY = 2000;
static const int CARDINALITY_DEFAULT_PRIORITY = 1000;

/**
 * Linear Constraint Propagator
 */

/**
 * Advise for changes in a variable appearing in the constraint
 * with a positive coefficient
 */

class PositiveLinearAdvisor : public AdvisorI
{
public:
	PositiveLinearAdvisor(LinearProp& p, int j, double coef) : AdvisorI(p, j), a(coef) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		p.minAct += a;
		p.dirty |= (lessThan(p.rhs, INFBOUND) && (p.minActInfCnt <= 1));
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		p.maxAct -= a;
		p.dirty |= (greaterThan(p.lhs, -INFBOUND) && (p.maxActInfCnt <= 1));
	}
	// events for other variables
	void tightenLb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		if (p.minActInfIdx == var) p.minActInfIdx = -1;
		p.minAct += delta * a;
		p.minActInfCnt -= decreaseInfCnt;
		p.dirty |= (propagate && lessThan(p.rhs, INFBOUND) && (p.minActInfCnt <= 1));
	}
	void tightenUb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		if (p.maxActInfIdx == var) p.maxActInfIdx = -1;
		p.maxAct += delta * a;
		p.maxActInfCnt -= decreaseInfCnt;
		p.dirty |= (propagate && greaterThan(p.lhs, -INFBOUND) && (p.maxActInfCnt <= 1));
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({}, +, idx={} coef={})", prop.getName(), var, a);
	}
protected:
	double a;
};

/**
 * Advise for changes in a variable appearing in the constraint
 * with a negative coefficient
 */

class NegativeLinearAdvisor : public AdvisorI
{
public:
	NegativeLinearAdvisor(LinearProp& p, int j, double coef) : AdvisorI(p, j), a(coef) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		p.maxAct += a;
		p.dirty |= (greaterThan(p.lhs, -INFBOUND) && (p.maxActInfCnt <= 1));
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		p.minAct -= a;
		p.dirty |= (lessThan(p.rhs, INFBOUND) && (p.minActInfCnt <= 1));
	}
	// events for other variables
	void tightenLb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		if (p.maxActInfIdx == var) p.maxActInfIdx = -1;
		p.maxAct += delta * a;
		p.maxActInfCnt -= decreaseInfCnt;
		p.dirty |= (propagate && greaterThan(p.lhs, -INFBOUND) && (p.maxActInfCnt <= 1));
	}
	void tightenUb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		LinearProp& p = getMyProp<LinearProp>();
		if (p.minActInfIdx == var) p.minActInfIdx = -1;
		p.minAct += delta * a;
		p.minActInfCnt -= decreaseInfCnt;
		p.dirty |= (propagate && lessThan(p.rhs, INFBOUND) && (p.minActInfCnt <= 1));
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({}, -, idx={} coef={})", prop.getName(), var, a);
	}
protected:
	double a;
};

class LinearPropState : public State
{
public:
	LinearPropState(LinearProp& p) : prop(p) {}
	std::shared_ptr<LinearPropState> clone() const
	{
		return std::make_shared<LinearPropState>(*this);
	}
	void dump()
	{
		minAct = prop.minAct;
		maxAct = prop.maxAct;
		minActInfCnt = prop.minActInfCnt;
		maxActInfCnt = prop.maxActInfCnt;
		state = prop.state;
	}
	void restore()
	{
		prop.minAct = minAct;
		prop.maxAct = maxAct;
		prop.minActInfCnt = minActInfCnt;
		prop.maxActInfCnt = maxActInfCnt;
		prop.minActInfIdx = -1;
		prop.maxActInfIdx = -1;
		prop.maxActDelta = -1.0;
		prop.state = state;
		prop.dirty = false;
	}
protected:
	LinearProp& prop;
	double minAct;
	double maxAct;
	int minActInfCnt;
	int maxActInfCnt;
	PropagatorState state;
};

LinearProp::LinearProp(Domain& d, Constraint* c) : Propagator(d)
{
	DOMINIQS_ASSERT( c );
	name = c->name;
	switch(c->sense)
	{
		case 'L':
			lhs = -INFBOUND;
			rhs = c->rhs;
			break;
		case 'E':
			lhs = c->rhs;
			rhs = c->rhs;
			break;
		case 'G':
			lhs = c->rhs;
			rhs = INFBOUND;
			break;
		case 'R':
			lhs = c->rhs - c->range;
			rhs = c->rhs;
			break;
		default:
			throw std::runtime_error("Unknown constraint sense!");
	}
	minAct = 0.0;
	maxAct = 0.0;
	minActInfCnt = 0;
	maxActInfCnt = 0;
	unsigned int rowSize = c->row.size();
	const int* idx = c->row.idx();
	const double* coef = c->row.coef();
	for (unsigned int k = 0; k < rowSize; k++)
	{
		int j = idx[k];
		double a = coef[k];
		if (isNull(a)) continue;
		if (a > 0.0)
		{
			if (domain.varType(j) == 'B')
			{
				maxAct += a;
				posBinIdx.push_back(j);
				posBinCoef.push_back(a);
			}
			else
			{
				if (lessThan(domain.varUb(j), INFBOUND)) maxAct += (domain.varUb(j) * a);
				else maxActInfCnt++;
				if (greaterThan(domain.varLb(j), -INFBOUND)) minAct += (domain.varLb(j) * a);
				else minActInfCnt++;
				posIdx.push_back(j);
				posCoef.push_back(a);
			}
		}
		else // a < 0.0
		{
			if (domain.varType(j) == 'B')
			{
				minAct += a;
				negBinIdx.push_back(j);
				negBinCoef.push_back(a);
			}
			else
			{
				if (lessThan(domain.varUb(j), INFBOUND)) minAct += (domain.varUb(j) * a);
				else minActInfCnt++;
				if (greaterThan(domain.varLb(j), -INFBOUND)) maxAct += (domain.varLb(j) * a);
				else maxActInfCnt++;
				negIdx.push_back(j);
				negCoef.push_back(a);
			}
		}
	}
	lastPosBin = posBinIdx.size();
	lastNegBin = negBinIdx.size();
	lastPos = posIdx.size();
	lastNeg = negIdx.size();
	updateState();
	// std::cout << name << " " << lastPosBin << " " << lastNegBin << " / " << lastPos << " " << lastNeg << std::endl;
	setPriority(LINEAR_DEFAULT_PRIORITY);
	// optimizations
	minActInfIdx = -1;
	maxActInfIdx = -1;
	maxActDelta = -1.0;
}

void LinearProp::createAdvisors(std::vector<AdvisorPtr>& advisors)
{
	unsigned int k;
	int j;
	double a;
	for (k = 0; k < lastPosBin; k++)
	{
		j = posBinIdx[k];
		a = posBinCoef[k];
		advisors.push_back(std::make_shared<PositiveLinearAdvisor>(*this, j, a));
	}
	for (k = 0; k < lastNegBin; k++)
	{
		j = negBinIdx[k];
		a = negBinCoef[k];
		advisors.push_back(std::make_shared<NegativeLinearAdvisor>(*this, j, a));
	}
	for (k = 0; k < lastPos; k++)
	{
		j = posIdx[k];
		a = posCoef[k];
		advisors.push_back(std::make_shared<PositiveLinearAdvisor>(*this, j, a));
	}
	for (k = 0; k < lastNeg; k++)
	{
		j = negIdx[k];
		a = negCoef[k];
		advisors.push_back(std::make_shared<NegativeLinearAdvisor>(*this, j, a));
	}
}

void LinearProp::updateState()
{
	if (state != CSTATE_UNKNOWN) return;
	if (((minActInfCnt == 0) && greaterThan(minAct, rhs)) || ((maxActInfCnt == 0) && lessThan(maxAct, lhs)))
	{
		state = CSTATE_INFEAS;
		dirty = false;
		return;
	}
	if (((minActInfCnt == 0) && greaterEqualThan(minAct, lhs)) && ((maxActInfCnt == 0) && lessEqualThan(maxAct, rhs)))
	{
		state = CSTATE_ENTAILED;
		dirty = false;
		return;
	}
}

void LinearProp::propagate()
{
	updateState();
	if (dirty == false) return;
	unsigned int k;
	int j;
	double a;
	double newB;
	dirty = false; // here at beginning: this is not monotonic
	// update maxActDelta if necessary
	if (maxActDelta < 0.0)
	{
		for (k = 0; k < lastPosBin; k++)
		{
			if (domain.isVarFixed(posBinIdx[k])) continue;
			if (maxActDelta < posBinCoef[k]) { maxActDelta = posBinCoef[k]; }
		}
		for (k = 0; k < lastNegBin; k++)
		{
			if (domain.isVarFixed(negBinIdx[k])) continue;
			if (maxActDelta < -negBinCoef[k]) { maxActDelta = -negBinCoef[k]; }
		}
		for (k = 0; k < lastPos; k++)
		{
			if (domain.isVarFixed(posIdx[k])) continue;
			double alpha = posCoef[k] * (domain.varUb(posIdx[k]) - domain.varLb(posIdx[k]));
			if (maxActDelta < alpha) { maxActDelta = alpha; }
		}
		for (k = 0; k < lastNeg; k++)
		{
			if (domain.isVarFixed(negIdx[k])) continue;
			double alpha = negCoef[k] * (domain.varLb(negIdx[k]) - domain.varUb(negIdx[k]));
			if (maxActDelta < alpha) { maxActDelta = alpha; }
		}
	}
	if (maxActDelta < 0.0)
	{
		 // may happen if all variables are fixed
		updateState();
		return;
	}
	double slack = lessThan(rhs, INFBOUND) ? rhs - minAct : INFBOUND;
	double surplus = greaterThan(lhs, -INFBOUND) ? maxAct - lhs : INFBOUND;
	if (lessEqualThan(maxActDelta, std::min(slack, surplus))) return; //< cannot deduce anything!
	maxActDelta = -1.0; //< will be recalculated on the next propagation!
	// one pass
	if (lessThan(rhs, INFBOUND))
	{
		double beta = rhs - minAct;
		if (minActInfCnt == 0)
		{
			for (k = 0; (k < lastPosBin) && (state == CSTATE_UNKNOWN); k++)
			{
				j = posBinIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = posBinCoef[k];
				if (greaterThan(a, beta)) domain.fixBinDown(j);
			}
			for (k = 0; (k < lastNegBin) && (state == CSTATE_UNKNOWN); k++)
			{
				j = negBinIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = negBinCoef[k];
				if (greaterThan(-a, beta)) domain.fixBinUp(j);
			}
			for (k = 0; (k < lastPos) && (state == CSTATE_UNKNOWN); k++)
			{
				j = posIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = posCoef[k];
				double alpha = a * (domain.varUb(j) - domain.varLb(j));
				if (greaterThan(alpha, beta))
				{
					newB = domain.varLb(j) + beta / a;
					if (domain.varType(j) != 'C') newB = floorEps(newB);
					domain.tightenUb(j, newB);
				}
			}
			for (k = 0; (k < lastNeg) && (state == CSTATE_UNKNOWN); k++)
			{
				j = negIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = negCoef[k];
				double alpha = a * (domain.varLb(j) - domain.varUb(j));
				if (greaterThan(alpha, beta))
				{
					newB = domain.varUb(j) + beta / a;
					if (domain.varType(j) != 'C') newB = ceilEps(newB);
					domain.tightenLb(j, newB);
				}
			}
		}
		if (minActInfCnt == 1)
		{
			if (minActInfIdx == -1) {
				for (k = 0; k < lastPos; k++) {
					j = posIdx[k];
					if (domain.isVarFixed(j)) continue;
					if (lessEqualThan(domain.varLb(j), -INFBOUND)) { minActInfIdx = j; minActInfCoef = posCoef[k]; break; }
				}
			}
			if (minActInfIdx == -1) {
				for (k = 0; k < lastNeg; k++) {
					j = negIdx[k];
					if (domain.isVarFixed(j)) continue;
					if (greaterEqualThan(domain.varUb(j), INFBOUND)) { minActInfIdx = j; minActInfCoef = negCoef[k]; break; }
				}
			}
			if ( minActInfIdx == -1 ) throw std::runtime_error("minActInfIdx to -1");
			if (minActInfCoef > 0.0)
			{
				newB = beta / minActInfCoef;
				if (domain.varType(minActInfIdx) != 'C') newB = floorEps(newB);
				if (lessThan(newB, domain.varUb(minActInfIdx)))
				{
					domain.tightenUb(minActInfIdx, newB);
				}
			}
			else // maxActInfCoef < 0.0
			{
				newB = beta / minActInfCoef;
				if (domain.varType(minActInfIdx) != 'C') newB = ceilEps(newB);
				if (greaterThan(newB, domain.varLb(minActInfIdx)))
				{
					domain.tightenLb(minActInfIdx, newB);
				}
			}
		}
	}
	if (state != CSTATE_UNKNOWN) return;
	if (greaterThan(lhs, -INFBOUND))
	{
		double beta = maxAct - lhs;
		if (maxActInfCnt == 0)
		{
			for (k = 0; (k < lastPosBin) && (state == CSTATE_UNKNOWN); k++)
			{
				j = posBinIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = posBinCoef[k];
				if (greaterThan(a, beta)) domain.fixBinUp(j);
			}
			for (k = 0; (k < lastNegBin) && (state == CSTATE_UNKNOWN); k++)
			{
				j = negBinIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = negBinCoef[k];
				if (greaterThan(-a, beta)) domain.fixBinDown(j);
			}
			for (k = 0; (k < lastPos) && (state == CSTATE_UNKNOWN); k++)
			{
				j = posIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = posCoef[k];
				double alpha = a * (domain.varUb(j) - domain.varLb(j));
				if (greaterThan(alpha, beta))
				{
					newB = domain.varUb(j) - beta / a;
					if (domain.varType(j) != 'C') newB = ceilEps(newB);
					domain.tightenLb(j, newB);
				}
			}
			for (k = 0; (k < lastNeg) && (state == CSTATE_UNKNOWN); k++)
			{
				j = negIdx[k];
				if (domain.isVarFixed(j)) continue;
				a = negCoef[k];
				double alpha = a * (domain.varLb(j) - domain.varUb(j));
				if (greaterThan(alpha, beta))
				{
					newB = domain.varLb(j) - beta / a;
					if (domain.varType(j) != 'C') newB = floorEps(newB);
					domain.tightenUb(j, newB);
				}
			}
		}
		if (maxActInfCnt == 1)
		{
			if (maxActInfIdx == -1) {
				for (k = 0; k < lastPos; k++) {
					j = posIdx[k];
					if (domain.isVarFixed(j)) continue;
					if (greaterEqualThan(domain.varUb(j), INFBOUND)) { maxActInfIdx = j; maxActInfCoef = posCoef[k]; break; }
				}
			}
			if (maxActInfIdx == -1) {
				for (k = 0; k < lastNeg; k++) {
					j = negIdx[k];
					if (domain.isVarFixed(j)) continue;
					if (lessEqualThan(domain.varLb(j), -INFBOUND)) { maxActInfIdx = j; maxActInfCoef = negCoef[k]; break; }
				}
			}
			if ( maxActInfIdx == -1 ) throw std::runtime_error("maxActInfIdx to -1");
			if (maxActInfCoef > 0.0)
			{
				newB = -beta / maxActInfCoef;
				if (domain.varType(maxActInfIdx) != 'C') newB = ceilEps(newB);
				if (greaterThan(newB, domain.varLb(maxActInfIdx)))
				{
					domain.tightenLb(maxActInfIdx, newB);
				}
			}
			else // maxActInfCoef < 0.0
			{
				newB = -beta / maxActInfCoef;
				if (domain.varType(maxActInfIdx) != 'C') newB = floorEps(newB);
				if (lessThan(newB, domain.varUb(maxActInfIdx)))
				{
					domain.tightenUb(maxActInfIdx, newB);
				}
			}
		}
	}
}

StatePtr LinearProp::getStateMgr()
{
	return std::make_shared<LinearPropState>(*this);
}

PropagatorFactoryPtr LinearFactory::clone() const
{
	return std::make_shared<LinearFactory>();
}

int LinearFactory::getPriority() const
{
	return LINEAR_DEFAULT_PRIORITY;
}

const char* LinearFactory::getName() const
{
	return "linear";
}

PropagatorPtr LinearFactory::analyze(Domain& d, Constraint* c)
{
	numCreated++;
	return std::make_shared<LinearProp>(d, c);
}

/**
 * Cardinality propagator
 */

class CardinalityAdvisor : public AdvisorI
{
public:
	CardinalityAdvisor(CardinalityProp& p, int j) : AdvisorI(p, j) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		CardinalityProp& p = getMyProp<CardinalityProp>();
		p.minAct++;
		if (p.minAct == p.rhs) p.dirty = true;
		else if (p.minAct > p.rhs)
		{
			p.state = CSTATE_INFEAS;
			p.dirty = false;
		}
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		CardinalityProp& p = getMyProp<CardinalityProp>();
		p.maxAct--;
		if (p.maxAct == p.lhs) p.dirty = true;
		else if (p.maxAct < p.lhs)
		{
			p.state = CSTATE_INFEAS;
			p.dirty = false;
		}
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({})", prop.getName());
	}
};

class CardinalityPropState : public State
{
public:
	CardinalityPropState(CardinalityProp& p) : prop(p) {}
	std::shared_ptr<CardinalityPropState> clone() const
	{
		return std::make_shared<CardinalityPropState>(*this);
	}

	void dump()
	{
		minAct = prop.minAct;
		maxAct = prop.maxAct;
		state = prop.state;
	}
	void restore()
	{
		prop.minAct = minAct;
		prop.maxAct = maxAct;
		prop.state = state;
		prop.dirty = false;
	}
protected:
	CardinalityProp& prop;
	int minAct;
	int maxAct;
	PropagatorState state;
};

CardinalityProp::CardinalityProp(Domain& d, Constraint* c) : Propagator(d)
{
	DOMINIQS_ASSERT( c );
	name = c->name;
	// assuming data is correct
	idx.resize(c->row.size());
	std::copy(c->row.idx(), c->row.idx() + c->row.size(), idx.begin());
	switch(c->sense)
	{
		case 'L':
			lhs = 0;
			rhs = double2int<int>(floorEps(c->rhs));
			break;
		case 'E':
			lhs = double2int<int>(c->rhs);
			rhs = lhs;
			break;
		case 'G':
			lhs = double2int<int>(ceilEps(c->rhs));
			rhs = idx.size();
			break;
		case 'R':
			lhs = double2int<int>(ceilEps(c->rhs - c->range));
			rhs = double2int<int>(floorEps(c->rhs));
			break;
		default:
			throw std::runtime_error("Unknown constraint sense!");
	}
	minAct = 0;
	int k = (int)c->row.size();
	maxAct = k;
	for (int i = 0; i < k; i++)
	{
		int j = idx[i];
		if (equal(domain.varLb(j), domain.varUb(j)))
		{
			if (isNull(domain.varLb(j))) maxAct--;
			else minAct++;
		}
	}
	updateState();
	setPriority(CARDINALITY_DEFAULT_PRIORITY);
}

void CardinalityProp::createAdvisors(std::vector<AdvisorPtr>& advisors)
{
	for (int i: idx) advisors.push_back(std::make_shared<CardinalityAdvisor>(*this, i));
}

void CardinalityProp::updateState()
{
	if (state != CSTATE_UNKNOWN) return;
	if ((minAct > rhs) || (maxAct < lhs))
	{
		dirty = false;
		state = CSTATE_INFEAS;
		return;
	}
	if ((minAct >= lhs) && (maxAct <= rhs))
	{
		dirty = false;
		state = CSTATE_ENTAILED;
		return;
	}
}

void CardinalityProp::propagate()
{
	if (dirty == false) return;
	if (minAct == rhs)
	{
		for (int j: idx)
		{
			if (!domain.isVarFixed(j)) domain.fixBinDown(j);
		}
	}
	if (maxAct == lhs)
	{
		for (int j: idx)
		{
			if (!domain.isVarFixed(j)) domain.fixBinUp(j);
		}
	}
	dirty = false; // at the end: this is monotonic
}

StatePtr CardinalityProp::getStateMgr()
{
	return std::make_shared<CardinalityPropState>(*this);
}

std::ostream& CardinalityProp::print(std::ostream& out) const
{
	return out << fmt::format("CardinalityProp({}, {}, {} <= sum <= {}, minAct={}, maxAct={})",
								name, PropagatorStateName[state], lhs, rhs, minAct, maxAct);
}

PropagatorFactoryPtr CardinalityFactory::clone() const
{
	return std::make_shared<CardinalityFactory>();
}

int CardinalityFactory::getPriority() const
{
	return CARDINALITY_DEFAULT_PRIORITY;
}

const char* CardinalityFactory::getName() const
{
	return "cardinality";
}

PropagatorPtr CardinalityFactory::analyze(Domain& d, Constraint* c)
{
	bool isCard = true;
	const int* idx = c->row.idx();
	const double* coef = c->row.coef();
	unsigned int size = c->row.size();
	for (unsigned int k = 0; (k < size) && isCard; k++)
	{
		if (d.varType(idx[k]) != 'B') isCard = false;
		if (different(coef[k], 1.0)) isCard = false;
	}
	if (isCard)
	{
		numCreated++;
		return std::make_shared<CardinalityProp>(d, c);
	}
	return nullptr;
}

/**
 * Knapsack  propagator
 */

class KnapsackPropState : public State
{
public:
	KnapsackPropState(KnapsackProp& p) : prop(p) {}
	KnapsackPropState* clone() const
	{
		return new KnapsackPropState(*this);
	}
	void dump()
	{
		minAct = prop.minAct;
		maxAct = prop.maxAct;
		state = prop.state;
	}
	void restore()
	{
		prop.minAct = minAct;
		prop.maxAct = maxAct;
		prop.maxActDelta = -1.0;
		prop.state = state;
		prop.dirty = false;
	}
protected:
	KnapsackProp& prop;
	double minAct;
	double maxAct;
	int minActInfCnt;
	int maxActInfCnt;
	PropagatorState state;
};

class KnapsackAdvisor : public AdvisorI
{
public:
	KnapsackAdvisor(KnapsackProp& p, int j, double coef) : AdvisorI(p, j), a(coef) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		KnapsackProp& p = getMyProp<KnapsackProp>();
		p.minAct += a;
		p.dirty |= lessThan(p.rhs, INFBOUND);
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		KnapsackProp& p = getMyProp<KnapsackProp>();
		p.maxAct -= a;
		p.dirty |= greaterThan(p.lhs, -INFBOUND);
	}
	// events for other variables
	void tightenLb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		KnapsackProp& p = getMyProp<KnapsackProp>();
		p.minAct += delta * a;
		p.dirty |= (propagate && lessThan(p.rhs, INFBOUND));
	}
	void tightenUb(double delta, bool decreaseInfCnt, bool propagate)
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		KnapsackProp& p = getMyProp<KnapsackProp>();
		p.maxAct += delta * a;
		p.dirty |= (propagate && greaterThan(p.lhs, -INFBOUND));
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({}, idx={} coef={})", prop.getName(), var, a);
	}
protected:
	double a;
};

KnapsackProp::KnapsackProp(Domain& d, Constraint* c) : Propagator(d)
{
	DOMINIQS_ASSERT( c );
	name = c->name;
	// assuming data is correct
	switch(c->sense)
	{
		case 'L':
			lhs = -INFBOUND;
			rhs = c->rhs;
			break;
		case 'E':
			lhs = c->rhs;
			rhs = c->rhs;
			break;
		case 'G':
			lhs = c->rhs;
			rhs = INFBOUND;
			break;
		case 'R':
			lhs = c->rhs - c->range;
			rhs = c->rhs;
			break;
		default:
			throw std::runtime_error("Unknown constraint sense!");
	}
	minAct = 0.0;
	maxAct = 0.0;
	unsigned int rowSize = c->row.size();
	const int* idx = c->row.idx();
	const double* coef = c->row.coef();
	for (unsigned int k = 0; k < rowSize; k++)
	{
		int j = idx[k];
		double a = coef[k];
		if (isNull(a)) continue;
		DOMINIQS_ASSERT( a > 0.0 );
		if (domain.varType(j) == 'B')
		{
			maxAct += a;
			posBinIdx.push_back(j);
			posBinCoef.push_back(a);
		}
		else
		{
			minAct += a * domain.varLb(j);
			maxAct += a * domain.varUb(j);
			posIdx.push_back(j);
			posCoef.push_back(a);
		}
	}
	lastPosBin = posBinIdx.size();
	lastPos = posIdx.size();
	maxActDelta = -1.0;
	updateState();
	setPriority(KNAPSACK_DEFAULT_PRIORITY);
}

void KnapsackProp::createAdvisors(std::vector<AdvisorPtr>& advisors)
{
	unsigned int k;
	int j;
	double a;
	for (k = 0; k < lastPosBin; k++)
	{
		j = posBinIdx[k];
		a = posBinCoef[k];
		advisors.push_back(std::make_shared<KnapsackAdvisor>(*this, j, a));
	}
	for (k = 0; k < lastPos; k++)
	{
		j = posIdx[k];
		a = posCoef[k];
		advisors.push_back(std::make_shared<KnapsackAdvisor>(*this, j, a));
	}
}

void KnapsackProp::updateState()
{
	if (state != CSTATE_UNKNOWN) return;
	if (greaterThan(minAct, rhs) || lessThan(maxAct, lhs))
	{
		state = CSTATE_INFEAS;
		dirty = false;
		return;
	}
	if (lessEqualThan(maxAct, rhs) && greaterEqualThan(minAct, lhs))
	{
		state = CSTATE_ENTAILED;
		dirty = false;
		return;
	}
}

void KnapsackProp::propagate()
{
	updateState();
	if (dirty == false) return;
	unsigned int k;
	int j;
	double a;
	double newB;
	dirty = false; // here at beginning because this is not monotonic
	// update maxActDelta if necessary
	if (maxActDelta < 0.0)
	{
		for (k = 0; k < lastPosBin; k++)
		{
			if (domain.isVarFixed(posBinIdx[k])) continue;
			if (maxActDelta < posBinCoef[k]) { maxActDelta = posBinCoef[k]; }
		}
		for (k = 0; k < lastPos; k++)
		{
			if (domain.isVarFixed(posIdx[k])) continue;
			double alpha = posCoef[k] * (domain.varUb(posIdx[k]) - domain.varLb(posIdx[k]));
			if (maxActDelta < alpha) { maxActDelta = alpha; }
		}
	}
	if (maxActDelta < 0.0)
	{
		 // may happen if all variables are fixed
		updateState();
		return;
	}
	double slack = lessThan(rhs, INFBOUND) ? (rhs - minAct) : INFBOUND;
	double surplus = greaterThan(lhs, -INFBOUND) ? (maxAct - lhs) : INFBOUND;
	if (lessEqualThan(maxActDelta, std::min(slack, surplus))) return; //< cannot deduce anything!
	maxActDelta = -1.0; //< will be recalculated on the next propagation!
	// one pass
	if (lessThan(rhs, INFBOUND))
	{
		double beta = rhs - minAct;
		for (k = 0; (k < lastPosBin) && (state == CSTATE_UNKNOWN); k++)
		{
			j = posBinIdx[k];
			if (domain.isVarFixed(j)) continue;
			a = posBinCoef[k];
			if (greaterThan(a, beta)) domain.fixBinDown(j);
		}
		for (k = 0; (k < lastPos) && (state == CSTATE_UNKNOWN); k++)
		{
			j = posIdx[k];
			if (domain.isVarFixed(j)) continue;
			a = posCoef[k];
			double alpha = a * (domain.varUb(j) - domain.varLb(j));
			if (greaterThan(alpha, beta))
			{
				newB = domain.varLb(j) + beta / a;
				if (domain.varType(j) != 'C') newB = floorEps(newB);
				domain.tightenUb(j, newB);
			}
		}
	}
	if (greaterThan(lhs, -INFBOUND))
	{
		double beta = maxAct - lhs;
		for (k = 0; (k < lastPosBin) && (state == CSTATE_UNKNOWN); k++)
		{
			j = posBinIdx[k];
			if (domain.isVarFixed(j)) continue;
			a = posBinCoef[k];
			if (greaterThan(a, beta)) domain.fixBinUp(j);
		}
		for (k = 0; (k < lastPos) && (state == CSTATE_UNKNOWN); k++)
		{
			j = posIdx[k];
			if (domain.isVarFixed(j)) continue;
			a = posCoef[k];
			double alpha = a * (domain.varUb(j) - domain.varLb(j));
			if (greaterThan(alpha, beta))
			{
				newB = domain.varUb(j) - beta / a;
				if (domain.varType(j) != 'C') newB = ceilEps(newB);
				domain.tightenLb(j, newB);
			}
		}
	}
}

StatePtr KnapsackProp::getStateMgr()
{
	return std::make_shared<KnapsackPropState>(*this);
}

std::ostream& KnapsackProp::print(std::ostream& out) const
{
	return out << fmt::format("KnapsackProp({}, {}, {} <= a^T x <= {}, minAct={}, maxAct={})",
								name, PropagatorStateName[state],
								(greaterThan(lhs, -INFBOUND) ? std::to_string(lhs) : "-inf"),
								(lessThan(rhs, INFBOUND) ? std::to_string(rhs) : "inf"),
								minAct, maxAct);
}

PropagatorFactoryPtr KnapsackFactory::clone() const
{
	return std::make_shared<KnapsackFactory>();
}

int KnapsackFactory::getPriority() const
{
	return KNAPSACK_DEFAULT_PRIORITY;
}

const char* KnapsackFactory::getName() const
{
	return "knapsack";
}

PropagatorPtr KnapsackFactory::analyze(Domain& d, Constraint* c)
{
	bool isKp = true;
	const int* idx = c->row.idx();
	const double* coef = c->row.coef();
	unsigned int size = c->row.size();
	for (unsigned int k = 0; (k < size) && isKp; k++)
	{
		if (isNegative(coef[k])) isKp = false;
		if (isNegative(d.varLb(idx[k]))) isKp = false;
		if (!lessThan(d.varUb(idx[k]), INFBOUND)) isKp = false;
	}
	if (isKp)
	{
		numCreated++;
		return std::make_shared<KnapsackProp>(d, c);
	}
	return nullptr;
}

// auto registration

class LINEAR_FACTORY_RECORDER
{
public:
	LINEAR_FACTORY_RECORDER()
	{
		std::cout << "Registering LinearPropagatorFactories...";
		PropagatorFactories::getInstance().registerClass<LinearFactory>("linear");
		PropagatorFactories::getInstance().registerClass<CardinalityFactory>("cardinality");
		PropagatorFactories::getInstance().registerClass<KnapsackFactory>("knapsack");
		std::cout << "done" << std::endl;
	}
};

LINEAR_FACTORY_RECORDER my_linear_factory_recorder;
