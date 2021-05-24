/**
 * @file lex_propagator.cpp
 * @brief Propagators for the lexicographic global constraint
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2010
 */

#include <iostream>

#include <fmt/format.h>
#include <utils/floats.h>

#include "propagator/lex_propagator.h"

using namespace dominiqs;

void LexLessEqualPropAdvisor::fixedUp()
{
	if (prop->getState() != CSTATE_UNKNOWN) return;
	LexLessEqualProp* p = getMyProp<LexLessEqualProp>();
	p->updateState(idx);
	p->dirty = (p->state == CSTATE_UNKNOWN);
}

void LexLessEqualPropAdvisor::fixedDown()
{
	if (prop->getState() != CSTATE_UNKNOWN) return;
	LexLessEqualProp* p = getMyProp<LexLessEqualProp>();
	p->updateState(idx);
	p->dirty = (p->state == CSTATE_UNKNOWN);
}

void LexLessEqualPropAdvisor::tightenLb(double delta, bool decreaseInfCnt, bool propagate)
{
	if (prop->getState() != CSTATE_UNKNOWN) return;
	LexLessEqualProp* p = getMyProp<LexLessEqualProp>();
	p->updateState(idx);
	p->dirty = propagate && (p->state == CSTATE_UNKNOWN);
}

void LexLessEqualPropAdvisor::tightenUb(double delta, bool decreaseInfCnt, bool propagate)
{
	if (prop->getState() != CSTATE_UNKNOWN) return;
	LexLessEqualProp* p = getMyProp<LexLessEqualProp>();
	p->updateState(idx);
	p->dirty = propagate && (p->state == CSTATE_UNKNOWN);
}

std::ostream& LexLessEqualPropAdvisor::print(std::ostream& out) const
{
	return out << fmt::format("adv({}, {})", prop->getName(), idx);
}

class LexLessEqualPropState : public State
{
public:
	LexLessEqualPropState(LexLessEqualProp* p) : prop(p) {}
	LexLessEqualPropState* clone() const
	{
		return new LexLessEqualPropState(*this);
	}
	void dump()
	{
		state = prop->state;
		alpha = prop->alpha;
		beta = prop->beta;
	}
	void restore()
	{
		prop->dirty = false;
		prop->state = state;
		prop->alpha = alpha;
		prop->beta = beta;
	}
protected:
	LexLessEqualProp* prop;
	PropagatorState state;
	int alpha;
	int beta;
};

LexLessEqualProp::LexLessEqualProp(ActiveDomain* d, const std::string& _name, const std::vector<int>& _x, const std::vector<int>& _y)
	: Propagator(d), x(_x), y(_y)
{
	name = _name;
	DOMINIQS_ASSERT( domain );
	DOMINIQS_ASSERT( x.size() == y.size() );
	for (int j: x)
	{
		if (domain->varType(j) == 'C') throw std::runtime_error("Cannot have continuous variables in lex constraint!");
	}
	for (int j: y)
	{
		if (domain->varType(j) == 'C') throw std::runtime_error("Cannot have continuous variables in lex constraint!");
	}
	n = x.size();
	alpha = 0;
	beta = -1;
	// setup alpha
	while ((alpha < n) && fixedEqual(alpha)) ++alpha;
	if ((alpha >= n) || checkLex(alpha))
	{
		state = CSTATE_ENTAILED;
		return;
	}
	// setup beta
	int i = alpha;
	for (;i < n; i++)
	{
		if (greaterThan(domain->varLb(x[i]), domain->varUb(y[i]))) break;
		if (equal(domain->varLb(x[i]), domain->varUb(y[i])))
		{
			if (beta == -1) beta = i;
		}
		else beta = -1;
	}
	if (i == n) beta = n + 1;
	else if (beta == -1) beta = i;
	// std::cout << fmt::format("alpha = {} beta = {}", alpha, beta) << std::endl;
	if (alpha >= beta)
	{
		state = CSTATE_INFEAS;
		return;
	}
	// setup advisors
	for (unsigned int i = 0; i < x.size(); i++)
	{
		domain->pushAdvisor(x[i], new LexLessEqualPropAdvisor(this, i));
		domain->pushAdvisor(y[i], new LexLessEqualPropAdvisor(this, i));
	}
}

void LexLessEqualProp::propagate()
{
	dirty = false;
	if (state == CSTATE_UNKNOWN)
	{
		int xA = x[alpha];
		int yA = y [alpha];
		double xLB = domain->varLb(xA);
		double xUB = domain->varUb(xA);
		double yLB = domain->varLb(yA);
		double yUB = domain->varUb(yA);
		if (alpha == beta - 1)
		{
			// enforce x_alpha < y_alpha
			if (greaterEqualThan(xUB, yUB))
			{
				// impose x_alpha.ub < y_alpha.ub
				double newUB = yUB - 1.0;
				if (lessThan(newUB, xLB)) { state = CSTATE_INFEAS; return; }
				if (domain->varType(xA) == 'B') domain->fixBinDown(xA);
				else domain->tightenUb(xA, newUB);
			}
			if (lessEqualThan(yLB, xLB))
			{
				// impose y_alpha.lb < x_alpha.lb
				double newLB = xLB + 1.0;
				if (greaterThan(newLB, yUB)) { state = CSTATE_INFEAS; return; }
				if (domain->varType(yA) == 'B') domain->fixBinUp(yA);
				else domain->tightenLb(yA, newLB);
			}
		}
		else if (alpha < (beta - 1))
		{
			// enforce x_alpha <= y_alpha
			if (greaterThan(xUB, yUB))
			{
				// impose x_alpha.ub <= y_alpha.ub
				double newUB = yUB;
				if (lessThan(newUB, xLB)) { state = CSTATE_INFEAS; return; }
				if (domain->varType(xA) == 'B') domain->fixBinDown(xA);
				else domain->tightenUb(xA, newUB);
			}
			if (lessThan(yLB, xLB))
			{
				// impose y_alpha.lb <= x_alpha.lb
				double newLB = xLB;
				if (greaterThan(newLB, yUB)) { state = CSTATE_INFEAS; return; }
				if (domain->varType(yA) == 'B') domain->fixBinUp(yA);
				else domain->tightenLb(yA, newLB);
			}
		}
	}
}

void LexLessEqualProp::increaseAlpha()
{
	int stop = (beta < n) ? beta : n;
	while ((alpha < stop) && fixedEqual(alpha)) ++alpha;
}

void LexLessEqualProp::decreaseBeta()
{
	while ((beta > alpha) && equal(domain->varLb(x[beta-1]), domain->varUb(y[beta-1]))) --beta;
}

StatePtr LexLessEqualProp::getStateMgr()
{
	return new LexLessEqualPropState(this);
}

bool LexLessEqualProp::fixedEqual(int i) const
{
	return (domain->isVarFixed(x[i]) && domain->isVarFixed(y[i]) && (domain->varLb(x[i]) == domain->varLb(y[i])));
}

bool LexLessEqualProp::checkLex(int i) const
{
	return ( lessThan(domain->varUb(x[i]), domain->varLb(y[i]))
		|| ((i == (n - 1)) && lessEqualThan(domain->varUb(x[i]), domain->varLb(y[i]))) );
}

void LexLessEqualProp::updateState(int i)
{
	DOMINIQS_ASSERT( (i >= 0) && (i < n) );
	if ((i < alpha) || (i >= beta)) return;
	// update alpha or beta
	if (i == alpha) increaseAlpha();
	else if (i == beta - 1)
	{
		DOMINIQS_ASSERT( (beta > 0) && (beta - 1) < n );
		if (equal(domain->varLb(x[beta-1]), domain->varUb(y[beta-1])))
		{
			beta = i;
			decreaseBeta();
		}
	}
	else // alpha < i < beta - 1
	{
		if (greaterThan(domain->varLb(x[i]), domain->varUb(y[i])))
		{
			beta = i;
			decreaseBeta();
		}
	}
	// detect entailment or failure
	if (alpha == beta) state = CSTATE_INFEAS;
	if (alpha == n) state = CSTATE_ENTAILED;
	if ((alpha == beta - 1) && checkLex(alpha)) state = CSTATE_ENTAILED;
}

std::ostream& LexLessEqualProp::print(std::ostream& out) const
{
	return out << fmt::format("LexLessEqualProp({}, {}, alpha={}, beta={}, n={})",
							name, PropagatorStateName[state], alpha, beta, n);
}
