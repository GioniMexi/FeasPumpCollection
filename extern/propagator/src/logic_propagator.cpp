/**
 * @file logic_propagator.cpp
 * @brief Propagators for logic constraints
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008-2012
 */

#include <iostream>

#include <utils/floats.h>
#include <fmt/format.h>

#include "propagator/logic_propagator.h"

using namespace dominiqs;

static const int FREE_IDX = 0;
static const int FALSE_IDX = 1;
static const int TRUE_IDX = 2;
static const int IMPLIES_STATE_CNT = 3;

static PropagatorState ImpliesState[IMPLIES_STATE_CNT][IMPLIES_STATE_CNT] =
{
	{CSTATE_UNKNOWN, CSTATE_UNKNOWN, CSTATE_ENTAILED},
	{CSTATE_ENTAILED, CSTATE_ENTAILED, CSTATE_STRONG_ENTAILED},
	{CSTATE_UNKNOWN, CSTATE_INFEAS, CSTATE_ENTAILED}
};

class ImpliesAntecedentAdvisor : public AdvisorI
{
public:
	ImpliesAntecedentAdvisor(ImpliesProp& p, int j) : AdvisorI(p, j) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		ImpliesProp& p = getMyProp<ImpliesProp>();
		p.antecedent = TRUE_IDX;
		p.state = ImpliesState[p.antecedent][p.consequent];
		p.dirty = (p.state == CSTATE_UNKNOWN);
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		ImpliesProp& p = getMyProp<ImpliesProp>();
		p.antecedent = FALSE_IDX;
		p.state = ImpliesState[p.antecedent][p.consequent];
		p.dirty = (p.state == CSTATE_UNKNOWN);
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({}, antecedent)", prop.getName());
	}
};

class ImpliesConsequentAdvisor : public AdvisorI
{
public:
	ImpliesConsequentAdvisor(ImpliesProp& p, int j) : AdvisorI(p, j) {}
	Propagator* getPropagator() const;
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		ImpliesProp& p = getMyProp<ImpliesProp>();
		p.consequent = TRUE_IDX;
		p.state = ImpliesState[p.antecedent][p.consequent];
		p.dirty = (p.state == CSTATE_UNKNOWN);
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		ImpliesProp& p = getMyProp<ImpliesProp>();
		p.consequent = FALSE_IDX;
		p.state = ImpliesState[p.antecedent][p.consequent];
		p.dirty = (p.state == CSTATE_UNKNOWN);
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({}, consequent)", prop.getName());
	}
};

class ImpliesPropState : public State
{
public:
	ImpliesPropState(ImpliesProp& p) : prop(p) {}
	std::shared_ptr<ImpliesPropState> clone() const
	{
		return std::make_shared<ImpliesPropState>(*this);
	}
	void dump()
	{
		antecedent = prop.antecedent;
		consequent = prop.consequent;
	}
	void restore()
	{
		prop.dirty = false;
		prop.antecedent = antecedent;
		prop.consequent = consequent;
		prop.state = ImpliesState[antecedent][consequent];
	}
protected:
	ImpliesProp& prop;
	int antecedent;
	int consequent;
};

ImpliesProp::ImpliesProp(Domain& d, const std::string& _name, int xIdx, int yIdx)
	: Propagator(d), anteIdx(xIdx), consIdx(yIdx)
{
	name = _name;
	antecedent = FREE_IDX;
	consequent = FREE_IDX;
	if (domain.isVarFixed(anteIdx))
	{
		if (isNull(domain.varLb(anteIdx))) antecedent = FALSE_IDX;
		else antecedent = TRUE_IDX;
	}
	if (domain.isVarFixed(consIdx))
	{
		if (isNull(domain.varLb(consIdx))) consequent = FALSE_IDX;
		else consequent = TRUE_IDX;
	}
	state = ImpliesState[antecedent][consequent];
}

void ImpliesProp::createAdvisors(std::vector<AdvisorPtr>& advisors)
{
	advisors.push_back(std::make_shared<ImpliesAntecedentAdvisor>(*this, anteIdx));
	advisors.push_back(std::make_shared<ImpliesConsequentAdvisor>(*this, consIdx));
}

void ImpliesProp::propagate()
{
	dirty = false;
	if (state == CSTATE_UNKNOWN)
	{
		if (antecedent == TRUE_IDX)
		{
			domain.fixBinUp(consIdx);
		}
		if (consequent == FALSE_IDX)
		{
			domain.fixBinDown(anteIdx);
		}
	}
}

StatePtr ImpliesProp::getStateMgr()
{
	return std::make_shared<ImpliesPropState>(*this);
}

std::ostream& ImpliesProp::print(std::ostream& out) const
{
	return out << fmt::format("ImpliesProp({}, {}, {} -> {})", name, PropagatorStateName[state], domain.varName(anteIdx), domain.varName(consIdx));
}

class EquivAdvisor : public AdvisorI
{
public:
	EquivAdvisor(EquivProp& p, int j) : AdvisorI(p, j) {}
	// events for binary variables
	void fixedUp()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		EquivProp& p = getMyProp<EquivProp>();
		p.updateState();
		p.dirty = (p.state == CSTATE_UNKNOWN);
	}
	void fixedDown()
	{
		if (prop.getState() != CSTATE_UNKNOWN) return;
		EquivProp& p = getMyProp<EquivProp>();
		p.updateState();
		p.dirty = (p.state == CSTATE_UNKNOWN);
	}
	// output
	std::ostream& print(std::ostream& out) const
	{
		return out << fmt::format("adv({})", prop.getName());
	}
};

class EquivPropState : public State
{
public:
	EquivPropState(EquivProp& p) : prop(p) {}
	std::shared_ptr<EquivPropState> clone() const
	{
		return std::make_shared<EquivPropState>(*this);
	}
	void dump()
	{
		state = prop.state;
	}
	void restore()
	{
		prop.dirty = false;
		prop.state = state;
	}
protected:
	EquivProp& prop;
	PropagatorState state;
};

EquivProp::EquivProp(Domain& d, const std::string& _name, int xIdx, int yIdx) : Propagator(d), firstIdx(xIdx), secondIdx(yIdx)
{
	name = _name;
	updateState();
}

void EquivProp::createAdvisors(std::vector<AdvisorPtr>& advisors)
{
	advisors.push_back(std::make_shared<EquivAdvisor>(*this, firstIdx));
	advisors.push_back(std::make_shared<EquivAdvisor>(*this, secondIdx));
}

void EquivProp::propagate()
{
	dirty = false;
	if (state == CSTATE_UNKNOWN)
	{
		if (domain.isVarFixed(firstIdx) && !domain.isVarFixed(secondIdx))
		{
			if (isNull(domain.varLb(firstIdx))) domain.fixBinDown(secondIdx);
			else domain.fixBinUp(secondIdx);
		}
	}
	if (state == CSTATE_UNKNOWN)
	{
		if (domain.isVarFixed(secondIdx) && !domain.isVarFixed(firstIdx))
		{
			if (isNull(domain.varLb(secondIdx))) domain.fixBinDown(firstIdx);
			else domain.fixBinUp(firstIdx);
		}
	}
}

StatePtr EquivProp::getStateMgr()
{
	return std::make_shared<EquivPropState>(*this);
}

std::ostream& EquivProp::print(std::ostream& out) const
{
	return out << fmt::format("EquivProp({}, {}, {} <-> {})", name, PropagatorStateName[state], domain.varName(firstIdx), domain.varName(secondIdx));
}

void EquivProp::updateState()
{
	bool firstFixed = domain.isVarFixed(firstIdx);
	bool secondFixed = domain.isVarFixed(secondIdx);
	if (firstFixed && secondFixed)
	{
		if (equal(domain.varLb(firstIdx), domain.varLb(secondIdx))) state = CSTATE_ENTAILED;
		else state = CSTATE_INFEAS;
	}
	else state = CSTATE_UNKNOWN;
}

PropagatorFactoryPtr LogicFactory::clone() const
{
	return std::make_shared<LogicFactory>();
}

int LogicFactory::getPriority() const
{
	return 10;
}

const char* LogicFactory::getName() const
{
	return "logic";
}

PropagatorPtr LogicFactory::analyze(Domain& d, Constraint* c)
{
	if (c->sense == 'R') return nullptr;
	bool isImplies = (c->row.size() == 2) && (c->sense != 'G') && isNull(c->rhs);
	const int* idx = c->row.idx();
	const double* coef = c->row.coef();
	if (isImplies)
	{
		if (d.varType(idx[0]) != 'B') return nullptr;
		if (d.varType(idx[1]) != 'B') return nullptr;
		double minCoef = std::min(coef[0], coef[1]);
		double maxCoef = std::max(coef[0], coef[1]);
		if (equal(minCoef, -1.0) && equal(maxCoef, 1.0))
		{
			numCreated++;
			if (c->sense == 'L')
			{
				if (coef[0] > 0) return std::make_shared<ImpliesProp>(d, c->name, idx[0], idx[1]);
				else return std::make_shared<ImpliesProp>(d, c->name, idx[1], idx[0]);
			}
			else if (c->sense == 'E') return std::make_shared<EquivProp>(d, c->name, idx[0], idx[1]);
		}
	}
	return nullptr;
}

// auto registration

class LOGIC_FACTORY_RECORDER
{
public:
	LOGIC_FACTORY_RECORDER()
	{
		std::cout << "Registering LogicFactory...";
		PropagatorFactories::getInstance().registerClass<LogicFactory>("logic");
		std::cout << "done" << std::endl;
	}
};

LOGIC_FACTORY_RECORDER my_logic_factory_recorder;
