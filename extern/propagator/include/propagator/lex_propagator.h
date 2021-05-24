/**
 * @file lex_propagator.h
 * @brief Propagators for the lexicographic global constraint
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2010
 */

#ifndef LEX_PROPAGATOR_H
#define LEX_PROPAGATOR_H

#include "prop_interface.h"

/**
 * Propagator for x <=_lex y
 * where x and y are vector of the same length n
 */

class LexLessEqualProp : public Propagator
{
public:
	LexLessEqualProp(ActiveDomain* d, const std::string& _name, const std::vector<int>& _x, const std::vector<int>& _y);
	const char* getTypeName() const { return "Lex<="; }
	void propagate();
	StatePtr getStateMgr();
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	friend class LexLessEqualPropAdvisor;
	friend class LexLessEqualPropState;
	int alpha;
	int beta;
	int n;
	std::vector<int> x;
	std::vector<int> y;
	// helpers
	bool fixedEqual(int i) const;
	bool checkLex(int i) const;
	void updateState(int i);
	void increaseAlpha();
	void decreaseBeta();
};

class LexLessEqualPropAdvisor : public AdvisorI
{
public:
	LexLessEqualPropAdvisor(LexLessEqualProp* p, int i) : AdvisorI(p), idx(i) {}
	// events for binary variables
	void fixedUp();
	void fixedDown();
	// events for other (general integer and continuous) variables
	void tightenLb(double delta, bool decreaseInfCnt, bool propagate);
	void tightenUb(double delta, bool decreaseInfCnt, bool propagate);
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	int idx;
};

#endif /* LEX_PROPAGATOR_H */
