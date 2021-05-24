/**
 * @file domain.h
 * @brief Variables' domains Header
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef DOMAIN_H
#define DOMAIN_H

#include <vector>
#include <string>
#include <functional>

#include <utils/floats.h>
#include <utils/asserter.h>

#include "history.h"

#define INFBOUND 1e20

// forward declaration
class DomainState;

/**
 * Stores the domains of a set of variables and their info
 */

class Domain
{
public:
	~Domain() { clear(); }
	/**
	 * Add a variable to the domain
	 * @param name variable name
	 * @param t variable type
	 * @param l variable lower bound
	 * @param u variable upper bound
	 */
	virtual void pushVar(const std::string& name, char t, double l, double u);
	virtual void clear();
	//@{
	// getters
	inline unsigned int size() const { return names.size(); }
	inline std::string varName(int j) const { return names[j]; }
	inline double varLb(int j) const { return lb[j]; }
	inline double varUb(int j) const { return ub[j]; }
	inline bool isVarFixed(int j) const { return fixed[j]; }
	inline char varType(int j) const { return type[j]; }
	//@}
	//@{
	// setters
	inline void fixBinUp(int j)
	{
		DOMINIQS_ASSERT( dominiqs::equal(ub[j], 1.0) );
		DOMINIQS_ASSERT( type[j] == 'B' );
		lb[j] = 1.0;
		fixed[j] = true;
		if (emitFixedBinUp) emitFixedBinUp(j);
	}
	inline void fixBinDown(int j)
	{
		DOMINIQS_ASSERT( dominiqs::equal(lb[j], 0.0) );
		DOMINIQS_ASSERT( type[j] == 'B' );
		ub[j] = 0.0;
		fixed[j] = true;
		if (emitFixedBinDown) emitFixedBinDown(j);
	}
	inline void tightenLb(int j, double newValue)
	{
		DOMINIQS_ASSERT( type[j] != 'B' );
		double oldValue = lb[j];
		newValue = std::min(newValue, ub[j]);
		if (dominiqs::greaterThan(newValue, oldValue))
		{
			lb[j] = newValue;
			if (dominiqs::isNull(ub[j] - lb[j])) fixed[j] = true;
			if (emitTightenedLb) emitTightenedLb(j, newValue, oldValue);
		}
	}
	inline void tightenUb(int j, double newValue)
	{
		DOMINIQS_ASSERT( type[j] != 'B' );
		double oldValue = ub[j];
		newValue = std::max(newValue, lb[j]);
		if (dominiqs::lessThan(newValue, oldValue))
		{
			ub[j] = newValue;
			if (dominiqs::isNull(ub[j] - lb[j])) fixed[j] = true;
			if (emitTightenedUb) emitTightenedUb(j, newValue, oldValue);
		}
	}
	//@}
	//@{
	// callbacks
	std::function<void (int)> emitFixedBinUp;
	std::function<void (int)> emitFixedBinDown;
	std::function<void (int, double, double)> emitTightenedLb;
	std::function<void (int, double, double)> emitTightenedUb;
	//@}
	StatePtr getStateMgr();
protected:
	friend class DomainState;
	std::vector<std::string> names;
	std::vector<double> lb;
	std::vector<double> ub;
	std::vector<bool> fixed;
	std::vector<char> type;
};

typedef std::shared_ptr<Domain> DomainPtr;

/**
 * State Management
 */

class DomainState : public State
{
public:
	DomainState(Domain& d) : domain(d) {}
	std::shared_ptr<DomainState> clone() const;
	void dump();
	void restore();
private:
	Domain& domain;
	std::vector<double> lb;
	std::vector<double> ub;
	std::vector<bool> fixed;
};

typedef std::shared_ptr<DomainState> DomainStatePtr;

#endif /* DOMAIN_H */
