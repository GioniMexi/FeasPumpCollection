/**
 * @file advisors.h
 * @brief Advisors Interface
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008-2012
 */

#ifndef ADVISORS_H
#define ADVISORS_H

#include <memory>
#include <iosfwd>

// forward declarations
class Propagator;
class PropagationEngine;

/**
 * @brief Propagator Advisor interface
 *
 * An advisor incrementally updates a propagator internal state
 * and eventually trigger its propagation in response to a variable domain change
 * Each advisor is linked to a single variable AND to a single propagator
 * Note that advisors are NOT allowed to change variable domains!
 */

class AdvisorI
{
public:
	inline AdvisorI(Propagator& p, int j) : prop(p), var(j) {}
	virtual ~AdvisorI() {}
	//@{
	inline Propagator& getPropagator() const { return prop; }
	inline int getVar() const { return var; }
	//@}

	//@{
	/**
	 * events for binary variables
	 */
	virtual void fixedUp() {}
	virtual void fixedDown() {}
	//@}
	//@{
	/**
	 * events for other (general integer and continuous) variables
	 * @param delta: bound increase (if it was finite) or new value if tightening from infinite
	 * @param decreaseInfCnt (this flag is true if the previous value for the bound was infinite)
	 * @param propagate: mask to trigger or not constraint propagation
	 */
	virtual void tightenLb(double delta, bool decreaseInfCnt, bool propagate) {}
	virtual void tightenUb(double delta, bool decreaseInfCnt, bool propagate) {}
	//@}
	// output (virtual output idiom)
	virtual std::ostream& print(std::ostream& out) const { return out; }
protected:
	Propagator& prop;
	int var;
	template<class T> inline T& getMyProp() { return static_cast<T&>(prop); }
};

inline std::ostream& operator<<(std::ostream& out, const AdvisorI& adv) { return adv.print(out); }

typedef std::shared_ptr<AdvisorI> AdvisorPtr;

#endif /* ADVISORS_H */
