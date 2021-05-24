/**
 * @file linear_propagator.h
 * @brief Linear propagators
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008-2012
 */

#ifndef LINEAR_PROPAGATOR_H
#define LINEAR_PROPAGATOR_H

#include <vector>

#include "propagator.h"

/**
 * @brief Propagator for the linear (ranged) constraint lhs <= a^T x <= rhs
 *
 * lhs or rhs can be infinite or even equal (equality constraint)
 * It is used as fallback propagator if we cannot detect some specific
 * structure in the linear constraint
 */

class LinearProp : public Propagator
{
public:
	LinearProp(Domain& d, dominiqs::Constraint* c);
	void createAdvisors(std::vector<AdvisorPtr>& advisors);
	void propagate();
	StatePtr getStateMgr();
protected:
	friend class PositiveLinearAdvisor;
	friend class NegativeLinearAdvisor;
	friend class LinearPropState;
	double lhs;
	double rhs;
	double minAct;
	double maxAct;
	int minActInfCnt;
	int maxActInfCnt;
	std::vector<int> posBinIdx;
	std::vector<double> posBinCoef;
	unsigned int lastPosBin;
	std::vector<int> negBinIdx;
	std::vector<double> negBinCoef;
	unsigned int lastNegBin;
	std::vector<int> posIdx;
	std::vector<double> posCoef;
	unsigned int lastPos;
	std::vector<int> negIdx;
	std::vector<double> negCoef;
	unsigned int lastNeg;
	// optimizations
	int minActInfIdx;
	double minActInfCoef;
	int maxActInfIdx;
	double maxActInfCoef;
	double maxActDelta;
	// helpers
	void updateState();
};

class LinearFactory : public PropagatorFactory
{
public:
	PropagatorFactoryPtr clone() const;
	int getPriority() const;
	const char* getName() const;
	PropagatorPtr analyze(Domain& d, dominiqs::Constraint* c);
};

/**
 * @brief Propagator for the cardinality constraint
 * lhs <= e x <= rhs, all x_j binary
 *
 * Covers set covering/packing/partitioning as special cases
 */

class CardinalityProp : public Propagator
{
public:
	CardinalityProp(Domain& d, dominiqs::Constraint* c);
	void createAdvisors(std::vector<AdvisorPtr>& advisors);
	void propagate();
	StatePtr getStateMgr();
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	friend class CardinalityAdvisor;
	friend class CardinalityPropState;
	int lhs;
	int rhs;
	std::vector<int> idx;
	int minAct;
	int maxAct;
	// helpers
	void updateState();
};

class CardinalityFactory : public PropagatorFactory
{
public:
	PropagatorFactoryPtr clone() const;
	int getPriority() const;
	const char* getName() const;
	PropagatorPtr analyze(Domain& d, dominiqs::Constraint* c);
};

/**
 * @brief This is a propagator for the knapsack constraint:
 * lhs <= a^T x <= rhs, a,lhs,rhs > 0, 0 <= l <= x <= u
 */

class KnapsackProp : public Propagator
{
public:
	KnapsackProp(Domain& d, dominiqs::Constraint* c);
	void createAdvisors(std::vector<AdvisorPtr>& advisors);
	void propagate();
	StatePtr getStateMgr();
	// output
	std::ostream& print(std::ostream& out) const;
protected:
	friend class KnapsackAdvisor;
	friend class KnapsackPropState;
	double lhs;
	double rhs;
	double minAct;
	double maxAct;
	std::vector<int> posBinIdx;
	std::vector<double> posBinCoef;
	unsigned int lastPosBin;
	std::vector<int> posIdx;
	std::vector<double> posCoef;
	unsigned int lastPos;
	// optimizations
	double maxActDelta;
	// helpers
	void updateState();
};

class KnapsackFactory : public PropagatorFactory
{
public:
	PropagatorFactoryPtr clone() const;
	int getPriority() const;
	const char* getName() const;
	PropagatorPtr analyze(Domain& d, dominiqs::Constraint* c);
};

#endif /* LINEAR_PROPAGATOR_H */
