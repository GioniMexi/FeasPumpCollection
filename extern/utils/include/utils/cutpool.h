/**
 * @file cutpool.h
 * @brief Cut & CutPool
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef CUTPOOL_H
#define CUTPOOL_H

#include <assert.h>
#include <unordered_set>

#include <stdint.h>
#include <list>
#include <iostream>

#include "maths.h"
#include "fileconfig.h"

namespace dominiqs {

/**
 * This represent a row (cut) in the constraint matrix Ax ~ b
 */

class Cut : public Constraint
{
public:
	Cut() : slackType('U'), removable(false), inUse(false), age(0), separator(0),
			rank(0), efficacy(0.0), ortho(1.0), sig(0), norm(0.0) {}
	Cut(const Constraint& other, char sType) : Constraint(other),
			slackType(sType), removable(false), inUse(false), age(0),
			rank(0), efficacy(0.0), ortho(1.0), sig(0), norm(0.0) {}
	Cut* clone() const;
	char slackType;
	bool removable;
	bool inUse;
	uint8_t age;
	uint8_t separator;
	int rank; //< upper bound on cut rank
	double efficacy; //< violation divided by cut norm (=geometric distance)
	double ortho;
	std::size_t sig; //< signature (hash) to speed up comparisons
	double norm; //< euclidian norm of the cut
#ifdef DOMINIQS_TRACE_CUTS
	std::vector<double> mult; //< multipliers that generated the cut
	std::vector<int> complemented; //< indexes of complemented variables
#endif // DOMINIQS_TRACE_CUTS
	/**
	 * Get the score of the cut (function of efficacy and ortho)
	 */
	double score()
	{
		return efficacy + ortho;
	}
	/**
	 * Get the score of the cut for purging (function of efficacy and age)
	 * (the highest the score, the worst the cut!)
	 */
	double purgeScore()
	{
		return age - efficacy / (fabs(efficacy) + 1) - 10.0 * inUse;
	}
	/**
	 * Computes the dynamism of the cut
	 */
	double dynamism() const;
	/**
	 * Computes signatures and norm of the cut
	 */
	void digest();
	/**
	 * Compute memory occupation of the cut
	 */
	int memoryUsed() const;
	/**
	 * Save cut to stream (text form)
	 */
	void save(std::ostream& out) const;
	/**
	 * Read cut from stream (text form)
	 */
	void read(std::istream& in);
};

typedef std::shared_ptr<Cut> CutPtr;

struct cutptr_equal : std::binary_function<CutPtr, CutPtr, bool>
{
	bool operator()(const CutPtr& x, const CutPtr& y) const
	{
		if ((x->sig == y->sig) && equal(x->rhs, y->rhs) && (x->row.size() == y->row.size()))
		{
			// trivial checks passed: need a full blown comparison!
			if (x->row == y->row) return true;
		}
		return false;
	}
};

struct cutptr_hash : std::unary_function<CutPtr, std::size_t>
{
	std::size_t operator()(const CutPtr& x) const { return x->sig; }
};


/**
 * A list of cuts
 */

typedef std::vector<CutPtr> CutList;

/**
 * This class implements a cut pool, with a corresponding selection strategy
 */

class CutPool
{
public:
	CutPool();
	void readConfig(FileConfig& config, const std::string& root);
	bool push(CutPtr cut);
	void merge(CutPool& other);
	bool isViolated(const double* x) const;
	int numViolated(const double* x) const;
	bool select(const std::vector<double>& x, CutList& cuts);
	int purge(const std::vector<double>& x);
	void clear() { pool.clear(); nadded = 0; nduplicate = 0; mem = 0; }
	unsigned int size() const { return pool.size(); }
	unsigned int numAdded() const { return nadded; }
	unsigned int numDuplicate() const { return nduplicate; }
	double memoryUsed() const { return mem; }
	// text I/O
	void save(std::ostream& out) const;
	void read(std::istream& in);
	// iterators
	typedef std::unordered_set<CutPtr, cutptr_hash, cutptr_equal>::const_iterator const_iterator;
	const_iterator begin() const { return pool.begin(); }
	const_iterator end() const { return pool.end(); }
	// properties / options
	std::string name;
	unsigned int maxSeparated;
	double minEfficacy;
	double minOrtho;
	double maxMemory;
	unsigned int maxSize;
protected:
	std::unordered_set<CutPtr, cutptr_hash, cutptr_equal> pool;
	unsigned int nadded;
	unsigned int nduplicate;
	double mem;
	// helpers
	bool isDuplicate(CutPtr cut) const;
};

} // namespace dominiqs

#endif /* CUTPOOL_H */

