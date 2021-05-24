/**
 * @file maths.h
 * @brief Mathematical utilities Header
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef MATHS_H
#define MATHS_H

#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <functional>

#include "asserter.h"
#include "floats.h"
#include "machine_utils.h"


namespace dominiqs {

/**
 * Sparse Vector
 * Stores a vector in sparse form
 * The internal vectors have 16-bit alignment (SSE2 compatible)
 */

class SparseVector
{
public:
	typedef unsigned int size_type;
	// constructors
	SparseVector() {}
	SparseVector(const SparseVector& other);
	// assignment operator
	SparseVector& operator=(const SparseVector& other);
	// destructor
	inline ~SparseVector() throw()
	{
		freeBlock(m_idx);
	}
	inline void swap(SparseVector& other) throw()
	{
		std::swap(m_idx, other.m_idx);
		std::swap(m_coef, other.m_coef);
		std::swap(length, other.length);
		std::swap(alloc, other.alloc);
	}
	// size
	inline size_type size() const { return length; }
	inline size_type capacity() const { return alloc; }
	void resize(size_type newSize);
	inline bool empty() const { return (length == 0); }
	void reserve(size_type n);
	inline void clear() { length = 0; }
	// push/pop
	inline void push(int i, double v)
	{
		if (length == alloc) reserve(std::max(2 * alloc, (size_type)8));
		push_unsafe(i, v);
	}
	inline void push_unsafe(int i, double v)
	{
		m_idx[length] = i;
		m_coef[length] = v;
		++length;
	}
	inline void pop()
	{
		if (length) --length;
	}
	void copy(const int* idx, const double* coef, size_type cnt);
	// get data
	//@{
	int* idx() { return m_idx; }
	const int* idx() const { return m_idx; }
	double* coef() { return m_coef; }
	const double* coef() const { return m_coef; }
	//@}
	/**
	* Compare two sparse vectors for equality
	* CAUTION: the procedure assumes the indexes to be in increasing order, i.e., sorted!!!
	*/
	bool operator==(const SparseVector& rhs) const;
	/** conversions to/from dense vectors */
	/** Shrink a dense vector into a sparse one */
	void gather(const double* in, int n, double eps = defaultEPS);
	/** Expand a sparse vector into a dense one */
	void scatter(double* out, int n, bool reset = false);
	/** Zero out the coefficients of a dense vector corresponding to the support of a sparse vector */
	void unscatter(double* out);
	/** Scale sparse vector by a given multiplier */
	void scale(double lambda = 1.0);
	/** Invert signs (e.g., multiply by -1) */
	void negate();
protected:
	int* m_idx = nullptr; //< indexes
	double* m_coef = nullptr; //< coefficients
	size_type length = 0;
	size_type alloc = 0;
};

typedef std::shared_ptr<SparseVector> SparseVectorPtr;

typedef std::shared_ptr<std::vector<double>> FloatVectorPtr;

/**
 * Linear Constraint
 */

class Constraint
{
public:
	/** @name Data */
	//@{
	std::string name; //< constraint name
	SparseVector row; /**< coefficient row \f$ a_i \f$ */
	double rhs; /**< right hand side \f$ b_i \f$ */
	double range = 0.0; /**< range value for ranged row: linear expression in [rhs-range,rhs] */
	char sense; /**< constraint sense */
	//@}
	Constraint* clone() const { return new Constraint(*this); }
	/**
	 * Check if this constraint is satisfied by assignment x, with tolerance eps
	 */
	bool satisfiedBy(const double* x, double eps = defaultEPS) const;
	/**
	 * Compute the violation of constraint by assignment x
	 * A positive value means a constraint violation,
	 * while a negative one means constraint is slack
	 */
	double violation(const double* x) const;
	/**
	 * Check if a constraint is slack w.r.t. assignment x, with tolerance eps
	 * @return true if constraint is slack, false otherwise
	 */
	bool isSlack(const double* x, double eps = defaultEPS) const;
};

inline bool Constraint::satisfiedBy(const double* x, double eps) const
{
	return (!isPositive(violation(x), eps));
}

inline bool Constraint::isSlack(const double* x, double eps) const
{
	return isNegative(violation(x), eps);
}

typedef std::shared_ptr<Constraint> ConstraintPtr;

/**
 * Performe the operation: v <- v + lambda w
 * where both v and w are dense vectors and lambda is a scalar
 */

void accumulate(double* v, const double* w, int n, double lambda = 1.0);

/**
 * Performe the operation: v <- v + lambda w
 * where v is a dense vector, w is sparse and lambda is a scalar
 */

void accumulate(double* v, const int* wIdx, const double* wCoef, int n, double lambda = 1.0);

inline void accumulate(double* v, const SparseVector& w, double lambda = 1.0)
{
	accumulate(v, w.idx(), w.coef(), w.size(), lambda);
}

/**
 * Scale a dense vector v multiplying it by lambda: v <- lambda v
 */

void scale(double* v, int n, double lambda = 1.0);

/**
 * Dot Product between dense vectors (manual loop unrolling)
 */

double dotProduct(const double* a1, const double* a2, int n);

/**
 * Dot Product between a sparse and a dense vector
 */

double dotProduct(const int* idx1, const double* a1, int n, const double* a2);

inline double dotProduct(const SparseVector& a1, const double* a2)
{
	return dotProduct(a1.idx(), a1.coef(), a1.size(), a2);
}

/**
 * Checks if two dense vectors have disjoint support
 */

bool disjoint(const double* a1, const double* a2, int n);

/**
 * Euclidian norm on a dense vector
 */

double euclidianNorm(const double* v, int n);

/**
 * Euclidian distance between two dense vectors
 */

double euclidianDistance(const double* a1, const double* a2, int n);

/**
 * Lexicographically compare two array of doubles
 * @param s1 first array
 * @param s2 second array
 * @param n arrays size
 * @return comparison result {-1, 0, 1} if {<, =, >} respectively
 */

int lexComp(const double* s1, const double* s2, int n);


/**
 * Sparse Matrix
 * Stores a vector in sparse form (either col or row major)
 */
class SparseMatrix
{
public:
	std::vector<int> matbeg;
	std::vector<int> matind;
	std::vector<double> matval;
	int k = 0; // size of vector matbeg
	int nnz = 0; // size of vectors matind and matval
};


/**
 * Unary predicate that incrementally compute the variance of a list of numbers
 * The results can be obtained with result()
 * The flag fromSample decides if the variance was from the entire population
 * or a subset of observation, i.e. if fromSample is true the sum of squares is
 * divided by n - 1, otherwise by n
 */

class IncrementalVariance : public std::unary_function<double, void>
{
public:
	IncrementalVariance() : cnt(0), mean(0.0), sumsq(0.0) {}
	void operator() (double x)
	{
		cnt++;
		double delta = x - mean;
		mean += delta / cnt;
		sumsq += delta * (x - mean);
	}
	int count() const { return cnt; }
	double result(bool fromSample = false) const { return (cnt > 1) ? (sumsq / (cnt - int(fromSample))) : 0.0; }
protected:
	int cnt;
	double mean;
	double sumsq;
};

/**
 * Return the arithmetic mean of a list of numbers
 */

template<typename ForwardIterator, typename T = double>
T mean(ForwardIterator first, ForwardIterator last)
{
	T sum = 0.0;
	int cnt = 0;
	std::for_each(first, last, [&](T value) { sum += value; ++cnt; });
	return (cnt > 0) ? (sum / cnt) : sum;
}

/**
 * Return the geometric mean of a list of positive numbers
 */

template<typename ForwardIterator, typename T = double>
T geomMean(ForwardIterator first, ForwardIterator last)
{
	T sum = 0.0;
	int cnt = 0;
	std::for_each(first, last, [&](T value) { sum += std::log(value); ++cnt; });
	return (cnt > 0) ? std::exp(sum / cnt) : sum;
}

/**
 * Return the variance of a list of numbers
 * The flag fromSample decides if the variance was from the entire population
 * or a subset of observation, i.e. if fromSample is true the sum of squares is
 * divided by n - 1, otherwise by n
 */

template<typename ForwardIterator, typename T = double>
T variance(ForwardIterator first, ForwardIterator last, bool fromSample = false)
{
	T m = mean(first, last);
	int cnt = 0;
	T sum = 0.0;
	std::for_each(first, last, [&](T value) { sum += std::pow(value - m, 2.0); ++cnt; });
	return (cnt > 1) ? (sum / (cnt - int(fromSample))) : sum;
}

/**
 * Return the standard deviation of a list of numbers
 * The flag fromSample decides if the variance was from the entire population
 * or a subset of observation, i.e. if fromSample is true the sum of squares is
 * divided by n - 1, otherwise by n
 */

template<typename ForwardIterator, typename T = double>
T stdev(ForwardIterator first, ForwardIterator last, bool fromSample = false)
{
	return std::sqrt(variance(first, last, fromSample));
}

/**
 * Compute the hash of a set of objects incrementally
 */
template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
	std::hash<T> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

/**
 * Compute the hash of the elements in the range [first,last)
 */
template <class ForwardIterator>
inline void hash_range(std::size_t& seed, ForwardIterator first, ForwardIterator last)
{
	while (first != last) hash_combine(seed, *first++);
}


} // namespace dominiqs

#endif /* MATHS_H */
