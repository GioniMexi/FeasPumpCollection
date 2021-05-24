/**
 * @file index_partition.h
 * @brief Indexed partition data structure
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef INDEX_PARTITION_H
#define INDEX_PARTITION_H

#include "asserter.h"
#include "maths.h"
#include "floats.h"
#include <numeric>
#include <boost/iterator/iterator_facade.hpp>

namespace dominiqs
{

/**
 * @brief Common Iterator Class for IndexPartition and IndexSet
 */

class IndexPartition;
class IndexSet;

template<typename IndexType, typename ValueType>
class IndexIteratorImpl : public boost::iterator_facade<IndexIteratorImpl<IndexType,ValueType>, ValueType, boost::bidirectional_traversal_tag>
{
public:
	IndexIteratorImpl() {}
	explicit IndexIteratorImpl(IndexType* p, int pos) : part(p), index(pos) {}
	void debug() const { std::cout << part << "[" << index << "]"; }
private:
	// interface
	friend class boost::iterator_core_access;
	void increment() { index = part->flink[index]; }
	void decrement() { index = part->blink[index]; }
	template <class OtherValue>
	bool equal(OtherValue const& other) const
	{
		return (index == other.index);
	}
	ValueType& dereference() const { return index; }
	// data
	IndexType* part = nullptr;
	int index = 0;
};


/**
 * @brief Efficient data structure for keeping linked sets of indices
 *
 * The data structure stores numIndices indices from 0 to numIndices - 1 
 * in up to numSets sets of indices 0 to numSets - 1
 * A given index can be in at most one set at a given time.
 */

class IndexPartition
{
public:
	/* Creates an index partition with indices up to @param n and sets up to @param ns */
	IndexPartition(int n = 0, int ns = 0) : numIndices(n), numSets(ns), flink(numIndices + numSets), blink(numIndices + numSets)
	{
		std::iota(flink.begin(), flink.end(), 0);
		std::iota(blink.begin(), blink.end(), 0);
	}
	/* Resizes the partition to indices up to @param n and to sets up to @param ns (previous content is reset!) */
	inline void resize(int n, int ns)
	{
		numIndices = n;
		numSets = ns;
		flink.resize(numIndices + numSets);
		std::iota(flink.begin(), flink.end(), 0);
		blink.resize(numIndices + numSets);
		std::iota(blink.begin(), blink.end(), 0);
	}
	/* Reset the content of the set, without resizing */
	inline void clear()
	{
		std::iota(flink.begin(), flink.end(), 0);
		std::iota(blink.begin(), blink.end(), 0);
	}
	/* add index @param index to set @param set (does NOT check whether the item was already in other set!!!) */
	inline void add(int index, int set)
	{
		DOMINIQS_ASSERT( index >= 0 && index < numIndices );
		DOMINIQS_ASSERT( set >= 0 && set < numSets );
		DOMINIQS_ASSERT( isInNoSet(index) );
		int k = flink[numIndices + set];
		flink[numIndices + set] = index;
		flink[index] = k;
		blink[k] = index;
		blink[index] = numIndices + set;
	}
	/* remove index @param index from its current set (safe even if the index was already in no set) */
	inline void remove(int index)
	{
		DOMINIQS_ASSERT( index >= 0 && index < numIndices );
		flink[blink[index]] = flink[index];
		blink[flink[index]] = blink[index];
		flink[index] = index;
		blink[index] = index;
	}
	inline bool isEmpty(int set) const
	{
		return (flink[numIndices + set] == (numIndices + set));
	}
	inline int top(int set) const
	{
		DOMINIQS_ASSERT( !isEmpty(set) );
		return flink[numIndices + set];
	}
	inline bool isInNoSet(int index) const
	{
		return (flink[index] == index);
	}
	inline int getNumIndices() const { return numIndices; }
	inline int getNumSets() const { return numSets; }
	/* Iterators */
	typedef IndexIteratorImpl<IndexPartition, int const> iterator;
	typedef IndexIteratorImpl<IndexPartition const, int const> const_iterator;
	friend class IndexIteratorImpl<IndexPartition, int const>;
	friend class IndexIteratorImpl<IndexPartition const, int const>;
	iterator begin(int set) { return iterator(this, flink[numIndices + set]); }
	iterator end(int set) { return iterator(this, numIndices + set); }
	const_iterator begin(int set) const { return const_iterator(this, flink[numIndices + set]); }
	const_iterator end(int set) const { return const_iterator(this, numIndices + set); }
protected:
	// data
	int numIndices;
	int numSets;
	std::vector<int> flink;
	std::vector<int> blink;
};

/**
 * @brief Efficient data structure for keeping a linked set of indices
 *
 * This is a special case of @class IndexPartition, in which there is only one set
 */

class IndexSet
{
public:
	/* Creates an index set with indices up to @param n */
	IndexSet(int n = 0) : numIndices(n), flink(numIndices + 1), blink(numIndices + 1)
	{
		std::iota(flink.begin(), flink.end(), 0);
		std::iota(blink.begin(), blink.end(), 0);
	}
	/* Resizes the partition to indices up to @param n (previous content is reset!) */
	inline void resize(int n)
	{
		numIndices = n;
		flink.resize(numIndices + 1);
		std::iota(flink.begin(), flink.end(), 0);
		blink.resize(numIndices + 1);
		std::iota(blink.begin(), blink.end(), 0);
		_count = 0;
	}
	/* Reset the content of the set, without resizing */
	inline void clear()
	{
		iterator itr = begin();
		iterator last = end();
		while (itr != last)
		{
			int elem = *itr++;
			remove(elem);
		}
	}
	/* add index @param index (does NOT check whether the item was already in other set!!!) */
	inline void add(int index)
	{
		push_back(index);
	}
	/* add index @param index to the back (does NOT check whether the item was already in other set!!!) */
	inline void push_back(int index)
	{
		DOMINIQS_ASSERT( index >= 0 && index < numIndices );
		DOMINIQS_ASSERT( !has(index) );
		int k = blink[numIndices];
		flink[k] = index;
		flink[index] = numIndices;
		blink[index] = k;
		blink[numIndices] = index;
		_count++;
	}
	/* add index @param index to the front (does NOT check whether the item was already in other set!!!) */
	inline void push_front(int index)
	{
		DOMINIQS_ASSERT( index >= 0 && index < numIndices );
		DOMINIQS_ASSERT( !has(index) );
		int k = flink[numIndices];
		flink[numIndices] = index;
		flink[index] = k;
		blink[k] = index;
		blink[index] = numIndices;
		_count++;
	}
	/* remove index @param index (safe even if the index was already not in) */
	inline void remove(int index)
	{
		DOMINIQS_ASSERT( index >= 0 && index < numIndices );
		_count -= (has(index));
		flink[blink[index]] = flink[index];
		blink[flink[index]] = blink[index];
		flink[index] = index;
		blink[index] = index;
	}
	inline bool isEmpty() const
	{
		return (flink[numIndices] == numIndices);
	}
	inline int top() const
	{
		DOMINIQS_ASSERT( !isEmpty() );
		return flink[numIndices];
	}
	inline bool has(int index) const
	{
		return (flink[index] != index);
	}
	inline int getNumIndices() const { return numIndices; }
	inline int count() const { return _count; }
	/* Iterators */
	typedef IndexIteratorImpl<IndexSet, int const> iterator;
	typedef IndexIteratorImpl<IndexSet const, int const> const_iterator;
	friend class IndexIteratorImpl<IndexSet, int const>;
	friend class IndexIteratorImpl<IndexSet const, int const>;
	iterator begin() { return iterator(this, flink[numIndices]); }
	iterator end() { return iterator(this, numIndices); }
	const_iterator begin() const { return const_iterator(this, flink[numIndices]); }
	const_iterator end() const { return const_iterator(this, numIndices); }
protected:
	// data
	int numIndices;
	std::vector<int> flink;
	std::vector<int> blink;
	int _count = 0;
};


/**
 * @brief Efficient data structure for keeping a sparse vector with a linked set of indices
 */

class IndexVector
{
public:
	// typedefs
	typedef typename std::vector<double>::value_type value_type;
	typedef typename std::vector<double>::pointer pointer;
	typedef typename std::vector<double>::const_pointer const_pointer;
	typedef typename std::vector<double>::reference reference;
	typedef typename std::vector<double>::const_reference const_reference;
	typedef typename std::vector<double>::size_type size_type;
	IndexVector(size_type n = 0) : data(n), support(n), changed(n) {}
	explicit IndexVector(const std::vector<double>& other) : data(other), support(other.size()), changed(other.size())
	{
		for (size_type j = 0; j < size(); j++)
		{
			if (isNotNull(data[j])) support.add(j);
		}
	}
	inline void write(std::vector<double>& other) const { other = data; }
	inline void set(size_type i, double val)
	{
		double oldval = data[i];
		data[i] = val;
		if (dominiqs::different(oldval, val))
		{
			changed.add(i);
			if (dominiqs::isNull(val)) support.remove(i);
			if (dominiqs::isNotNull(val) && !support.has(i)) support.add(i);
		}
	}
	inline const_reference operator[](size_type i) const { return data[i]; }
	inline size_type size() const { return data.size(); }
	inline bool empty() const { return data.empty(); }
	inline void resize(size_type n)
	{
		data.resize(n);
		std::fill(data.begin(), data.end(), 0.0);
		support.resize(n);
		changed.resize(n);
	}
	inline void clear()
	{
		data.clear();
		support.resize(0);
		changed.resize(0);
	}
	// dense iterators
	typedef typename std::vector<double>::iterator iterator;
	typedef typename std::vector<double>::const_iterator const_iterator;
	inline iterator begin() { return data.begin(); }
	inline iterator end() { return data.end(); }
	inline const_iterator begin() const { return data.begin(); }
	inline const_iterator end() const { return data.end(); }
	// support iterators
	typedef IndexSet::iterator index_iterator;
	typedef IndexSet::const_iterator index_const_iterator;
	inline index_iterator beginS() { return support.begin(); }
	inline index_iterator endS() { return support.end(); }
	inline index_const_iterator beginS() const { return support.begin(); }
	inline index_const_iterator endS() const { return support.end(); }
	// changed iterators
	inline index_iterator beginC() { return changed.begin(); }
	inline index_iterator endC() { return changed.end(); }
	inline index_const_iterator beginC() const { return changed.begin(); }
	inline index_const_iterator endC() const { return changed.end(); }
	// reset changed
	inline void resetChanged() { changed.clear(); }
protected:
	std::vector<double> data;
	IndexSet support;
	IndexSet changed;
};

} // namespace

#endif /** INDEX_PARTITION_H */
