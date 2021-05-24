/**
 * @file pq.h
 * @brief Priority Queue for a set of integers
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2014 Domenico Salvagnin
 */

#ifndef PQ_H
#define PQ_H

#include <utils/asserter.h>
#include <utils/floats.h>
#include <vector>
#include <iostream>

template<typename ScoreType>
class PriorityQueue
{
public:
	/**
	 * Construct a priority queue for the integer form 0 to @param _n - 1
	 */
	PriorityQueue(int _n) : n(_n), cnt(0), data(n), prior(n), position(n, -1) {}
	/**
	 * Return the integer with minimal priority (throws exception if empty)
	 */
	int top() const
	{
		DOMINIQS_ASSERT( !isEmpty() );
		return data[0];
	}
	/**
	 * Check whether the queue is empty
	 */
	bool isEmpty() const { return (cnt == 0); }
	/**
	 * Checks whether an integer @param j is in the queue
	 */
	bool has(int j) const { return (position[j] >= 0); }
	/**
	 * Clear content
	 */
	void clear()
	{
		std::fill(position.begin(), position.begin() + cnt, -1);
		cnt = 0;
	}
	/**
	 * Restore heap structure
	 */
	void heapify()
	{
		int start = (cnt - 2) / 2;
		while (start >= 0)
		{
			int j = data[start];
			position[j] = -1;
			ScoreType p = prior[j];
			int gap = start;
			siftDown(gap, p);
			data[gap] = j;
			position[j] = gap;
			start--;
		}
	}
	/**
	 * Insert integer @param j into the queue with a priority @param p
	 */
	void push(int j, ScoreType p, bool maintainHeap = true)
	{
		// std::cout << " >Adding " << j << std::endl;
		DOMINIQS_ASSERT( (j >= 0) && (j < n) );
		DOMINIQS_ASSERT( position[j] == -1 );
		prior[j] = p;
		// put gap at last position
		int gap = cnt++;
		if (maintainHeap)
		{
			// move gap down to the proper position
			siftUp(gap, p);
		}
		// fill gap with new element
		data[gap] = j;
		position[j] = gap;
		// DOMINIQS_ASSERT( checkConsistency() );
	}
	/**
	 * Removes the integer with minimal priority
	 */
	void pop()
	{
		DOMINIQS_ASSERT( !isEmpty() );
		removeAt(0, true);
	}
	/**
	 * Removes integer @param j from the queue
	 */
	void remove(int j, bool maintainHeap = true)
	{
		DOMINIQS_ASSERT( position[j] >= 0 );
		removeAt(position[j], maintainHeap);
		DOMINIQS_ASSERT( position[j] == -1 );
	}
	/**
	 * Changes the score of an integer @param j already in the queue to @param p
	 */
	void change(int j, ScoreType p, bool maintainHeap = true)
	{
		DOMINIQS_ASSERT( position[j] >= 0 );
		int gap = position[j];
		if (maintainHeap)
		{
			ScoreType oldp = prior[gap];
			// if the priority didn't change, return immediately
			if (dominiqs::equal(oldp, p)) return;
			if (maintainHeap)
			{
				// percolate up or down depending on new priority value
				if (p < oldp) siftUp(gap, p);
				else siftDown(gap, p);
			}
			// fill gap with last element
			data[gap] = j;
			position[j] = gap;
		}
		prior[j] = p;
	}
	/**
	 * Debug
	 */
	bool checkConsistency() const
	{
		// check heap property
		for (int k = 0; k < cnt; k++)
		{
			int left = 2*k+1;
			int right = 2*k+2;
			if ((left < cnt) && (prior[data[left]] < prior[data[k]])) return false;
			if ((right < cnt) && (prior[data[right]] < prior[data[k]])) return false;
		}
		// check positions
		for (int k = 0; k < cnt; k++)
		{
			int i = data[k];
			if (position[i] != k) return false;
		}
		// check positions
		for (int i = 0; i < n; i++)
		{
			int pos = position[i];
			if ((pos >= 0) && (pos < cnt) && (data[pos] != i)) return false;
			if (pos >= cnt) return false;
		}
		return true;
	}
	void print() const
	{
		for (int k = 0; k < cnt; k++)
		{
			std::cout << "\t" << k << "\t" << data[k] << "\t" << prior[data[k]] << std::endl;
		}
	}
protected:
	int n; //< size of arrays heap and position
	int cnt; //< number of elements in the queue
	std::vector<int> data;
	std::vector<ScoreType> prior;
	std::vector<int> position;
	/**
	 * Moves the gap up from the current position @param gap to the proper place for a priority of value @param p
	 * The final position of the gap is stored in @param gap
	 */
	inline void siftUp(int& gap, ScoreType p)
	{
		while (gap > 0)
		{
			int parent = (gap - 1) / 2;
			if (p < prior[data[parent]])
			{
				data[gap] = data[parent];
				position[data[gap]] = gap;
				gap = parent;
			}
			else break;
		}
	}
	/**
	 * Moves the gap down from the current position @param gap
	 * The final position of the gap is stored in @param gap
	 */
	inline void siftDown(int& gap, ScoreType p)
	{
		int newgap;
		int left;
		int right;
		while (true)
		{
			left = 2 * gap + 1;
			right = left + 1;
			if (right < cnt)
			{
				if (prior[data[left]] < prior[data[right]]) newgap = left;
				else newgap = right;
				if (p < prior[data[newgap]]) break;
				else
				{
					data[gap] = data[newgap];
					position[data[gap]] = gap;
					gap = newgap;
				}
			}
			else if (right == cnt)
			{
				newgap = left;
				if (p < prior[data[newgap]]) break;
				data[gap] = data[newgap];
				position[data[gap]] = gap;
				gap = newgap;
			}
			else break;
		}
	}
	/**
	 * Removes the integer in position gap
	 */
	void removeAt(int gap, bool maintainHeap)
	{
		int last = --cnt;
		position[data[gap]] = -1;
		if (gap == last) return;
		ScoreType gapP = prior[data[gap]];
		// get last element priority
		ScoreType lastP = prior[data[last]];
		if (maintainHeap)
		{
			// we have a gap in position pos,
			// that we will fill with an element of priority lastP
			// Depending this value, we may need to percolate up or down!!!
			if (lastP < gapP) siftUp(gap, lastP);
			else siftDown(gap, lastP);
		}
		// fill gap with last element
		data[gap] = data[last];
		position[data[gap]] = gap;
		// DOMINIQS_ASSERT( checkConsistency() );
	}
};

#endif /* PQ_H */
