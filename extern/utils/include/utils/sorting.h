/**
 * @file sorting.h
 * @brief Sorting utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2011-2014
 */

#ifndef SORTING_H
#define SORTING_H

#include <vector>
#include <numeric>
#include <functional>

namespace dominiqs {

/**
 * @brief Sort a given vector in place using a modified Shell Sort algorithm O(n^(3/2))
 *
 * @tparam T vector value type
 * @tparam StrictWeakOrdering comparison operator
 *
 * @param v input vector (passed as C array)
 * @param n vector size
 * @param comp operator instance
 */

template<typename T, typename StrictWeakOrdering>
void shellSort(T* v, size_t n, StrictWeakOrdering comp)
{
	if (n == 0) return;
	size_t gap = n / 2;
	if (gap % 2 == 0) gap++;
	while (true)
	{
		for (size_t k = gap; k < n; k++)
		{
			T tmp = v[k];
			size_t j = k;
			while (j >= gap && comp(tmp, v[j - gap]))
			{
				v[j] = v[j - gap];
				j -= gap;
			}
			v[j] = tmp;
		}
		if (gap == 1) break;
		else
		{
			gap = gap / 2;
			if (gap % 2 == 0) gap++;
		}
	}
}

/**
 * @brief Syntactic sugar for using operator< as default comparison operator
 *
 * @tparam T vector type
 *
 * @param v input vector
 */

template<typename T>
void shellSort(T* v, size_t n)
{
	shellSort(v, n, std::less<T>());
}

/**
 * @brief Syntactic sugar for std::vectors
 *
 * @tparam T vector type
 * @tparam StrictWeakOrdering comparison operator
 *
 * @param v input vector
 * @param comp operator instance
 */

template<typename T, typename StrictWeakOrdering>
void shellSort(T& v, StrictWeakOrdering comp)
{
	if (v.size() > 1) shellSort(&v[0], v.size(), comp);
}

/**
 * @brief Syntactic sugar for std::vectors and operator< as comparison operator
 *
 * @tparam T vector type
 *
 * @param v input vector
 */

template<typename T>
void shellSort(T& v)
{
	if (v.size() > 1) shellSort(&v[0], v.size(), std::less<typename T::value_type>());
}

/**
 * @brief Sort a given vector using a modified Shell Sort algorithm O(n^(3/2))
 *
 * The algorithm does not rearrange the vector items, rather it computes
 * the permutation that does the sorting
 *
 * @tparam T value type of vector
 * @tparam P value type of permutation vector
 * @tparam StrictWeakOrdering comparison operator
 *
 * @param v input vector (passed as C array)
 * @param perm permutation vector (passed as C array of the appropriate size!)
 * @param n vectors size (both input and permutation)
 * @param comp operator instance
 * @param inputPerm if true perm stores a starting permutation of v
 * if false perm is initialized with the identity permutation [0, 1, ..., |v| - 1]
 */

template<typename T, typename P, typename StrictWeakOrdering>
void permShellSort(const T* v, P* perm, size_t n, StrictWeakOrdering comp, bool inputPerm = false)
{
	if (n == 0) return;
	if (!inputPerm) std::iota(perm, perm + n, 0);
	size_t gap = n / 2;
	if (gap % 2 == 0) gap++;
	while (true)
	{
		for (size_t k = gap; k < n; ++k)
		{
			P tmp = perm[k];
			size_t j = k;
			while (j >= gap && comp(v[tmp], v[perm[j - gap]]))
			{
				perm[j] = perm[j - gap];
				j -= gap;
			}
			perm[j] = tmp;
		}
		if (gap == 1) break;
		else
		{
			gap = gap / 2;
			if (gap % 2 == 0) gap++;
		}
	}
}

/**
 * @brief Sort a given vector in place using the insertion sort algorithm O(n^2)
 * Insertion sort is asymptotically worst than shell short, but it is stable.
 *
 * @tparam T vector value type
 * @tparam StrictWeakOrdering comparison operator
 *
 * @param v input vector (passed as C array)
 * @param n vector size
 * @param comp operator instance
 */

template<typename T, typename StrictWeakOrdering>
void insertionSort(T* v, size_t n, StrictWeakOrdering comp)
{
	if (n == 0) return;
	for (size_t i = 1; i < n; ++i)
	{
		T k = v[i];
		size_t j = i;
		while ((j > 0) && !comp(v[j-1], k))
		{
			v[j] = v[j-1];
			--j;
		}
		v[j] = k;
	}
}

/**
 * @brief Permute the elements of a vector according to a given permutation
 *
 * @tparam T value type of vector
 * @tparam P value type of permutation vector
 *
 * @param v input/output vector (passed as C array)
 * @param perm permutation vector (passed as C array of the appropriate size!)
 * @param n vectors size (both vector and permutation)
 */
template<typename T, typename P>
void permute(T* v, const P* perm, size_t n)
{
	std::vector<bool> done(n, false);
	for (size_t i = 0; i < n; ++i)
	{
		if (!done[i])
		{
			T tmp = v[i];
			for (size_t j = i;;)
			{
				done[j] = true;

				if (perm[j] != i)
				{
					v[j] = v[perm[j]];
					j = perm[j];
				}
				else
				{
					v[j] = tmp;
					break;
				}
			}
		}
	}
}

} // namespace dominiqs

#endif /* SORTING_H */
