/**
 * @file maths.cpp
 * @brief Mathematical utilities Source
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#include "utils/maths.h"
#include "utils/floats.h"
#include <cstring>

static const int DEF_GATHER_SIZE = 1024;

namespace dominiqs {

SparseVector::SparseVector(const SparseVector& other)
{
	if (other.length)
	{
		int* tmp_idx = nullptr;
		double* tmp_coef = nullptr;
		void* mem = allocateBlock(&tmp_idx, other.length, &tmp_coef, other.length);
		DOMINIQS_ASSERT( mem );
		std::memcpy(tmp_idx, other.m_idx, other.length * sizeof(int));
		std::memcpy(tmp_coef, other.m_coef, other.length * sizeof(double));
		// nothing can throw now
		m_idx = tmp_idx;
		m_coef = tmp_coef;
		length = other.length;
		alloc = other.length;
	}
	else
	{
		m_idx = 0;
		m_coef = 0;
		length = 0;
		alloc = 0;
	}
}

SparseVector& SparseVector::operator=(const SparseVector& other)
{
	if (this != &other)
	{
		SparseVector temp(other);
		swap(temp);
	}
	return *this;
}

void SparseVector::resize(size_type newSize)
{
	if (newSize)
	{
		if (newSize > alloc) reserve(newSize);
		length = newSize;
	}
	else
	{
		freeBlock(m_idx);
		m_idx = 0;
		m_coef = 0;
		length = 0;
		alloc = 0;
	}
}

void SparseVector::reserve(size_type n)
{
	if (n > alloc)
	{
		int* tmp_idx = nullptr;
		double* tmp_coef = nullptr;
		void* mem = allocateBlock(&tmp_idx, n, &tmp_coef, n);
		DOMINIQS_ASSERT( mem );
		std::memcpy(tmp_idx, m_idx, alloc * sizeof(int));
		std::memcpy(tmp_coef, m_coef, alloc * sizeof(double));
		// nothing can throw now
		freeBlock(m_idx);
		m_idx = tmp_idx;
		m_coef = tmp_coef;
		alloc = n;
	}
}

void SparseVector::copy(const int* idx, const double* coef, size_type cnt)
{
	resize(cnt);
	DOMINIQS_ASSERT( length == cnt );
	std::memcpy(m_idx, idx, cnt * sizeof(int));
	std::memcpy(m_coef, coef, cnt * sizeof(double));
}

bool SparseVector::operator==(const SparseVector& rhs) const
{
	if (size() != rhs.size()) return false;
	for (size_type i = 0; i < size(); i++)
	{
		if (m_idx[i] != rhs.m_idx[i]) return false;
		if (different(m_coef[i], rhs.m_coef[i])) return false;
	}
	return true;
}

void SparseVector::gather(const double* in, int n, double eps)
{
	clear();
	int fastN = std::min(n, DEF_GATHER_SIZE);
	reserve(fastN);
	for (int i = 0; i < fastN; ++i) if (isNotNull(in[i], eps)) push_unsafe(i, in[i]);
	for (int i = fastN; i < n; ++i) if (isNotNull(in[i], eps)) push(i, in[i]);
}

void SparseVector::scatter(double* out, int n, bool reset)
{
	if (reset) std::memset((void*)out, 0, n * sizeof(double));
	SparseVector::size_type cnt = size();
	for (SparseVector::size_type i = 0; i < cnt; ++i) out[m_idx[i]] = m_coef[i];
}

void SparseVector::unscatter(double* out)
{
	SparseVector::size_type cnt = size();
	for (SparseVector::size_type i = 0; i < cnt; ++i) out[m_idx[i]] = 0.0;
}

void SparseVector::scale(double lambda)
{
	SparseVector::size_type cnt = size();
	for (SparseVector::size_type i = 0; i < cnt; ++i) m_coef[i] *= lambda;
}

void SparseVector::negate()
{
	SparseVector::size_type cnt = size();
	for (SparseVector::size_type i = 0; i < cnt; ++i) m_coef[i] = -m_coef[i];
}

double Constraint::violation(const double* x) const
{
	double slack = rhs - dotProduct(row, x);
	if (sense == 'L') return -slack;
	if (sense == 'G') return slack;
	if (sense == 'R') return std::max(-slack, slack-range);
	return fabs(slack);
}

void accumulate(double* v, const double* w, int n, double lambda)
{
	for (int i = 0; i < n; ++i) v[i] += (lambda * w[i]);
}

void accumulate(double* v, const int* wIdx, const double* wCoef, int n, double lambda)
{
	for (int i = 0; i < n; ++i) v[wIdx[i]] += (lambda * wCoef[i]);
}

void scale(double* v, int n, double lambda)
{
	for (int i = 0; i < n; ++i) v[i] *= lambda;
}

double dotProduct(const double* x, const double* y, int n)
{
	double ans = 0.0;
	int i;
	if ( n >= 8 )
	{
		for ( i = 0; i < ( n >> 3 ); ++i, x += 8, y += 8 )
			ans += x[0] * y[0] + x[1] * y[1] +
				x[2] * y[2] + x[3] * y[3] +
				x[4] * y[4] + x[5] * y[5] +
				x[6] * y[6] + x[7] * y[7];
		n -= i << 3;
	}
	for ( i = 0; i < n; ++i ) ans += (x[i] * y[i]);
	return ans;
}

double dotProduct(const int* idx, const double* x, int n, const double* y)
{
	double ans = 0.0;
	for (int i = 0; i < n; i++) ans += (x[i] * y[idx[i]]);
	return ans;
}

bool disjoint(const double* x, const double* y, int n)
{
	int cnt = n;
	cnt = cnt >> 1;
	for (int i = cnt; i > 0 ; --i)
	{
		double p1 = fabs(x[0] * y[0]);
		double p2 = fabs(x[1] * y[1]);
		if (isPositive(p1 + p2)) return false;
		x += 2;
		y += 2;
	}
	if (n % 2) return (isPositive(fabs(x[0] * y[0])));
	return true;
}

double euclidianNorm(const double* v, int n)
{
	double ans = 0.0;
	for (int i = 0; i < n; i++) ans += (v[i] * v[i]);
	return sqrt(ans);
}

double euclidianDistance(const double* a1, const double* a2, int n)
{
	double ans = 0.0;
	for (int i = 0; i < n; i++) ans += ((a1[i] - a2[i]) * (a1[i] - a2[i]));
	return sqrt(ans);
}

int lexComp(const double* s1, const double* s2, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (lessThan(s1[i], s2[i])) return -1;
		if (greaterThan(s1[i], s2[i])) return 1;
	}
	return 0;
}

} // namespace dominiqs
