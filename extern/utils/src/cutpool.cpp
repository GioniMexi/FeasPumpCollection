/**
 * @file cutpool.cpp
 * @brief Cut & Cut Pool implementation
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * Copyright 2009 Domenico Salvagnin
 */

#include <algorithm>

#include "utils/cutpool.h"
#include "utils/maths.h"

namespace dominiqs {

/**
 * Class Cut
 */

Cut* Cut::clone() const
{
	return new Cut(*this);
}

double Cut::dynamism() const
{
	double largest = std::numeric_limits<double>::min();
	double smallest = std::numeric_limits<double>::max();
	for (SparseVector::size_type i = 0; i < row.size(); i++)
	{
		double tmp = fabs(row.coef()[i]);
		largest = std::max(largest, tmp);
		smallest = std::min(smallest, tmp);
	}
	return (largest / smallest);
}

void Cut::digest()
{
	sig = 0;
	hash_range(sig, row.idx(), row.idx() + row.size());
	norm = euclidianNorm(row.coef(), row.size());
}

int Cut::memoryUsed() const
{
	return sizeof(Cut) + row.capacity() * (sizeof(int) + sizeof(double));
}

void Cut::save(std::ostream& out) const
{
	out << name << " " << row.size() << " ";
	for (unsigned int k = 0; k < row.size(); k++) out << row.idx()[k] << " " << row.coef()[k] << " ";
	out << sense << " " << rhs << " " << slackType;
}

void Cut::read(std::istream& in)
{
	in >> name;
	unsigned int count = 0;
	in >> count;
	row.clear();
	for (unsigned int k = 0; k < count; k++)
	{
		int idx;
		double coef;
		in >> idx >> coef;
		row.push(idx, coef);
	}
	in >> sense >> rhs >> slackType;
	digest();
}

/**
 * Helpers
 */

struct CutCmpByScore
{
public:
	bool operator()(CutPtr c1, CutPtr c2) const
	{
		return (c1->score() > c2->score());
	}
};

struct CutOrthoFilter
{
public:
	CutOrthoFilter(double minOrtho) : thr(minOrtho) {}
	bool operator()(CutPtr c) const
	{
		return (c->ortho < thr);
	}
protected:
	double thr;
};

struct CutCmpForPurge
{
public:
	bool operator()(CutPtr c1, CutPtr c2) const
	{
		return (c1->purgeScore() < c2->purgeScore());
	}
};

/**
 * Class CutPool
 */

static const int MAX_SEPARATED_DEF = 100;
static const double MIN_EFFICACY_DEF = 0.0001;
static const double MIN_ORTHO_DEF = 0.1;
static const double MAX_MEMORY_DEF = 1E+9;
static const int MAX_SIZE_DEF = 1000000;

CutPool::CutPool()
		: name("CutPool"), maxSeparated(MAX_SEPARATED_DEF), minEfficacy(MIN_EFFICACY_DEF), minOrtho(MIN_ORTHO_DEF),
		maxMemory(MAX_MEMORY_DEF), maxSize(MAX_SIZE_DEF),
		nadded(0), nduplicate(0), mem(0) {}

bool CutPool::push(CutPtr cut)
{
	if (isDuplicate(cut))
	{
		nduplicate++;
		return false;
	}
	nadded++;
	pool.insert(cut);
	mem += (double)(cut->memoryUsed());
	return true;
}

void CutPool::merge(CutPool& other)
{
	for (CutPtr c: other.pool)
	{
		push(c);
	}
	other.clear();
}

bool CutPool::isViolated(const double* x) const
{
	for (CutPtr c: pool)
	{
		c->efficacy = c->violation(&x[0]) / c->norm; // can be negative for slack constraints!
		if (c->efficacy > minEfficacy) return true;
	}
	return false;
}

int CutPool::numViolated(const double* x) const
{
	int ret = 0;
	for (CutPtr c: pool)
	{
		c->efficacy = c->violation(&x[0]) / c->norm; // can be negative for slack constraints!
		if (c->efficacy > minEfficacy) ret++;
	}
	return ret;
}

bool CutPool::select(const std::vector<double>& x, CutList& cuts)
{
	cuts.clear();
	CutList reducedPool;
	// first loop: collect violated cuts and compute initial efficacy/ortho/score
	for (CutPtr c: pool)
	{
		c->efficacy = c->violation(&x[0]) / c->norm; //< can be negative for slack constraints!
		if ((!c->inUse) && (c->efficacy > minEfficacy))
		{
			c->ortho = 1.0;
			reducedPool.push_back(c);
			c->age = 0;
		}
		else c->age++;
	}
	if (isPositive(minOrtho))
	{
		std::vector<double> ar(x.size());
		// second loop: add cut with greatest score and recompoute ortho/score
		while ( (cuts.size() < maxSeparated) && reducedPool.size() )
		{
			// select best cut and remove it from reduced pool
			CutList::iterator bestItr = min_element(reducedPool.begin(), reducedPool.end(), CutCmpByScore());
			CutPtr best = *bestItr;
			cuts.push_back(best);
			*bestItr = reducedPool.back();
			reducedPool.pop_back();
			// update ortho and recompute score
			std::fill(ar.begin(), ar.end(), 0);
			best->row.scatter(&ar[0], ar.size());
			for (CutPtr c: reducedPool)
			{
				double newOrtho = 1.0 - dotProduct(c->row, &ar[0]) / (best->norm * c->norm);
				c->ortho = std::min(newOrtho, c->ortho);
			}
			CutList::iterator newEnd = remove_if(reducedPool.begin(), reducedPool.end(), CutOrthoFilter(minOrtho));
			reducedPool.erase(newEnd, reducedPool.end());
		}
	}
	else
	{
		if (reducedPool.size() > maxSeparated)
		{
			sort(reducedPool.begin(), reducedPool.end(), CutCmpByScore());
			reducedPool.resize(maxSeparated);
		}
		cuts.swap(reducedPool);
	}
	return cuts.size();
}

int CutPool::purge(const std::vector<double>& x)
{
	if ((memoryUsed() < maxMemory) && (size() < maxSize)) return 0;
	CutList::size_type oldSize = pool.size();
	std::vector<CutPtr> tmpPool;
	tmpPool.reserve(oldSize);
	// sort by increasing cut badness
	for (CutPtr c: pool)
	{
		c->efficacy = c->violation(&x[0]) / c->norm;
		tmpPool.push_back(c);
	}
	sort(tmpPool.begin(), tmpPool.end(), CutCmpForPurge());
	// remove last half of the pool
	for (CutList::size_type i = oldSize / 2; i < oldSize; i++)
	{
		pool.erase(tmpPool[i]);
		mem -= (double)(tmpPool[i]->memoryUsed());
	}
	CutList::size_type newSize = pool.size();
	return (oldSize - newSize);
}

bool CutPool::isDuplicate(CutPtr cut) const
{
	DOMINIQS_ASSERT( cut );
	return (pool.count(cut) > 0);
}

void CutPool::readConfig(FileConfig& config, const std::string& root)
{
	maxSeparated = config.get(root + ".maxSeparated", MAX_SEPARATED_DEF);
	minEfficacy = config.get(root + ".minEfficacy", MIN_EFFICACY_DEF);
	minOrtho = config.get(root + ".minOrtho", MIN_ORTHO_DEF);
	maxMemory = config.get(root + ".maxMemory", MAX_MEMORY_DEF);
	maxSize = config.get(root + ".maxSize", MAX_SIZE_DEF);
}

void CutPool::save(std::ostream& out) const
{
	out << name << " " << pool.size() << std::endl;
	for (auto cut: pool)
	{
		cut->save(out);
		out << std::endl;
	}
}

void CutPool::read(std::istream& in)
{
	clear();
	unsigned int count = 0;
	in >> name >> count;
	for (unsigned int k = 0; k < count; k++)
	{
		CutPtr cut(new Cut);
		cut->read(in);
		push(cut);
	}
}

} // namespace dominiqs
