/**
 * @file trail.h
 * @brief Trail to track and undo changes on a std::vector
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * Copyright 2017 Domenico Salvagnin
 */

#ifndef TRAIL_H
#define TRAIL_H

#include <vector>
#include <limits>

namespace dominiqs
{

/**
 * @brief Trail to track and undo changes on a std::vector.
 *
 * @tparam T value type of the associated std::vector
 */

template<typename T>
class trail
{
public:
	/**
	 * Create an empty trail
	 */
	trail(std::vector<T>& v) : data(v), posmap(data.size(), npos) {}

	/* Disable copy constructor & assignment operator */
	trail(const trail& that) = delete;
	trail& operator=(const trail& other) = delete;

	/* Get the number of changes made so far */
	inline size_t nChanged() const { return changes.size(); }

	/* Commit the changes to the vector (they won't be undone) */
	inline void clear()
	{
		for (const Change& c: changes) posmap[c.pos] = npos;
		changes.clear();
	}

	/* Was the value in position pos changed? */
	bool hasChanged(size_t pos) const { return (posmap[pos] != npos); }

	/* Get the original value at position pos */
	T getOriginal(size_t pos) const
	{
		if (hasChanged(pos)) return changes[posmap[pos]].oldval;
		return data[pos];
	}

	/* Set method to change a value int the vector (and keep track of the old value) */
	void set(size_t pos, T value)
	{
		if (posmap[pos] == npos)
		{
			// the value was never changed before
			posmap[pos] = changes.size();
			changes.emplace_back(pos, data[pos]);
			data[pos] = value;
		}
		else
		{
			// value was changed already: we do not need to keep track
			data[pos] = value;
		}
	}

	/* Undo all changes applied to the vector */
	void restore()
	{
		for (const Change& c: changes)
		{
			data[c.pos] = c.oldval;
			posmap[c.pos] = npos;
		}
		changes.clear();
	}

	/* Data structure for changed item */
	struct Change
	{
	public:
		Change(size_t p, T v) : pos(p), oldval(v) {}
		size_t pos;
		T oldval;
	};
	const std::vector<Change>& getChanges() const { return changes; }

private:
	static constexpr size_t npos = std::numeric_limits<std::size_t>::max(); //< marker for unchanged positions
	std::vector<T>& data; //< vector we are tracking changes of
	std::vector<size_t> posmap; //< flag positions that were changed
	std::vector<Change> changes; //< trail of changes
};

/* Definition (as opposed to declaration) for trail::npos */
template<typename T> constexpr size_t trail<T>::npos;

} // namespace

#endif /** TRAIL_H */
