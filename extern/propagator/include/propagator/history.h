/**
 * @file history.h
 * @brief State Managemente/Undo System Base Classes
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef HISTORY_H
#define HISTORY_H

#include <memory>

/**
 * @brief State Management Base Class
 *
 * This class provides the simple interface dump & restore
 * used for state managament.
 */

class State
{
public:
	/**
	 * @brief copy some state inside the object
	 *
	 * Derived classes are responsible for providing a link to the object we are
	 * storing the state for
	 */
	virtual void dump() = 0;
	/**
	 * @brief restore the state of the object
	 *
	 * Derived classes are responsible for providing a link to the object we are
	 * restoring the state to
	 */
	virtual void restore() =0;
	// virtual destructor
	virtual ~State() {}
};

typedef std::shared_ptr<State> StatePtr;

#endif /* HISTORY_H */
