/**
 * @file singleton.h
 * @brief Singleton Pattern
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2011
 */

#ifndef SINGLETON_H
#define SINGLETON_H

namespace dominiqs {

/**
 * Class implementing the singleton pattern for type T
 * @param T type of object
 * Note that the implementation below is NOT thread-safe
 */

template<class T>
class SingletonHolder
{
public:
	static T& getInstance()
	{
		static T theInstance;
		return theInstance;
	}
};

} // namespace dominiqs

#endif /* SINGLETON_H */
