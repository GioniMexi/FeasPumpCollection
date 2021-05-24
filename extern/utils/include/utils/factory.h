/**
 * @file factory.h
 * @brief Factory Pattern
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2013
 */

#ifndef FACTORY_H
#define FACTORY_H

#include <map>
#include <utility>

namespace dominiqs {

template<typename BaseClassType, typename ClassType, typename... Args>
BaseClassType* Creator(Args... args)
{
	return new ClassType(args...);
}

template<typename BaseClassType, typename UniqueIdType, typename... Args>
class Factory
{
private:
	using CreatorType = BaseClassType* (*)(Args...);
	using CMap = ::std::map<UniqueIdType, CreatorType>;
public:
	/**
	 * Register class ClassType with name id
	 */
	template<typename ClassType>
	bool registerClass(UniqueIdType id)
	{
		if (m.find(id) != m.end()) return false;
		m[id] = &Creator<BaseClassType, ClassType, Args...>;
		return true;
	}
	/**
	 * Unregister class with name id
	 */
	bool unregisterClass(UniqueIdType id)
	{
		return (m.erase(id) == 1);
	}
	/**
	 * Instantiate a class of name id
	 */
	BaseClassType* create(UniqueIdType id, Args... args)
	{
		auto iter = m.find(id);
		if (iter == m.end()) return 0;
		return ((*iter).second)(std::forward<Args>(args)...);
	}
	/**
	 * Returns the lists of currently registered names (ids)
	 */
	template<typename OutputIterator>
	void getIDs(OutputIterator out) const
	{
		for (auto itr: m) *out++ = itr.first;
	}
protected:
	CMap m;
};

} // namespace dominiqs

#endif // FACTORY_H

