/**
 * @file serialization.h
 * @brief Serialization routines
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2011
 */

#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include <vector>
#include <ostream>

#include "base64.h"

namespace dominiqs {

/**
 * Class used to serialize binary data to text files, using base64 encoding
 * Specialize the serialize method to support custom types
 * The class provides the buffer to write to (if needed)
 * Beware. You need this only for binary serialization!
 */

class Serializer
{
public:
	/**
	 * Serialize a binary element
	 * @tparam T element type
	 * @param value element to serialize
	 * @param out output stream
	 */
	template<class T>
	void serialize(const T& value, std::ostream& out)
	{
		b64_encode(value, buffer);
		out.write(buffer.c_str(), buffer.size());
	}
	/**
	 * Serialize a std::vector of elements
	 * @tparam T vector value type
	 * @param value vector of elements to serialize
	 * @param out output stream
	 */
	template<class T>
	void serialize(const std::vector<T>& value, std::ostream& out)
	{
		if (value.empty()) return;
		const uint8_t* begin = reinterpret_cast<const uint8_t*>(&value[0]);
		b64_encode(begin, begin + sizeof(T) * value.size(), buffer);
		out.write(buffer.c_str(), buffer.size());
	}
	/**
	 * Serialize a C array of elements
	 * @tparam T element type
	 * @param value C array to serialize
	 * @param n array size
	 * @param out output stream
	 */
	template<class T>
	void serialize(const T* value, size_t n, std::ostream& out)
	{
		if (n == 0) return;
		const uint8_t* begin = reinterpret_cast<const uint8_t*>(value);
		b64_encode(begin, begin + sizeof(T) * n, buffer);
		out.write(buffer.c_str(), buffer.size());
	}
	/**
	 * Serialize a binary element to std::string
	 * @tparam T element type
	 * @param value element to serialize
	 */
	template<class T>
	std::string toString(const T& value)
	{
		b64_encode(value, buffer);
		return buffer;
	}
	/**
	 * Serialize a std::vector of elements to std::string
	 * @tparam T vector value type
	 * @param value vector of elements to serialize
	 * @param out output stream
	 */
	template<class T>
	std::string toString(const std::vector<T>& value)
	{
		if (value.empty()) return "";
		const uint8_t* begin = reinterpret_cast<const uint8_t*>(&value[0]);
		b64_encode(begin, begin + sizeof(T) * value.size(), buffer);
		return buffer;
	}
	/**
	 * Serialize a C array of elements to std::string
	 * @tparam T element type
	 * @param value C array to serialize
	 * @param n array size
	 * @param out output stream
	 */
	template<class T>
	std::string toString(const T* value, size_t n)
	{
		if (n == 0) return "";
		const uint8_t* begin = reinterpret_cast<const uint8_t*>(value);
		b64_encode(begin, begin + sizeof(T) * n, buffer);
		return buffer;
	}
protected:
	std::string buffer;
};

} // namespace dominiqs

#endif /* SERIALIZATION_H */
