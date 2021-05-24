/**
 * @file base64.h
 * @brief Fast Base64 encoding/decoding
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#ifndef BASE64_H
#define BASE64_H

#include <stdint.h> //< for uint8_t
#include <stddef.h> //< for ptrdiff_t and size_t
#include <type_traits>
#include <string>
#include <stdexcept>

namespace dominiqs {

/**
 * Encode a binary range [begin, end) into a buffer out of size out_size
 * @return the number of characters written (if successful)
 * If out_size is to small return -need_size
 */

ptrdiff_t b64_encode(const uint8_t* begin, const uint8_t* end, uint8_t* out, size_t out_size);

/**
 * Convenience functions accepting std::string for output and custom types in input (beware: use only PODs)
 */

void b64_encode(const uint8_t* begin, const uint8_t* end, std::string& out);

template<class T>
void b64_encode(const T& obj, std::string& out)
{
	static_assert(std::is_pod<T>::value, "Cannot Base64 encode non-POD");
	const uint8_t* begin = reinterpret_cast<const uint8_t*>(&obj);
	const uint8_t* end = begin + sizeof(T);
	b64_encode(begin, end, out);
}

/**
 * Decode a binary range [begin, end) into buffer out of size out_size
 * @return the number of bytes written (if successful)
 * If out_size is to small return -need_size
 * There cannot be bad characters in the stream (e.g. newlines, spaces, etc...): no checks!
 */

ptrdiff_t b64_decode(const uint8_t* begin, const uint8_t* end, uint8_t* out, size_t out_size);

/**
 * Convenience functions accepting std::string for input and custom types in output (beware: use only PODs)
 */

ptrdiff_t b64_decode(const std::string& in, uint8_t* out, size_t out_size);

template<class T>
ptrdiff_t b64_decode(const std::string& in, T& out)
{
	static_assert(std::is_pod<T>::value, "Cannot Base64 decode non-POD");
	ptrdiff_t ret = b64_decode(in, reinterpret_cast<uint8_t*>(&out), sizeof(T));
	if (ret < 0) throw std::runtime_error("Bad base64 decoding: maybe wrong type?");
	return ret;
}

} // namespace dominiqs

#endif /* BASE64_H */
