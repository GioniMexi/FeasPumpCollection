/**
 * @file base64.cpp
 * @brief Fast Base64 encoding/decoding Source
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#include <iostream>

#include "utils/base64.h"

namespace dominiqs {

/**
 * Base64 Conversion Tables
 */

static const char to_table[] =
{
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
	'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
	'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
	'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
	'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
	'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
	'w', 'x', 'y', 'z', '0', '1', '2', '3',
	'4', '5', '6', '7', '8', '9', '+', '/'
};

static const char from_table[] =
{
	-1, -1, -1, -1, -1, -1, -1, -1, // 0
	-1, -1, -1, -1, -1, -1, -1, -1, // 8
	-1, -1, -1, -1, -1, -1, -1, -1, // 16
	-1, -1, -1, -1, -1, -1, -1, -1, // 24
	-1, -1, -1, -1, -1, -1, -1, -1, // 32
	-1, -1, -1, 62, -1, -1, -1, 63, // 40
	52, 53, 54, 55, 56, 57, 58, 59, // 48
	60, 61, -1, -1, -1,  0, -1, -1, // 56
	-1,  0,  1,  2,  3,  4,  5,  6, // 64
	 7,  8,  9, 10, 11, 12, 13, 14, // 72
	15, 16, 17, 18, 19, 20, 21, 22, // 80
	23, 24, 25, -1, -1, -1, -1, -1, // 88
	-1, 26, 27, 28, 29, 30, 31, 32, // 96
	33, 34, 35, 36, 37, 38, 39, 40, // 104
	41, 42, 43, 44, 45, 46, 47, 48, // 112
	49, 50, 51, -1, -1, -1, -1, -1  // 120
};

ptrdiff_t b64_encode(const uint8_t* begin, const uint8_t* end, uint8_t* out, size_t out_size)
{
	ptrdiff_t input_size = end - begin;
	ptrdiff_t needed_size = 4 * (input_size / 3);
	if (input_size % 3) needed_size += 4;
	if ((ptrdiff_t)out_size < needed_size) return -needed_size;

	const uint8_t* it = begin;
	uint32_t input;

	ptrdiff_t nblocks = input_size / 3;
	// fast all blocks except the last one!
	for (ptrdiff_t i = 0; i < nblocks; i++)
	{
		input = 0;
		// read 3 bytes
		input <<= 8;
		input += static_cast<uint8_t>(*it++);
		input <<= 8;
		input += static_cast<uint8_t>(*it++);
		input <<= 8;
		input += static_cast<uint8_t>(*it++);
		// convert to base64
		*out++ = to_table[((input >> 18) & 0x3F)];
		*out++ = to_table[((input >> 12) & 0x3F)];
		*out++ = to_table[((input >> 6) & 0x3F)];
		*out++ = to_table[(input & 0x3F)];
	}
	// last block
	int bytes = 0;
	input = 0;
	// get the last bytes (and count them)
	while (it != end)
	{
		input <<= 8;
		input += static_cast<uint8_t>(*it++);
		bytes++;
	}
	// convert to base64
	int bits = bytes*8;
	while (bits > 0)
	{
		bits -= 6;
		const uint8_t index = ((bits < 0) ? input << -bits : input >> bits) & 0x3F;
		*out++ = to_table[index];
	}
	// add pad characters if necessary
	if (bytes == 1) { *out++ = '='; *out++ = '='; }
	if (bytes == 2) *out++ = '=';
	return needed_size;
}

void b64_encode(const uint8_t* begin, const uint8_t* end, std::string& out)
{
	ptrdiff_t needed = b64_encode(begin, end, 0, 0);
	out.resize((size_t)(-needed));
	b64_encode(begin, end, reinterpret_cast<uint8_t*>(&out[0]), out.size());
}

ptrdiff_t b64_decode(const uint8_t* begin, const uint8_t* end, uint8_t* out, size_t out_size)
{
	ptrdiff_t input_size = end - begin;
	ptrdiff_t needed_size = 3 * (input_size / 4);
	if (*(end - 2) == '=') needed_size -= 2;
	else if (*(end - 1) == '=') needed_size -= 1;
	if ((ptrdiff_t)out_size < needed_size) return -needed_size;

	const uint8_t* it = begin;
	uint8_t input[4];
	ptrdiff_t nblocks = input_size / 4;
	// fast all blocks except the last one!
	for (ptrdiff_t i = 0; i < nblocks - 1; i++)
	{
		// read 4 bytes
		input[0] = from_table[*it++];
		input[1] = from_table[*it++];
		input[2] = from_table[*it++];
		input[3] = from_table[*it++];
		// decode
		*out++ = static_cast<uint8_t>((input[0] << 2) + (input[1] >> 4));
		*out++ = static_cast<uint8_t>((input[1] << 4) + (input[2] >> 2));
		*out++ = static_cast<uint8_t>((input[2] << 6) + input[3]);
	}
	// last block
	input[0] = 0;
	input[1] = 0;
	input[2] = 0;
	input[3] = 0;
	// get the last bytes (and count them)
	int chars = 0;
	while(it != end)
	{
		uint8_t c = (*it++);
		if (c == '=') break; // pad character marks the end of the stream
		input[chars++] = from_table[c];
	}
	// output the binary data
	if (chars >= 2)
	{
		*out++ = static_cast<uint8_t>((input[0] << 2) + (input[1] >> 4));
		if (chars >= 3)
		{
			*out++ = static_cast<uint8_t>((input[1] << 4) + (input[2] >> 2));
			if (chars >= 4)
			{
				*out++ = static_cast<uint8_t>((input[2] << 6) + input[3]);
			}
		}
	}
	return needed_size;
}

ptrdiff_t b64_decode(const std::string& in, uint8_t* out, size_t out_size)
{
	const uint8_t* begin = reinterpret_cast<const uint8_t*>(&in[0]);
	const uint8_t* end = begin + in.size();
	return b64_decode(begin, end, out, out_size);
}

} // namespace dominiqs
