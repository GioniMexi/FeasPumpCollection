/**
 * @file str_utils.h
 * @brief String Utilities
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2011-2013
 */

#ifndef STR_UTILS_H
#define STR_UTILS_H

#include <string>
#include <iterator>
#include <cctype>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>


namespace dominiqs {

/**
 * Convert from string to an arbitrary type
 */

template<typename T>
inline T from_string(const std::string& str)
{
	// default (slow implementation)
	T ret;
	std::istringstream iss(str);
	iss >> ret;
	return ret;
}

/* Template specializations */
template<>
inline std::string from_string(const std::string& str) { return str; }

template<>
inline int from_string(const std::string& str) { return std::stoi(str); }

template<>
inline long from_string(const std::string& str) { return std::stol(str); }

template<>
inline long long from_string(const std::string& str) { return std::stoll(str); }

template<>
inline unsigned long from_string(const std::string& str) { return std::stoul(str); }

template<>
inline unsigned long long from_string(const std::string& str) { return std::stoull(str); }

template<>
inline float from_string(const std::string& str) { return std::stof(str); }

template<>
inline double from_string(const std::string& str) { return std::stod(str); }

template<>
inline long double from_string(const std::string& str) { return std::stold(str); }

/**
 * Convert a numerical type to string with a specified precision
 */
template<typename T>
inline std::string to_string(T num, const int p)
{
	std::ostringstream out;
	out << std::setprecision(p) << num;
	return out.str();
}

/**
 * Return true if string str starts with prefix, false otherwise
 */

inline bool starts_with(const std::string& str, const std::string& prefix)
{
	return ((prefix.size() <= str.size()) && std::equal(prefix.begin(), prefix.end(), str.begin()));
}

/**
 * Return true if string str ends with suffix, false otherwise
 */

inline bool ends_with(const std::string& str, const std::string& suffix)
{
	return ((suffix.size() <= str.size()) && std::equal(suffix.rbegin(), suffix.rend(), str.rbegin()));
}

/**
 * Trim whitespace from a string (*not* inplace)
 */
inline std::string trim(const std::string &s)
{
	auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c){return std::isspace(c);});
	auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c){return std::isspace(c);}).base();
	return ((wsback <= wsfront) ? std::string() : std::string(wsfront, wsback));
}

/**
 * Escape character data, replacing &<>"' with the corresponding
 * XML entities
 */

template<typename InputIterator>
void xmlEscape(InputIterator first, InputIterator end, std::ostream& out)
{
	while (first != end)
	{
		char ch = *first++;
		if (ch == '&') out << "&amp;";
		else if (ch == '<') out << "&lt;";
		else if (ch == '>') out << "&gt;";
		else if (ch == '\'') out << "&apos;";
		else if (ch == '"') out << "&quot;";
		else out << ch;
	}
}

/**
 * Print Python-like unnested list and maps
 */

namespace implementation
{

/**
 * Type-dependent quote character
 */

template<typename T>
inline
const char* quote() { return ""; }

template<>
inline
const char* quote<const std::string>() { return "\'"; }

template<>
inline
const char* quote<std::string>() { return "\'"; }

template<>
inline
const char* quote<const char*>() { return "\'"; }

template<>
inline
const char* quote<char*>() { return "\'"; }

} // namespace implementation

/**
 * Print range [first,last) of a container as a Python-like unnested list
 * @param out output stream
 * @param first beginnning of the range
 * @param last end of the range
 * @param precision precision for numeric values
 * @param nl true to add a newline at the end, false (default) otherwise
 */

template<typename BidirectionalIter>
void printList(std::ostream& out,
	BidirectionalIter first, BidirectionalIter last,
	std::streamsize precision = 15, bool nl = false)
{
	if (first != last)
	{
		using namespace implementation;
		typedef typename std::iterator_traits<BidirectionalIter>::value_type VT;
		// set precision
		std::streamsize p = out.precision();
		out.precision(precision);
		// print list
		out << "[";
		BidirectionalIter beforeLast = last;
		beforeLast--;
		while (first != beforeLast)
		{
			out << quote<VT>() << *first << quote<VT>() << ", ";
			first++;
		}
		out << quote<VT>() << *first << quote<VT>() << "]";
		// restore precision
		out.precision(p);
	}
	else out << "[]";
	if (nl) out << std::endl;
}

/**
 * Convert range [first,last) of a container to a string with a Python-like syntax
 * @param first beginnning of the range
 * @param last end of the range
 * @param precision precision for numeric values
 * @return string encoding of the range
 */

template<typename BidirectionalIter>
std::string list2str(BidirectionalIter first, BidirectionalIter last,
	std::streamsize precision = 15)
{
	std::ostringstream os;
	printList(os, first, last, precision, false);
	return os.str();
}

/**
 * Print possibly long vector with name
 * @param out output stream
 * @param container random access container
 * @param precision precision for numeric values
 * @param nl true to add a newline at the end, false (default) otherwise
 */

template<typename T>
void printLongVector(std::ostream& out, const std::string& name, const std::vector<T>& container,
	std::streamsize precision = 3, size_t howmany = 5, bool nl = false)
{
	if (container.empty())
	{
		out << name << ": []";
		if (nl) out << std::endl;
		return;
	}
	if (container.size() <= 3*howmany)
	{
		out << name << ": ";
		printList(std::cout, container.begin(), container.end(), precision, nl);
		return;
	}
	// set precision
	std::streamsize p = out.precision();
	out.precision(precision);
	out << name << ": [";
	for (size_t i = 0; i < howmany; i++) out << container[i] << ", ";
	out << " ... ";
	size_t n = container.size();
	for (size_t i = n - howmany; i < n - 1; i++) out << container[i] << ", ";
	out << container[n-1] << "]";
	if (nl) out << std::endl;
	// restore precision
	std::cout.precision(p);
}


/**
 * Print range [first,last) of an associative container as a Python-like map
 * @param out output stream
 * @param first beginnning of the range
 * @param last end of the range
 * @param precision precision for numeric values
 * @param nl true to add a newline at the end, false (default) otherwise
 */

template<typename BidirectionalIter>
void printMap(std::ostream& out,
	BidirectionalIter first, BidirectionalIter last,
	std::streamsize precision = 15, bool nl = false)
{
	if (first != last)
	{
		using namespace implementation;
		typedef typename std::iterator_traits<BidirectionalIter>::value_type PairType;
		typedef typename PairType::first_type KT;
		typedef typename PairType::second_type VT;
		// set precision
		std::streamsize p = out.precision();
		out.precision(precision);
		// print map
		out << "{";
		BidirectionalIter beforeLast = last;
		beforeLast--;
		while (first != beforeLast)
		{
			out << quote<KT>() << first->first << quote<KT>() << ": " << quote<VT>() << first->second << quote<VT>() << ", ";
			++first;
		}
		out << quote<KT>() << first->first << quote<KT>() << ": " << quote<VT>() << first->second << quote<VT>() << "}";
		// restore precision
		out.precision(p);
	}
	else out << "{}";
	if (nl) out << std::endl;
}

/**
 * Print vector as a matrix of size m x n
 * @param out output stream
 * @param name name of the container
 * @param first iterator to beginning of sequence holding the matrix
 * @param m number of rows
 * @param n number of columns
 * @param precision precision for numeric values
 */

template<typename BidirectionalIter>
void printMatrix(std::ostream& out, const std::string& name, BidirectionalIter first,
	int m, int n, std::streamsize width = 8, std::streamsize precision = 3)
{
	// set precision/width
	out << name << ":" << std::endl;
	// print matrix
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			out << std::setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(width) << *first++;
		}
		out << std::endl;
	}
}

/**
 * Python-like string splitting + type conversion
 * @param str string to split
 * @param delim string of delimiters
 * @return vector of values of type T
 */

template<typename T>
std::vector<T> split(const std::string& str, const std::string& delim)
{
	size_t start, end = 0;
	std::vector<T> ret;
	while (end < str.size())
	{
		start = end;
		// skip initial whitespace
		while (start < str.size() && (delim.find(str[start]) != std::string::npos)) start++;
		// skip to end of word
		end = start;
		while (end < str.size() && (delim.find(str[end]) == std::string::npos)) end++;
		// just ignore zero-length strings.
		if (end-start != 0) ret.push_back(from_string<T>(std::string(str, start, end-start)));
	}
	return ret;
}

} // namespace dominiqs

#endif /* STR_UTILS_H */
