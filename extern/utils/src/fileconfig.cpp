/**
* @file fileconfig.cpp
* @brief Config files reader/writer
*
* @author Domenico Salvagnin dominiqs@gmail.com
*/

#include "utils/fileconfig.h"
#include "utils/asserter.h"
#include "utils/str_utils.h"
#include <fstream>

namespace dominiqs {


bool FileConfig::load(const std::string& filename, bool merge)
{
	if (!merge) entries.clear();

	std::ifstream in(filename);
	std::string line;
	while (std::getline(in, line))
	{
		line = trim(line);
		if (line.empty())  continue;
		if (starts_with(line, "#"))  continue;

		size_t eqIdx = line.find_first_of('=');
		DOMINIQS_ASSERT( eqIdx != std::string::npos );
		std::string entry(line.begin(), line.begin() + eqIdx);
		std::string value(line.begin() + eqIdx + 1, line.end());
		entries[entry] = value;
	}
	return true;
}


void FileConfig::save(const std::string& filename) const
{
	std::ofstream out(filename);
	for (const auto& kv : entries)
	{
		out << kv.first << " = " << kv.second << std::endl;
	}
}


std::string FileConfig::getAsString(const std::string& entry, const std::string& defValue) const
{
	auto itr = entries.find(entry);
	if (itr == entries.end())  return defValue;
	return itr->second;
}


void FileConfig::setAsString(const std::string& entry, const std::string& value)
{
	entries[entry] = value;
}


template<> std::string FileConfig::get(const std::string& entry, const std::string& defValue) const
{
	return getAsString(entry, defValue);
}


template<> void FileConfig::set(const std::string& entry, const std::string& value)
{
	setAsString(entry, value);
}


FileConfig& gConfig()
{
	static FileConfig theFile;
	return theFile;
}

} // namespace dominiqs
