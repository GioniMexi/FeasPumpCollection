/**
 * \file path.cpp
 * \brief Path management Source
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2013
 */

#include <unistd.h> //< for getcwd
#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h> //< for mkdir
#endif
#include <string.h>
#include <errno.h>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <sstream>

#include "utils/path.h"
#include "utils/str_utils.h"


namespace dominiqs {

std::string getCwd()
{
	static const size_t DEFAULT_BUFFER_SIZE = 512;
	size_t size = DEFAULT_BUFFER_SIZE;
	while (true)
	{
		std::string s;
		char* buffer = (char*)malloc(size);
		if (getcwd(buffer, size))
		{
			s = std::string(buffer);
			free(buffer);
			return s;
		}
		else
		{
			if (errno != ERANGE) return s;
			else size *= 2;
		}
		free(buffer);
	}
}

void Path::read(const char* path)
{
	data = std::string(path);
	cleanPath(data);
}

void Path::read(const std::string& path)
{
	data = path;
	cleanPath(data);
}

Path Path::operator/(const std::string& chunk)
{
	Path ret(data);
	ret /= chunk;
	return ret;
}

void Path::operator/=(const std::string& chunk)
{
	if (chunk.size())
	{
		if ((chunk[0] == '/') || (data.size() == 0)) data += chunk;
		else data += "/" + chunk;
		cleanPath(data);
	}
}

void Path::mkdir(bool recursive)
{
	if (data.empty()) return;
	if (!recursive)
	{
		int res = ::mkdir(data.c_str(), 0777);
		if (res && (errno != EEXIST)) throw std::runtime_error(strerror(errno));
		return;
	}
	else
	{
		std::string::size_type index = 1;
		while (true)
		{
			std::string::size_type found = data.find_first_of('/', index);
			std::string token = data.substr(0, found);
			int res = ::mkdir(token.c_str(), 0777);
			if (res && (errno != EEXIST)) throw std::runtime_error(strerror(errno));
			if (found == std::string::npos || found == (data.size() - 1)) return;
			index = found + 1;
		}
	}
}

bool Path::isEmpty() const
{
	return data.empty();
}

bool Path::isAbsolute() const
{
	return (!isEmpty() && data[0] == '/');
}

bool Path::isRelative() const
{
	return (!isEmpty() && data[0] != '/');
}

std::string Path::getAbsolutePath() const
{
	if (isAbsolute()) return data;
	// if path is relative, use current working dir to build absolute path
	std::string apath = getCwd() + "/" + data;
	cleanPath(apath);
	return apath;
}

std::string Path::getBasename() const
{
	if (isEmpty()) return std::string("");
	std::string::size_type start = data.find_last_of("/");
	if (start == std::string::npos) return data;
	return std::string(data.begin() + start + 1, data.end());
}

void Path::cleanPath(std::string& path) const
{
	if (isEmpty()) return;
	bool hasRoot = (path[0] == '/');
	// tokenize
	std::vector<std::string> tokens = split<std::string>(path, "/");
	// sanity check
	if (hasRoot && tokens.size() && (*tokens.begin()) == std::string("..")) throw std::runtime_error(std::string("Bad path: ") + data);
	// remove single dots
	tokens.erase(std::remove_if(tokens.begin(), tokens.end(), [](const std::string& s) { return (s == std::string(".")); }), tokens.end());
	// remove double dots
	auto itr = std::find_if(tokens.begin(), tokens.end(), [](const std::string& s) { return (s != std::string("..")); });
	auto end = tokens.end();
	while (itr != end)
	{
		if (*itr == std::string(".."))
		{
			auto tmp = itr++;
			--tmp;
			if (*tmp != std::string(".."))
			{
				tmp = tokens.erase(tmp);
				tmp = tokens.erase(tmp);
			}
		}
		else ++itr;
	}
	// rebuild string
	std::ostringstream out;
	if (hasRoot) out << "/";
	itr = tokens.begin();
	end = tokens.end();
	while (itr != end)
	{
		out << *itr++;
		if (itr != end) out << "/";
	}
	path = out.str();
}

std::string getProbName(const std::string& fileName, const std::vector<std::string>& exts)
{
	std::string probName = Path(fileName).getBasename();
	for (const std::string& ext: exts)
	{
		if (ends_with(probName, ext)) probName.resize(probName.size() - ext.size());
	}
	return probName;
}

std::string getProbName(const std::string& fileName)
{
	std::vector<std::string> exts;
	exts.push_back(".gz");
	exts.push_back(".bz2");
	exts.push_back(".mps");
	exts.push_back(".lp");
	return getProbName(fileName, exts);
}

} // namespace dominiqs
