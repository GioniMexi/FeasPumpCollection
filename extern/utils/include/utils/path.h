/**
 * @file path.h
 * @brief Path management Header
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2011
 */

#ifndef PATH_H
#define PATH_H

#include <string>

namespace dominiqs {

/**
 * @brief Get The Current Working Directory
 * Convenience wrapper around getcwd
 */

std::string getCwd();

/**
 * @brief Path Class
 * Convenience path handling class
 * Its main usage is to recursively create directories.
 * Eventually there will be a portable way to do it...
 */

class Path
{
public:
	Path(const char* path) { read(path); }
	Path(const std::string& path) { read(path); }
	/**
	 * Read path from C-string or std::string
	 */
	void read(const char* path);
	void read(const std::string& path);
	/**
	 * Concatenate current path with chunk @param chunk
	 * @return new path
	 */
	Path operator/(const std::string& chunk);
	/**
	 * Concatenate current path with chunk @param chunk inplace
	 */
	void operator/=(const std::string& chunk);
	/**
	 * Make directory from path (optionally recursive)
	 */
	void mkdir(bool recursive = true);
	/**
	 * Test methods
	 */
	bool isEmpty() const;
	bool isAbsolute() const;
	bool isRelative() const;
	std::string getPath() const { return data; }
	/**
	 * Get absolute path corresponding to current path
	 */
	std::string getAbsolutePath() const;
	/**
	 * Get basename (including extensions if present) of current path
	 */
	std::string getBasename() const;
protected:
	// normalize path
	void cleanPath(std::string& path) const;
	// path string
	std::string data;
};

/**
 * @brief Get the problem name from a LP/MPS file
 */

std::string getProbName(const std::string& fileName, const std::vector<std::string>& exts);
std::string getProbName(const std::string& fileName);

} // namespace dominiqs

#endif /* PATH_H */
